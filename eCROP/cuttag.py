#!/usr/bin/env python3
"""
CUTtag 分析流程 - 智能样本发现版
自动从文件名中提取样本信息，根据关键词分类：
  - MFLAG/IgG → Control
  - RFLAG / MFLAG-1 / MFLAG-2 → ChIP (根据Target和Batch自动配对Control)
  - PC → Positive Control
  - NA → 跳过
"""

import glob
import os
import re
import subprocess
import sys
from pathlib import Path

# ==================== 路径配置 ====================
# 请根据实际环境修改以下路径
BASE_DIR = os.environ.get("PROJECT_DIR", "/path/to/project")
RAW_DIR = f"{BASE_DIR}/rawdata/CUTtag/Results/Lane00"
REF_FA = "/path/to/reference/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_INDEX = "/path/to/reference/hg38/hg38_bt2"
RESULT_DIR = f"{BASE_DIR}/res/cuttag"
LOG_DIR = f"{RESULT_DIR}/logs"

N_CORES = 8

# CUTtag文件名正则: ...CUTtag_<SAMPLE>_L00_R[12].fq.gz
RE_FILE = re.compile(r"CUTtag_(.+?)_(L00)_R([12])\.fq\.gz")

# 目标蛋白关键词
TARGETS = ["SCMH1", "SMYD3", "ZNF555"]


# ==================== 工具函数 ====================

def check_bowtie2_index():
    """检查并构建Bowtie2索引 (.bt2 为 bowtie2>=2.5, .ebwt 为旧版)"""
    for ext in ["bt2", "ebwt"]:
        indices = [f"{REF_INDEX}.{i}.{ext}" for i in ["1", "2", "3", "4", "rev.1", "rev.2"]]
        missing = [i for i in indices if not os.path.exists(i)]
        if not missing:
            print("Bowtie2索引已存在，跳过构建")
            return
    # 所有格式都不完整，需要构建
    print("Bowtie2索引不存在，正在构建...")
    subprocess.run(["bowtie2-build", "--threads", str(N_CORES), REF_FA, REF_INDEX], check=True)
    print("Bowtie2索引构建完成")


def find_file(sample_name):
    """从全局样本映射中查找R1/R2文件"""
    info = SAMPLE_FILES.get(sample_name, {})
    return info.get("R1"), info.get("R2")


def run_cmd(cmd, label=""):
    """运行命令并打印日志"""
    print(f"[{label}] 执行: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 and result.stderr:
        print(f"[{label}] STDERR: {result.stderr[:500]}")
    return result


def run_fastqc(sample_name):
    r1, r2 = find_file(sample_name)
    if not r1 or not r2:
        print(f"[WARNING] {sample_name}: R1或R2文件不存在，跳过")
        return
    out_dir = f"{RESULT_DIR}/qc/{sample_name}"
    os.makedirs(out_dir, exist_ok=True)
    log_file = f"{LOG_DIR}/fastqc_{sample_name}.log"
    cmd = ["fastqc", "-t", str(N_CORES), "-o", out_dir, "--nogroup", r1, r2]
    with open(log_file, "a") as log:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log)
        proc.wait()
    print(f"[QC] FastQC完成: {sample_name}")


def run_bowtie2(sample_name):
    r1, r2 = find_file(sample_name)
    if not r1 or not r2:
        print(f"[WARNING] {sample_name}: R1或R2文件不存在，跳过")
        return

    out_dir = f"{RESULT_DIR}/aligned/{sample_name}"
    os.makedirs(out_dir, exist_ok=True)
    bam_file = f"{out_dir}/{sample_name}.raw.bam"
    sam_file = f"{out_dir}/{sample_name}.sam"
    log_file = f"{LOG_DIR}/bowtie2_{sample_name}.log"

    with open(log_file, "a") as log:
        result = subprocess.run(
            ["bowtie2", "-x", REF_INDEX, "-1", r1, "-2", r2, "-p", str(N_CORES),
             "--very-sensitive", "-X", "2000", "-S", sam_file],
            stdout=log, stderr=log
        )
        if result.returncode != 0:
            print(f"[ALIGN] {sample_name}: bowtie2返回码{result.returncode}")
            return

    # SAM -> BAM -> 排序 -> 索引
    run_cmd(["samtools", "view", "-@", str(N_CORES), "-bS", sam_file, "-o", bam_file], "SAM2BAM")
    run_cmd(["samtools", "index", "-@", str(N_CORES), bam_file], "INDEX")
    os.remove(sam_file)

    # 统计比对率
    total = int(subprocess.check_output(
        ["samtools", "view", "-c", "-@", str(N_CORES), bam_file], text=True).strip())
    mapped = int(subprocess.check_output(
        ["samtools", "view", "-c", "-F", "4", "-@", str(N_CORES), bam_file], text=True).strip())
    rate = f"{mapped/total*100:.2f}%" if total else "0%"
    print(f"[ALIGN] 比对完成: {sample_name} | 总数: {total} | 比对: {mapped} | 率: {rate}")
    return total, mapped, rate


def run_dedup(sample_name):
    bam_file = f"{RESULT_DIR}/aligned/{sample_name}/{sample_name}.raw.bam"
    if not os.path.exists(bam_file):
        print(f"[WARNING] {sample_name}: BAM文件不存在，跳过")
        return

    out_dir = f"{RESULT_DIR}/filtered/{sample_name}"
    os.makedirs(out_dir, exist_ok=True)
    dedup_bam = f"{out_dir}/{sample_name}.dedup.bam"

    # 过滤: 正确配对、MAPQ>=30
    filt_bam = f"{out_dir}/{sample_name}.filtered.bam"
    run_cmd(["samtools", "view", "-@", str(N_CORES), "-b", "-f", "2", "-F", "772", "-q", "30",
             bam_file, "-o", filt_bam], "FILTER")

    # 去重
    try:
        run_cmd(["samtools", "markdup", "-@", str(N_CORES), filt_bam, dedup_bam], "MARKDUP")
    except subprocess.CalledProcessError:
        print("[DEDUP] markdup失败，尝试rmdup...")
        run_cmd(["samtools", "rmdup", "-@", str(N_CORES), filt_bam, dedup_bam], "RMDUP")

    run_cmd(["samtools", "index", "-@", str(N_CORES), dedup_bam], "INDEX")

    # 统计
    filtered = int(subprocess.check_output(
        ["samtools", "view", "-c", "-@", str(N_CORES), filt_bam], text=True).strip())
    dedup = int(subprocess.check_output(
        ["samtools", "view", "-c", "-@", str(N_CORES), dedup_bam], text=True).strip())
    dup_rate = f"{(1-dedup/filtered)*100:.2f}%" if filtered else "0%"
    print(f"[DEDUP] 去重完成: {sample_name} | 过滤后: {filtered} | 去重后: {dedup} | 重复率: {dup_rate}")


def run_macs3(target, sample, control):
    sample_bam = f"{RESULT_DIR}/filtered/{sample}/{sample}.dedup.bam"
    control_bam = f"{RESULT_DIR}/filtered/{control}/{control}.dedup.bam"

    if not os.path.exists(sample_bam):
        print(f"[PEAK] {sample}: BAM文件不存在，跳过")
        return
    if not os.path.exists(control_bam):
        print(f"[PEAK] {sample}: Control({control})不存在，跳过")
        return

    out_dir = f"{RESULT_DIR}/peaks/{target}/{sample}"
    os.makedirs(out_dir, exist_ok=True)
    log_file = f"{LOG_DIR}/macs3_{sample}.log"

    with open(log_file, "a") as log:
        result = subprocess.run(
            ["macs3", "callpeak",
             "-t", sample_bam, "-c", control_bam,
             "-f", "BAM", "-g", "hs",
             "-n", f"{target}_{sample}",
             "-q", "0.01", "--extsize", "200",
             "-o", out_dir],
            stdout=log, stderr=log
        )
    print(f"[PEAK] Peak calling完成: {sample} vs {control}")


# ==================== 全局样本映射 ====================
SAMPLE_FILES = {}  # sample_name -> {"R1": path, "R2": path}

def discover_samples():
    """从目录自动发现所有样本文件，返回样本名到文件路径的映射"""
    samples = {}
    fq_files = glob.glob(f"{RAW_DIR}/*_R1.fq.gz")
    fq_files += glob.glob(f"{RAW_DIR}/*_R2.fq.gz")

    for f in fq_files:
        m = RE_FILE.search(f)
        if not m:
            continue
        sample = m.group(1)
        pair = m.group(3)  # "1" or "2"
        if sample not in samples:
            samples[sample] = {"R1": None, "R2": None}
        if pair == "1":
            samples[sample]["R1"] = f
        else:
            samples[sample]["R2"] = f

    # 过滤掉只有单端的样本
    valid = {k: v for k, v in samples.items() if v["R1"] and v["R2"]}
    return valid


def classify_samples(samples):
    """根据文件名关键词分类样本"""
    chip = []       # ChIP样本
    controls = []   # Control样本 (IgG)
    qc_special = [] # 特殊QC (PC)
    skipped = []    # 跳过 (NA等)

    for name in sorted(samples.keys()):
        if "NA" in name:
            skipped.append(name)
        elif "PC" in name:
            qc_special.append(name)
        elif any(t in name for t in TARGETS):
            # 有Target关键词 -> 提取Target
            target = None
            for t in TARGETS:
                if t in name:
                    target = t
                    break

            # 所有MFLAG和RFLAG都是ChIP样本
            if "RFLAG" in name or "MFLAG" in name:
                chip.append((target, name))
            else:
                chip.append((target, name))
        elif "IgG" in name:
            # IgG是Control
            controls.append(name)
        else:
            print(f"[INFO] 未识别的样本: {name}")
            skipped.append(name)

    return chip, controls, qc_special, skipped


def get_control(chip_target, chip_name):
    """根据ChIP样本名称自动找到对应的Control"""
    # 根据样本前缀匹配对应的Control
    # 请根据实际样本命名规则修改此处的前缀
    if chip_name.startswith("SampleA"):
        control = "SampleA-IgG"
    elif chip_name.startswith("SampleB"):
        control = "SampleB-IgG"
    else:
        control = None
    return control


# ==================== 主流程 ====================

def main():
    print("=" * 60)
    print("  CUTtag 分析流程启动")
    print(f"  {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    # 创建目录
    os.makedirs(RESULT_DIR, exist_ok=True)
    for d in ["qc", "aligned", "filtered", "peaks", "logs"]:
        os.makedirs(f"{RESULT_DIR}/{d}", exist_ok=True)

    # 检查Bowtie2索引
    check_bowtie2_index()

    # 发现样本
    print("\n========== 样本发现 ==========")
    samples = discover_samples()
    print(f"发现 {len(samples)} 个样本")

    # 填充全局样本映射
    global SAMPLE_FILES
    SAMPLE_FILES = samples

    # 分类
    chip, controls, qc_special, skipped = classify_samples(samples)

    print(f"  ChIP样本: {len(chip)}")
    for t, n in chip:
        ctrl = get_control(t, n)
        print(f"    {n} (target={t}, control={ctrl})")
    print(f"  Control样本: {len(controls)} -> {controls}")
    print(f"  特殊QC: {qc_special}")
    print(f"  跳过: {skipped}")

    # 去重Control
    controls = list(set(controls))

    # 初始化对齐统计
    with open(f"{RESULT_DIR}/alignment_summary.csv", "w") as f:
        f.write("sample,total,mapped,rate\n")

    # ============ Step 1: QC ============
    print("\n========== Step 1: FastQC Quality Control ==========")

    # 所有需要QC的样本
    all_qc = [n for _, n in chip] + controls + qc_special
    for sample in all_qc:
        run_fastqc(sample)
    print("[QC] 全部样本FastQC完成")

    # ============ Step 2: Bowtie2 ============
    print("\n========== Step 2: Bowtie2 Alignment ==========")

    # 先比对Control
    for control in controls:
        run_bowtie2(control)

    # 再比对ChIP
    for _, sample in chip:
        run_bowtie2(sample)
    print("[ALIGN] 全部样本比对完成")

    # ============ Step 3: 去重 ============
    print("\n========== Step 3: Filtering & Deduplication ==========")
    for control in controls:
        run_dedup(control)
    for _, sample in chip:
        run_dedup(sample)
    print("[DEDUP] 全部样本去重完成")

    # ============ Step 4: MACS3 Peak Calling ============
    print("\n========== Step 4: MACS3 Peak Calling ==========")
    for target, sample in chip:
        control = get_control(target, sample)
        run_macs3(target, sample, control)
    print("[PEAK] 全部样本Peak calling完成")

    # ============ 汇总 ============
    print("\n========== 结果汇总 ==========")
    print(f"结果目录: {RESULT_DIR}/")
    print("  ├── qc/              (FastQC)")
    print("  ├── aligned/         (BAM)")
    print("  ├── filtered/        (去重BAM)")
    print("  ├── peaks/           (Peak calling)")
    print("  ├── alignment_summary.csv")
    print("  └── logs/            (日志)")

    # MultiQC
    print("\n========== MultiQC汇总 ==========")
    if os.path.exists("/usr/bin/multiqc"):
        print("[MultiQC] 生成汇总报告...")
        run_cmd(["multiqc", RESULT_DIR + "/qc", "-o", RESULT_DIR + "/qc", "--force"], "MULTIQC")
        if os.path.exists(f"{RESULT_DIR}/qc/multiqc_report.html"):
            print(f"[MultiQC] 报告: {RESULT_DIR}/qc/multiqc_report.html")
    else:
        print("[MultiQC] 未找到multiqc，跳过")

    # Peak统计
    print("\n========== Peak统计 ==========")
    print(f"{'Target':<12} {'Sample':<35} {'Peaks'}")
    print("-" * 65)

    peaks_dir = f"{RESULT_DIR}/peaks"
    if os.path.exists(peaks_dir):
        for target_dir in sorted(Path(peaks_dir).iterdir()):
            if not target_dir.is_dir():
                continue
            target = target_dir.name
            for sample_dir in sorted(target_dir.iterdir()):
                if not sample_dir.is_dir():
                    continue
                peak_file = sample_dir / f"{target}_{sample_dir.name}_peaks.narrowPeak"
                if not peak_file.exists():
                    # 尝试其他可能的命名
                    peak_file = list(sample_dir.glob("*_peaks.narrowPeak"))
                    if not peak_file:
                        continue
                    peak_file = peak_file[0]
                peak_count = sum(1 for line in open(peak_file) if not line.startswith("track"))
                print(f"{target:<12} {sample_dir.name:<35} {peak_count}")

    print()
    print("=" * 60)
    print("  CUTtag分析流程完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
