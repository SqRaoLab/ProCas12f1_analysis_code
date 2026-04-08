package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"os/exec"
	"sync"
	"time"

	"github.com/k0kubun/go-ansi"

	"github.com/schollz/progressbar/v3"
)

// --- 文件读取和并发处理 ---
// progressReader 是一个包装器，用于在读取时更新进度条
type progressReader struct {
	reader io.Reader // 原始 Reader (例如 *os.File)
	bar    *progressbar.ProgressBar
}

// Read 实现 io.Reader 接口
func (pr *progressReader) Read(p []byte) (n int, err error) {
	// 调用原始 Reader 的 Read 方法
	n, err = pr.reader.Read(p)
	// 更新进度条
	if n > 0 {
		pr.bar.Add(n)
	}
	return n, err
}

func progressBar(size int64) *progressbar.ProgressBar {
	if size > 0 {
		return progressbar.NewOptions64(
			size,
			progressbar.OptionSetDescription("Reading"),       // 设置描述
			progressbar.OptionSetWriter(ansi.NewAnsiStderr()), // 使用ansi防止奇怪的换行等显示错误
			progressbar.OptionUseANSICodes(true),              // avoid progressbar downsize error
			progressbar.OptionEnableColorCodes(true),
			progressbar.OptionSetElapsedTime(true),
			progressbar.OptionSetPredictTime(true),
			progressbar.OptionSetWriter(os.Stderr), // 避免 stdout 冲突
			progressbar.OptionShowBytes(true),      // 显示已读/总字节数
			progressbar.OptionShowCount(),          // 显示计数（可选）
			progressbar.OptionFullWidth(),          // 全宽显示（可选）
			progressbar.OptionThrottle(1*time.Second),
		)
	}
	return progressbar.NewOptions(-1,
		progressbar.OptionSetDescription("Reading"),
		progressbar.OptionSetWriter(ansi.NewAnsiStderr()), // 使用ansi防止奇怪的换行等显示错误
		progressbar.OptionUseANSICodes(true),              // avoid progressbar downsize error
		progressbar.OptionEnableColorCodes(true),
		progressbar.OptionSetItsString("items"),
		progressbar.OptionSpinnerType(14),      // 旋转光标
		progressbar.OptionSetWriter(os.Stderr), // 避免 stdout 冲突
		progressbar.OptionSetItsString("files"),
		progressbar.OptionThrottle(1*time.Second),
	)
}

func readFromCmd(filename, cmd_ string, recordChan chan<- []string, wg *sync.WaitGroup, showPprogress bool) {
	defer wg.Done()
	defer close(recordChan) // 读取完成后关闭channel

	cmd := exec.Command("bash", "-c", fmt.Sprintf("%s %s", cmd_, filename))

	// 获取命令的标准输出管道
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		sugar.Fatalf("Error creating StdoutPipe: %v", err)
	}

	// 开始执行命令
	if err := cmd.Start(); err != nil {
		sugar.Fatalf("Error starting command: %v", err)
	}

	defer stdout.Close()

	scanner := bufio.NewScanner(stdout)
	buffer := make([]string, 0, 4) // 缓冲4行

	var bar *progressbar.ProgressBar
	if showPprogress {
		bar = progressBar(-1)
	}

	for scanner.Scan() {
		line := scanner.Text()
		buffer = append(buffer, line)

		if bar != nil {
			bar.Add(1)
		}

		if len(buffer) == 4 { // 一个完整的FASTQ记录
			// 发送记录的副本到channel
			recordCopy := make([]string, 4)
			copy(recordCopy, buffer)
			recordChan <- recordCopy
			buffer = buffer[:0] // 重置缓冲区
		}
	}
	if err := scanner.Err(); err != nil {
		sugar.Fatalf("Error reading from %s: %v", filename, err)
	}
	if bar != nil {
		bar.Finish()
	}

	// 处理文件末尾可能存在的不完整记录（理论上不应存在）
	if len(buffer) > 0 {
		sugar.Infof("Warning: Incomplete FASTQ record found at end of %s", filename)
	}
}

// readFastqInChannel 从gzip压缩的FASTQ文件读取记录，并发送到channel
func readFastqInChannel(filename string, recordChan chan<- []string, wg *sync.WaitGroup, showProgress bool) {
	defer wg.Done()
	defer close(recordChan) // 读取完成后关闭channel

	file, err := os.Open(filename)
	if err != nil {
		sugar.Fatalf("Error opening file %s: %v", filename, err)
	}

	var gzReader *gzip.Reader
	var bar *progressbar.ProgressBar
	if showProgress {
		if stat, err := os.Stat(filename); !os.IsNotExist(err) {
			bar = progressBar(stat.Size())

			// --- 4. 创建进度条包装器 ---
			progressReaderInstance := &progressReader{
				reader: file, // 包装原始文件 Reader
				bar:    bar,  // 关联进度条
			}
			gzReader, err = gzip.NewReader(progressReaderInstance)
			if err != nil {
				sugar.Fatalf("Error creating gzip reader for %s: %v", filename, err)
			}
		} else {
			sugar.Fatal(err)
		}
	} else {
		gzReader, err = gzip.NewReader(file)
		if err != nil {
			sugar.Fatalf("Error creating gzip reader for %s: %v", filename, err)
		}
	}

	scanner := bufio.NewScanner(gzReader)
	buffer := make([]string, 0, 4) // 缓冲4行

	for scanner.Scan() {
		line := scanner.Text()

		if line == "" {
			continue
		}

		buffer = append(buffer, line)

		if len(buffer) == 4 { // 一个完整的FASTQ记录
			// 发送记录的副本到channel
			recordCopy := make([]string, 4)
			copy(recordCopy, buffer)
			recordChan <- recordCopy
			buffer = buffer[:0] // 重置缓冲区
		}
	}

	if err := scanner.Err(); err != nil {
		sugar.Fatalf("Error reading from %s: %v", filename, err)
	}

	if bar != nil {
		_ = bar.Finish()
	}

	if gzReader != nil {
		_ = gzReader.Close()
	}
	_ = file.Close()

	// 处理文件末尾可能存在的不完整记录（理论上不应存在）
	if len(buffer) > 0 {
		sugar.Infof("Warning: Incomplete FASTQ record found at end of %s", filename)
	}
}
