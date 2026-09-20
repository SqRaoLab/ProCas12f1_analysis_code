# ProCas12f1 analysis code

This repository contains the analysis code for ProCas12f1.

The `editing_analysis`, `eCROP` folders store the code for sequencing read splitting, editing efficiency calculation, and eCROP data analysis, respectively.

---

## 1. Editing Analysis

A comprehensive pipeline for analyzing CRISPR-Cas gene editing experiments from high-throughput sequencing data. This toolkit processes raw sequencing reads, identifies editing events, and quantifies editing efficiencies. For detailed usage, refer to the [splitFq README](./editing_analysis/splitFq/README.md) and [decodeFq README](./editing_analysis/decodeFq/README.md).

### 1.1 Installation

[Go version >= 1.20](https://go.dev/) is required to compile these tools. Please install Golang following the official instructions.

### 1.2 Argument Details

#### (1) splitFq

The `editing_analysis/splitFq` directory contains the tool used to split raw sequencing data by barcodes.

```bash
cd editing_analysis/splitFq

# install requirements
go mod tidy

# compile this tool
go build .

# print help info
./splitFq --help
Required parameters:

bash
Usage: splitFq [global options] 

Global options:
        -b, --barcode the path to json file with barcodes
        -i, --input   the path to json with file path
        -o, --output  the output directory (default: output.txt)
            --or      using or instead of and
            --debug   enable debug log
        -v, --version print version
        -h, --help    Show this help
```

#### (2) decodeFq

The editing_analysis/decodeFq directory contains the tool used to decode PAM, spacer, target, and editing events from FASTQ files.

```bash
cd editing_analysis/decodeFq

# install requirements
go mod tidy

# compile this tool
go build .

# print help info
./decodeFq --help
Required parameters:

bash
Usage: decodeFq [global options] 

Global options:
        -1, --r1              R1 FASTQ (gzipped)
        -2, --r2              R2 FASTQ (gzipped)
        -c, --cas             the Cas protein
        -l, --library         the path to json with pool library information
        -o, --output          the output file path (default: output.txt)
        -p, --process         the number of goroutines to use (default: 10)
            --target-length   length of gDNA target (default: 20)
            --primer-length   length of reverse primer (default: 20)
            --pam-length      length of PAM (default: 4)
            --spacer-length   length of spacer (default: 20)
            --spacer-distance length between spacer and before target (default: 6)
            --before-length   length of before target (default: 6)
            --behind-length   length of before target (default: 6)
            --cmd             used to read fastq file from quip format
        -U, --umi-anchor      the designed umi
        -u, --umi             the length of designed umi (default: 28)
        -r, --reverse         complement reverse reads
            --debug           enable debug
        -v, --version         print version
        -h, --help            Show this help
```

### 1.3 Example Pipeline

```bash
# split raw data by barcodes
./editing_analysis/splitFq/splitFq -b example/barcodes.json -i example/files.json -o example

# decode editing events from genome and plasmid samples
./editing_analysis/decodeFq/decodeFq -1 example/OsCas12f1-genome_R1.fq.gz -2 example/OsCas12f1-genome_R2.fq.gz -c OsCas12f1 -l example/library.json -o example/genome.tsv.gz
./editing_analysis/decodeFq/decodeFq -1 example/OsCas12f1-plasmid_R1.fq.gz -2 example/OsCas12f1-plasmid_R2.fq.gz -c OsCas12f1 -l example/library.json -o example/plasmid.tsv.gz

# calculate editing frequency (set the root path inside the script first)
Rscript editing_analysis/calculate_freq.R
```

---

## 2. eCROP

The eCROP folder contains the code for analyzing eCROP-seq data from the manuscript.

### 2.1 Installation

In addition to a Python environment (for sgRNA counting), R 4.5+ and the following R packages are required:

Package|Package|Package
---|---|---
corrgram|cowplot|data.table
dplyr|ggpubr|ggplot2
glue|glmnet|harmony
Matrix|reshape2|scales
Seurat|SeuratWrappers|stringr
tibble|tidyr|trqwe

The `count_sgrna_10x.py` script additionally requires the Python packages `pandas`, `pysam`, `click` and `rich`.

### 2.2 Usage Examples

Script|Description
---|---
count_sgrna_10x.py|Extract and quantify sgRNAs from 10x single-cell FASTQ/BAM files. Run `python count_sgrna_10x.py --help` for details.
cite_seq.Rmd|End-to-end eCROP-seq / CITE-seq analysis: Seurat object construction, QC, Harmony integration and clustering, cell-type annotation, per-cell sgRNA assignment, Fisher enrichment testing, STC guide correlation, sgRNA composition and UMAP density plots.
run_ridge.R|Ridge-regression model of gene expression on sgRNA assignments (standalone). Adjust the paths at the top of the script before running.

## License

Please refer to the repository for licensing information.

## Citation

If you use this code in your research, please cite the corresponding manuscript.

For questions or issues, please open an issue on GitHub.
