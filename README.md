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

# calculate editing frequency
./editing_analysis/calFreq/calFreq -i example/genome.tsv.gz -b example/plasmid.tsv.gz -o example/indel_freq.tsv.gz
```

---

## 2. eCROP

The eCROP folder contains the code for analyzing eCROP-seq data from the manuscript.

### 2.1 Installation

In addition to the Python environment, R 4.5+ and the following R libraries are required:

Package|Package|Package
---|---|---
corrplot|data.table|dplyr
forcasts|ggplot2|ggsci
GGally|gridExtra|Matrix
pheatmap|purrr|qs
readr|scales|Seurat
SingleR|tibble|tidyr
tidyverse|harmony| 

### 2.2 Usage Examples

Script|Description
---|---
crop_count.py|Count sgRNAs from eCROP-seq experiments. Use python crop_count.py -h for help.
sgrna_assign.R|Integrate gRNA assignments from read and transcriptome data, filter low-quality cells based on QC metrics, and append final assignments to the gene-count matrix.
eCROP_UMAPs.R|Load aggregated scRNA-seq data into R using Seurat, integrate sgRNA assignments into the object metadata, and perform data processing to generate Figure 6 of the manuscript.
eCROP_Seurat.R|Additional Seurat-based analysis for eCROP data.
eCROP_ClusterEnrichments.R|Cluster enrichment analysis for eCROP data.
cite_seq.Rmd|Rmarkdown file for CITE-seq data analysis.

## License

Please refer to the repository for licensing information.

## Citation

If you use this code in your research, please cite the corresponding manuscript.

For questions or issues, please open an issue on GitHub.
