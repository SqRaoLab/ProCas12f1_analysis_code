# decodeFq

Used to demultiplex each read pair into distinct components according to the designed library structure.

The input FASTQ files must be processed with splitFq to ensure that reads are correctly assigned to their respective files.


## Installation

```bash
# Go 1.24+ required
go build .
```

## Run

The running parameters listed:

```bash
decodeFq --help

Usage: decodeFq [global options] 

Global options:
        -1, --r1              Path to R1 FASTQ file (gzip)
        -2, --r2              Path to R2 FASTQ (gzip)
        -c, --cas             The cas protein name
        -l, --library         Path to csv (or gzipped csv) file of pool library
        -o, --output          The output path (default: output.txt)
        -p, --process         The number of goroutine to use (default: number of CPUs) (default: 10)
            --target-length   The length of gDNA target (default: 20)
            --primer-length   The length of Reverse Primer (default: 20)
            --pam-length      The length of PAM (default: 4)
            --spacer-length   The length of spacer (default: 20)
            --spacer-distance The distance between spacer and before target (default: 6)
            --before-length   The length of before target (default: 6)
            --behind-length   The length of behind target (default: 6)
            --cmd             The command line tools for read special compressed input file
        -U, --umi-anchor      The designed umi sequence
        -u, --umi             The length of umi (default: 28)
            --debug           Enable debug mode
        -v, --version         Print version
        -h, --help            Show this help
```
