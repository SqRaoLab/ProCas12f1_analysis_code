# splitFq

This tool splits mixed FASTQ files into separate files based on barcodes, 
assigning reads with the **upstream (forward) barcode** to **R1.fq.gz** and 
those with the **downstream (reverse) barcode** to **R2.fq.gz**.

## barcode

The barcode json, the key is output file prefix, inner key used to match input files, then list of both end barcodes

```json
{
  "OsCas12f1-plasmid":  {"0202-cas12f": [["ACGCTTAC", "GGGAAAGA"]]},
  "OsCas12f1-genome":   {"0202-cas12f": [["CTAATACT", "TGCCCGCC"]]},
  "enAsCas12f-plasmid": {"0821-Cas12f": [["ACGCTTAC", "GGGAAAGA"]]},
  "enAsCas12f-genome":  {"0821-Cas12f": [["TCCCAGAA", "AGAACGAT"]]}
}
```

### files

```json
{
    "0202-cas12f": "/rawdata/0202-cas12f/*.fastq.gz",
    "0821-Cas12f": "/rawdata/0821-Cas12f/*.fq.gz"
}
```

### Installation

```bash
# Go 1.24+ required
go build .
```

### Run

```bash
splitFq --help

Usage: splitFq [global options] 

Global options:
        -b, --barcode  Path to json file records barcode information
        -i, --input    Path to json file records fastq file paths
        -s, --scaffold Path to json file records scaffold, (for more checks, not necessary)
        -o, --output   Path to output directory (default: output.txt)
            --debug    Enable debug mode
        -v, --version  Print version
        -h, --help     Show this help
```

