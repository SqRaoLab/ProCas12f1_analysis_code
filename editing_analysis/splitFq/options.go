package main

import "github.com/voxelbrain/goptions"

type params struct {
	Barcode       string        `goptions:"-b, --barcode, description='Path to json file records barcode information'"`
	Pool          string        `goptions:"-i, --input, description='Path to json file records fastq file paths'"`
	Output        string        `goptions:"-o, --output, description='Path to output directory'"`
	Merge         bool          `goptions:"-m, --merge, description='Merge records into single fq file'"`
	Trim          bool          `goptions:"-t, --trim, description='Trim barcode'"`
	MinMeanQual   int           `goptions:"--min-mean-quality, description='min mean quality of reads'"`
	MinQ30Percent float32       `goptions:"--min-q20-perc, description='min percentage of valid Q30 score in single reads'"`
	Debug         bool          `goptions:"--debug, description='Enable debug mode'"`
	Version       bool          `goptions:"-v, --version, description='Print version'"`
	Help          goptions.Help `goptions:"-h, --help, description='Show this help'"`
}

func defaultParams() *params {
	return &params{
		Output:        "output",
		MinMeanQual:   30,
		MinQ30Percent: 90,
	}
}
