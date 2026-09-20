package main

import "github.com/voxelbrain/goptions"

type params struct {
	R1             string        `goptions:"-1, --r1, description='Path to R1 FASTQ file (gzip)'"`
	R2             string        `goptions:"-2, --r2, description='Path to R2 FASTQ (gzip)'"`
	Cas            string        `goptions:"-c, --cas, description='The cas protein name'"`
	Pool           string        `goptions:"-l, --library, description='Path to csv (or gzipped csv) file of pool library'"`
	Output         string        `goptions:"-o, --output, description='The output path'"`
	Procs          int           `goptions:"-p, --process, description='The number of goroutine to use (default: number of CPUs)'"`
	TargetLen      int           `goptions:"--target-length, description='The length of gDNA target'"`
	PrimerLen      int           `goptions:"--primer-length, description='The length of Reverse Primer'"`
	PAMLen         int           `goptions:"--pam-length, description='The length of PAM'"`
	SpacerLen      int           `goptions:"--spacer-length, description='The length of spacer'"`
	SpacerDistance int           `goptions:"--spacer-distance, description='The distance between spacer and before target'"`
	BeforeLen      int           `goptions:"--before-length, description='The length of before target'"`
	BehindLen      int           `goptions:"--behind-length, description='The length of behind target'"`
	ReadCmd        string        `goptions:"--cmd, description='The command line tools for read special compressed input file'"`
	UMIAnchor      string        `goptions:"-U, --umi-anchor, description='The designed umi sequence'"`
	UMILen         int           `goptions:"-u, --umi, description='The length of umi'"`
	Debug          bool          `goptions:"--debug, description='Enable debug mode'"`
	Version        bool          `goptions:"-v, --version, description='Print version'"`
	Help           goptions.Help `goptions:"-h, --help, description='Show this help'"`
}

func defaultParams() *params {
	return &params{
		Output:         "output.txt",
		Procs:          10,
		TargetLen:      20,
		PrimerLen:      20,
		PAMLen:         4,
		SpacerLen:      20,
		SpacerDistance: 6,
		BeforeLen:      6,
		BehindLen:      6,
		UMILen:         28,
	}
}
