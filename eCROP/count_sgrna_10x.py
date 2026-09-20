#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract sgRNAs from 10x data and assign them to cells and reads.

Requires the transcriptome R2 reads and the aligned BAM file in order to
recover the correct cell barcodes.
"""
import os
import csv
import re
import gzip
import pandas as pd

import click

import pysam

from rich import print
from rich.progress import track


# The first part is a custom scaffold used to locate the sgRNA.
# The captured sgRNA is 20 bp long.
PATTERN = re.compile(r"GTGCATGAGCCGCGAAAGCGGCTTGAAGG(?P<sgrna>[ATCG]{20})")
BARCODE_LEN = 16

def load_barcodes(path: str):
    res = set()

    with gzip.open(path, "rt") as r:
        for line in r:
            res.add(line.strip())
    return res


def complement(seq: str):
    codon = {
        "A": "T", "C": "G",
        "T": "A", "G": "C",
        "N": "N"
    }
    res = [codon[x] for x in seq]
    return "".join(res[::-1])


def extract_barcode_and_umi_from_bam(path: str,  barcodes: set, sgrnas: dict, output: str):
    """
    Read the sequence name together with the cell barcode and UMI barcode from a BAM file.
    CB:Z:CACCATTTCTATGGAT-1  corrected cell barcode, -1 is the GEM well index
    UB:Z:CGCACCGCAAGA         corrected UMI sequence, used for molecule counting
    """
    res = {}
    with pysam.AlignmentFile(path, "rc") as bam:
        for rec in track(bam, description="Reading bam..."):
            if rec.has_tag("CB") and rec.has_tag("UB"):

                if rec.get_tag("CB") in barcodes and rec.query_name in sgrnas.keys():
                    temp = sgrnas[rec.query_name]
                    temp.update({"cb": rec.get_tag("CB")})
                    key = temp["cb"] + "|" + temp["sgrna"] + "|" + temp["gene"]
                    res[key] = res.get(key, 0) + 1


    with gzip.open(output, "wt+", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["cb", "gene", "sgrna", "value"])
        writer.writeheader()

        for key, val in res.items():
            cb, gene, sgrna = key.split("|")
            temp = {
                "cb": cb, "gene": gene, "sgrna": sgrna, "value": val
            }
            writer.writerow(temp)

def load_sgrna_library(path: str):
    """
    Read the TF sgRNA library.
    """
    sgrna_libraries = {}
    ref = pd.read_excel(path, index_col=0)
    for _, row in ref.iterrows():
        sgrna_libraries[row["sgRNA sequence"].strip()] = {"sgrna": row["sgRNA"], "gene": row["gene"]}

    return sgrna_libraries


def extract_sgrna_from_fastq(path: str, sgrna_libraries: dict, tempfile: str):
    u"""
    :params path: path to transcriptome fastq file
    """

    with gzip.open(tempfile, 'wt+', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["name", "seq", "gene", "sgrna"])
        writer.writeheader()  # write header

        with pysam.FastxFile(path) as fin:
            for rec in track(fin, description="Reading fastq..."):
                m = PATTERN.search(rec.sequence)
                if m:
                    m = m.groupdict()

                    temp = {"name": rec.name, "seq": m["sgrna"]}

                    sgrna = sgrna_libraries.get(temp["seq"])
                    if sgrna:
                        temp["gene"] = sgrna["gene"]
                        temp["sgrna"] = sgrna["sgrna"]
                        writer.writerow(temp)


def decode_sequence(seq1: str, seq2: str):
    """ Extract the barcode and sgRNA from a read pair. """
    barcode = seq1[:BARCODE_LEN]
    m = PATTERN.search(seq2)
    if m:
        m = m.groupdict()
        m = m["sgrna"]

    return barcode, m


def read_pair(r1: str, r2: str):
    """ Read paired-end reads. """
    r1s, r2s = {}, {}
    with pysam.FastxFile(r1) as f1, pysam.FastxFile(r2) as f2:
        for read1, read2 in zip(f1, f2):
            if read1.name == read2.name:
                yield read1.sequence, read2.sequence
            elif read2.name in r1s:
                r1s[read1.name] = read1
                read1 = r1s.pop(read2.name)
                yield read1.sequence, read2.sequence
            elif read1.name in r2s:
                r2s[read2.name] = read2
                read2 = r2s.pop(read1.name)
                yield read1.sequence, read2.sequence


def extract_barcode_and_sgrna_from_fastq(r1: str, r2: str, sgrna_libraries: dict, barcodes: set, output: str):
    """
    Quantify sgRNAs from the 10x sgRNA fastq files following the library layout:
    - the first 16 bp of R1 is the cell barcode
    - the 20 bp of R2 after the scaffold is the sgRNA
    """

    res = {}
    for seq1, seq2 in track(read_pair(r1, r2)):
        barcode, sgrna = decode_sequence(seq1, seq2)
        if sgrna and sgrna in sgrna_libraries and barcode + "-1" in barcodes:
            key = f"{barcode}|{sgrna}"
            res[key] = res.get(key, 0) + 1

    with open(output, "w+") as w:
        writer = csv.DictWriter(w, fieldnames=["barcode", "sgrna", "gene", "value"])
        writer.writeheader()

        for key, val in res.items():
            barcode, seq = key.split("|")

            temp = sgrna_libraries[seq]
            temp["barcode"] = barcode
            temp["value"] = val
            writer.writerow(temp)


@click.command()
@click.option("-b", "--bam", type=click.Path(exists=True), help="bam file")
@click.option("-1", "--r1", type=click.Path(exists=True), help="r2 file")
@click.option("-2", "--r2", type=click.Path(exists=True), help="r2 file")
@click.option("-s", "--sgrna", type=click.Path(exists=True), help="sgrna library")
@click.option("-o", "--output", type=click.Path(), help="output path")
@click.option("-c", "--barcode", type=click.Path(exists=True), help="barcode whitelist file")
def main(bam: str, r1: str, r2: str, sgrna: str, output: str, barcode: str):
    """
    Extract the sgRNA (and related information) for each read from the Cell Ranger
    BAM file and the R2 reads of the raw fastq files.
    """

    print("Loading sgRNA library")
    libraries = load_sgrna_library(sgrna)
    barcodes = load_barcodes(barcode)

    if bam:
        print("Loading read names and their sgRNAs")
        tempfile = bam.replace(".cram", "_sgrna.tsv.gz")
        if not os.path.exists(tempfile) or os.path.getsize(barcode) < 100:
            extract_sgrna_from_fastq(r2, libraries, tempfile)

        sgrnas = {}
        with gzip.open(tempfile, "rt", encoding="utf-8") as r:
            reader = csv.DictReader(r, delimiter=",")
            for rec in reader:
                sgrnas[rec["name"]] = rec

        print(f"Load UB and CB from bam")
        extract_barcode_and_umi_from_bam(bam, barcodes, sgrnas, output)

    elif r1:
        print(f"Load UB and CB from r1 and r2")
        extract_barcode_and_sgrna_from_fastq(r1, r2, libraries, barcodes, output)
    pass



if __name__ == "__main__":
    main()
