#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="FASTQ file")
parser.add_argument("-o", "--output", help="FASTQ file with removed duplicates")
args = parser.parse_args()

def sort_fastq_by_sequence(input_fastq, output_fastq):
    # Read the input FASTQ file
    records = list(SeqIO.parse(input_fastq, "fastq"))
    
    # Sort records by sequence
    records.sort(key=lambda record: str(record.seq))
    
    # Write the sorted records to the output FASTQ file
    with open(output_fastq, "w") as output_handle:
        SeqIO.write(records, output_handle, "fastq")

sort_fastq_by_sequence(args.input, args.output)