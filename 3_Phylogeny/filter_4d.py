#!/usr/bin/env python3
# filter_vcf_by_4d.py
from sys import argv
from pysam import VariantFile

file_in = argv[1]
file_out = argv[2]
sites = argv[3]

sets = set([line.split('\t')[0] + line.split('\t')[1] for line in open(sites, "r")])

bcf_in = VariantFile(file_in)
bcf_out = VariantFile(file_out, "w", header = bcf_in.header)
for rec in bcf_in.fetch():
    cur = rec.chrom + str(rec.pos) + '\n'
    if cur not in sets:
        continue
    bcf_out.write(rec)
