#!/home/comp/cschaoyang/SOFTWARE/ANACONDA3/V2019/bin/python

import sys


if len(sys.argv) <= 2:
    print("usage: ", sys.argv[0], " <fasta file or index file> <region size> <outfile>")
    print("generates a list of freebayes/bamtools region specifiers on stdout")
    print("intended for use in creating cluster jobs")
    exit(1)

fasta_index_file = sys.argv[1]
if not fasta_index_file.endswith(".fai"):
    fasta_index_file = fasta_index_file + ".fai"

fasta_index_file = open(fasta_index_file)
region_size = int(sys.argv[2])
outfile = sys.argv[3]

OUT=open(outfile,"wt")
for line in fasta_index_file:

    fields = line.strip().split("\t")

    chrom_name = fields[0]
    chrom_length = int(fields[1])
    region_start = 0

    while region_start < chrom_length:
        start = region_start
        end = region_start + region_size
        if end > chrom_length:
            end = chrom_length
        OUT.write(chrom_name + ":" + str(region_start) + "-" + str(end) + "\n")
        region_start = end

OUT.close()

