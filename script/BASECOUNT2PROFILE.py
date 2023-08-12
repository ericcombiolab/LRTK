import sys
import argparse

parser = argparse.ArgumentParser(description='Process base counts and generate a profile.')
parser.add_argument('-i', '--infile', help='Input file (base count)', required=True)
parser.add_argument('-o', '--outfile', help='Output file (profile)', required=True)
args = parser.parse_args()

infile = args.infile
outfile = args.outfile


with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
    f_out.write("\t".join(["chrom", "pos", "ref_allele", "major_allele", "minor_allele", "ref_freq", "major_freq", "minor_freq", "A", "C", "G", "T"]) + "\n")

    for line in f_in:
        line = line.strip()
        a = line.split("\t")

        chrom = a[0]
        pos = a[1]
        ref_allele = a[2]
        major_allele = ref_allele
        minor_allele = ref_allele
        major_count = 0
        minor_count = 0

        basecount = {}
        for i in range(5, 9):
            b = a[i].split(":")
            basecount[b[0]] = b[1]

        for base in sorted(basecount, key=lambda x: basecount[x], reverse=True):
            if int(basecount[base]) > int(major_count):
                major_allele = base
                major_count = basecount[base]
            elif int(basecount[base]) > int(minor_count):
                minor_allele = base
                minor_count = basecount[base]

        total_count = int(major_count) +  int(minor_count)
        ref_count = int(basecount[ref_allele])

        ref_freq = 0
        major_freq = 0
        minor_freq = 0
        if total_count > 0:
            ref_freq = round(int(ref_count) / total_count, 2)
            major_freq = round(int(major_count) / total_count, 2)
            minor_freq = round(int(minor_count) / total_count, 2)

        f_out.write("\t".join([chrom, pos, ref_allele, major_allele, minor_allele, str(ref_freq), str(major_freq), str(minor_freq), str(basecount.get("A", 0)), str(basecount.get("C", 0)), str(basecount.get("G", 0)), str(basecount.get("T", 0))]) + "\n")
