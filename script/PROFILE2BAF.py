import sys
import argparse

parser = argparse.ArgumentParser(description="Calculate BAF from profile list")
parser.add_argument("-i", dest="infile", required=True, help="Input profile list file")
parser.add_argument("-o", dest="outfile", required=True, help="Output BAF file")
parser.add_argument("-g", dest="genome", required=True, help="Genome")
args = parser.parse_args()

infile = args.infile
outfile = args.outfile
genome = args.genome

hash = {}
with open(infile, "r") as list_file:
    for line in list_file:
        barcode = line.strip().split("/")[-1].replace(".profile.xls", "")
        hash[barcode] = line.strip()

snps = {}
MinorAlleleFreq = {}
MajorAllele = {}
MinorAllele = {}

for barcode in sorted(hash.keys()):
    infile = hash[barcode]

    with open(infile, "r") as in_file:
        for line in in_file:
            line = line.strip()
            if line.startswith("chrom") or not line.startswith(genome):
                continue

            a = line.split("\t")
            pos = a[0] + ":" + a[1]
            ref = a[2]
            major = a[3]
            minor = a[4]

            MajorAllele.setdefault(pos, {}).setdefault(major, 0)
            MajorAllele[pos][major] += 1

            MinorAllele.setdefault(pos, {}).setdefault(minor, 0)
            MinorAllele[pos][minor] += 1

            total = int(a[8]) + int(a[9]) + int(a[10]) + int(a[11])
            if total > 0:
                MinorAlleleFreq.setdefault(pos, {}).setdefault(barcode, {})
                MinorAlleleFreq[pos][barcode]["A"] = int(a[8]) / total
                MinorAlleleFreq[pos][barcode]["C"] = int(a[9]) / total
                MinorAlleleFreq[pos][barcode]["G"] = int(a[10]) / total
                MinorAlleleFreq[pos][barcode]["T"] = int(a[11]) / total

bases = ["A", "C", "T", "G"]
with open(outfile, "w") as out_file:
    out_file.write("Position\t")
    out_file.write("\t".join(sorted(hash.keys())))
    out_file.write("\n")

    for pos in sorted(MajorAllele.keys()):
        major_base = "A"
        major_count = 0

        for base in bases:
            if base in MajorAllele[pos]:
                count = MajorAllele[pos][base]
                if count > major_count:
                    major_count = count
                    major_base = base

        minor_base = "A"
        minor_count = 0

        for base in bases:
            if base in MinorAllele[pos]:
                count = MinorAllele[pos][base]
                if count > minor_count:
                    minor_count = count
                    minor_base = base

        BAFvec = "NA"
        BAFcheck = 0

        for barcode in sorted(hash.keys()):
            if minor_base in MinorAlleleFreq[pos].get(barcode, {}):
                freq = MinorAlleleFreq[pos][barcode][minor_base]
            else:
                freq = 0

            if BAFvec == "NA":
                BAFvec = str(freq)
            else:
                BAFvec += "\t" + str(freq)

            BAFcheck += freq

        if minor_base != major_base and BAFcheck > 0:
            out_file.write(pos + "\t" + BAFvec + "\n")
