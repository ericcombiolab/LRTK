import sys
import getopt

infile = ""
outfile = ""
genome = ""

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:g:")
except getopt.GetoptError:
    print("Usage: python script.py -i infile -o outfile -g genome")
    sys.exit(2)

for opt, arg in opts:
    if opt == "-i":
        infile = arg
    elif opt == "-o":
        outfile = arg
    elif opt == "-g":
        genome = arg

if not (infile and outfile and genome):
    print("Usage: python script.py -i infile -o outfile -g genome")
    sys.exit(2)

files = {}
with open(infile, "r") as file:
    for line in file:
        barcode = line.strip().split("/")[-2]
        files[barcode] = line.strip()

SNVs = {}
for barcode in sorted(files.keys()):
    infile = files[barcode]
    with open(infile, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith("#") and genome in line:
                columns = line.split("\t")
                key = columns[0] + "\t" + columns[1] + "\t" + columns[1]
                SNVs[key] = SNVs.get(key, 0) + 1

with open(outfile, "w") as file:
    for pos in sorted(SNVs.keys()):
        file.write(pos + "\n")
