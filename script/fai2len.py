import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as in_file, open(outfile, 'w') as out_file:
    for line in in_file:
        a = line.strip().split()
        out_file.write('\t'.join([a[0], a[1]]) + '\n')

