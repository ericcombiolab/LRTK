import sys

dbfile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

hash_dict = {}
cutoff = 0.001

with open(infile, 'r') as file:
    for line in file:
        line = line.strip()
        a = line.split()
        hash_dict[a[0]] = float(a[1])

with open(dbfile, 'r') as db_file, open(outfile, 'w') as out_file:
    for record in db_file.read().split('>'):
        if not record.strip():
            continue
        a = record.strip().split('\n')
        b = a[0].split('_')
        ref = '_'.join([b[0], b[1]])
        if ref in hash_dict:
            out_file.write('>' + '\n'.join(a) + '\n')

