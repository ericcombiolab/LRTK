import sys




dbfile = "/datahome/comp/ericteam/cschaoyang/SOFTWARE/MIDAS2/DATABASE/uhgg_merge/UHGG.reference.fa"
infile = sys.argv[1]
outfile = sys.argv[2]
refanno = "/datahome/comp/ericteam/cschaoyang/SOFTWARE/MIDAS2/DATABASE/uhgg_merge/UHGG.reference.txt"
newanno = sys.argv[3]

hash_dict = {}
cutoff = 0.001

with open(infile, 'r') as file:
    for line in file:
        line = line.strip()
        a = line.split()
        if float(a[1]) > cutoff:
            hash_dict[a[0]] = float(a[1])

with open(refanno, 'r') as ref_file, open(newanno, 'w') as new_file:
    for line in ref_file:
        line = line.strip()
        a = line.split('\t')
        b = a[0].split('_')
        ref = '_'.join([b[0], b[1]])
        if ref in hash_dict:
            new_file.write(line + '\n')

with open(dbfile, 'r') as db_file, open(outfile, 'w') as out_file:
    for record in db_file.read().split('>'):
        if not record.strip():
            continue
        a = record.strip().split('\n')
        b = a[0].split('_')
        ref = '_'.join([b[0], b[1]])
        if ref in hash_dict:
            out_file.write('>' + '\n'.join(a) + '\n')

