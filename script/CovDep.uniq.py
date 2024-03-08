import sys

faifile = sys.argv[1]
depfile = sys.argv[2]
outfile = sys.argv[3]

hash = {}
with open(faifile, 'r') as file:
    for line in file:
        a = line.strip().split('\t')
        b = a[0].split('_')
        genome = '_'.join(b[:2])

        if genome not in hash:
            hash[genome] = int(a[1])
        else:
            hash[genome] += int(a[1])

coverage = {}
base = {}
with open(depfile, 'r') as dep_file, open(outfile, 'w') as out_file:
    for line in dep_file:
        a = line.strip().split('\t')
        b = a[0].split('_')
        genome = '_'.join(b[:2])

        if genome not in coverage:
            coverage[genome] = 1
            base[genome] = int(a[2])
        else:
            coverage[genome] += 1
            base[genome] += int(a[2])

    for genome in sorted(coverage.keys()):
        ratio = coverage[genome] / hash[genome]
        m_dep = base[genome] / coverage[genome]
        if ratio > 0.4 and coverage[genome] > 500000:
            out_file.write(f"{genome}\t{coverage[genome]}\t{hash[genome]}\t{ratio}\t{m_dep}\n")

