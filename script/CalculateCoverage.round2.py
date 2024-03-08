import sys

infile = sys.argv[1]
outfile = sys.argv[2]

hash = {}
coverage = {}
read_length = 100

with open(infile, 'r') as in_file:
    for line in in_file:
        line = line.strip()
        if line.startswith('#'):
            continue
        a = line.split('\t')
        b = a[0].split('_')
        ref = b[0] + '_' + b[1]
        length = int(a[3])
        avg_depth = float(a[6]) * read_length / length

        if ref in hash:
            hash[ref] += ',' + str(avg_depth)
        else:
            hash[ref] = str(avg_depth)

total_cov = 0
for ele in sorted(hash.keys()):
    array = list(map(float, hash[ele].split(',')))
    array_sort = sorted(array)
    cutoff_min = int(0.2 * len(array))
    cutoff_max = int(0.8 * len(array))
    depth_total, count = 0, 0
    for i in range(cutoff_min, cutoff_max + 1):
        depth_total += array_sort[i]
        count += 1
    depth_mean = depth_total / count
    coverage[ele] = depth_mean
    total_cov += depth_mean

with open(outfile, 'w') as out_file:
    out_file.write("#Genome\tDepth\tRelativeAbundance\n")
    for ele in sorted(coverage.keys()):
        if coverage[ele] == 0:
            continue
        ratio = coverage[ele] / total_cov
        out_file.write(f"{ele}\t{coverage[ele]}\t{ratio}\n")

