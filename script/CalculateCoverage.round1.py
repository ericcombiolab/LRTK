import sys

hash_dict = {}
coverage_dict = {}
read_len = 150

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('#'):
            continue
        
        cols = line.split('\t')
        ref = cols[0].split('_')[0] + '_' + cols[0].split('_')[1]
        length = int(cols[2])
        avg_depth = float(cols[6]) * read_len / length
        
        if ref in hash_dict:
            hash_dict[ref].append(avg_depth)
        else:
            hash_dict[ref] = [avg_depth]

for ele in sorted(hash_dict):
    array = hash_dict[ele]
    array.sort()
    cutoff_min = int(0.2 * (len(array) - 1))
    cutoff_max = int(0.8 * (len(array) - 1))
    depth_total = sum(array[cutoff_min:cutoff_max+1])
    count = cutoff_max - cutoff_min + 1
    depth_mean = depth_total / count
    coverage_dict[ele] = depth_mean

total_cov = 0
for ele in sorted(coverage_dict):
    if coverage_dict[ele] == 0 or coverage_dict[ele] < 10:
        continue
    total_cov += coverage_dict[ele]

with open(outfile, 'w') as file:
    for ele in sorted(coverage_dict):
        if coverage_dict[ele] == 0 or coverage_dict[ele] < 10:
            continue
        relative_abundance = coverage_dict[ele] / total_cov
        if relative_abundance < 0.0001:
            continue
        file.write(f"{ele}\t{relative_abundance}\n")

