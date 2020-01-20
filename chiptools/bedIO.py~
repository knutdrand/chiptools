from collections import defaultdict
dtype= [('start', np.unicode_, 16), ('grades', np.float64, (2,))]

def get_chroms(lines,  buffer_size=10000):
    chroms = defaultdict(list)
    for line in lines:
        chrom, start, end, _ = line.split("\t", 3)
        chroms[chrom].append((int(start), int(end), attrs))
    for chrom, coords in chroms.items():
        chroms[chrom] = np.array(coords, dtype="int")
    return {chrom: coords[coords[:, 0].argsort()] for chrom, coords in chroms.items()}

def print_chroms(chroms):
    for chrom, lines in chroms.items():
        for line in lines:
            print(chrom+"\t"+"\t".join(line))
