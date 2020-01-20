from collections import defaultdict
import numpy as np
dtype= [('start', np.unicode_, 16), ('grades', np.float64, (2,))]

def get_chroms(lines):
    chroms = defaultdict(list)
    for line in lines:
        chrom, start, end = line.split("\t", 3)[:3]
        chroms[chrom].append((int(start), int(end)))
    for chrom, coords in chroms.items():
        chroms[chrom] = np.array(coords, dtype="int")
    return {chrom: coords[coords[:, 0].argsort()] for chrom, coords in chroms.items()}

def print_chroms(chroms):
    for chrom, lines in chroms.items():
        for line in lines:
            print(chrom+"\t"+"\t".join(str(c) for c in line))
