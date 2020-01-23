from collections import defaultdict
from .bedgraph import BedGraph
from .regions import Regions

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

def read_bedfile(lines):
    chroms = defaultdict(list)
    for line in lines:
        chrom, start, end, _, _, direction = line.split("\t", 6)[:6]
        chroms[chrom].append((int(start), int(end), 1 if direction=="+" else -1))
    for chrom, coords in chroms.items():
        coords = np.array(coords, dtype="int")
        chroms[chrom] = coords[coords[:, 0].argsort()]
    return {chrom: Regions(coords[:, 0], coords[:, 1], coords[:, 2])
            for chrom, coords in chroms.items()}

def read_bedgraphs(lines):
    cur_chrom = None
    cur_idxs, cur_values = ([], [])
    for line in lines:
        chrom, start, _, value = line.split("\t", 4)[:4]
        if chrom != cur_chrom:
            if cur_chrom is not None:
                yield (cur_chrom, BedGraph(cur_idxs, cur_values))
            cur_chrom = chrom
            cur_idxs, cur_values = ([], [])
        cur_idxs.append(int(start))
        cur_values.append(float(value))
    if cur_idxs:
        yield (cur_chrom, BedGraph(cur_idxs, cur_values))
