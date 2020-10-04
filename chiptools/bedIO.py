from collections import defaultdict, namedtuple
from itertools import repeat, takewhile, chain, groupby
import logging
from operator import itemgetter
from .bedgraph import BedGraph
from .regions import Regions

import numpy as np
import pandas as pd
dtype= [('start', np.unicode_, 16), ('grades', np.float64, (2,))]

BedEntry=namedtuple("BedEntry", ["chrom", "start", "end", "strand"])
log = logging
def vanillabed(lines):
    parted = (line.split() for line in lines)
    return (BedEntry(chrom, int(start), int(end), strand.strip())
            for chrom, start, end, _, _, strand, *_ in parted)

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
        chrom, start, end, _, _, direction = line.strip().split("\t", 6)[:6]
        assert direction in ("+", "-"), (line, direction)
        chroms[chrom].append((int(start), int(end), 1 if direction=="+" else -1))
    for chrom, coords in chroms.items():
        coords = np.array(coords, dtype="int")
        chroms[chrom] = coords[coords[:, 0].argsort()]
    return {chrom: Regions(coords[:, 0], coords[:, 1], coords[:, 2])
            for chrom, coords in chroms.items()}

def read_peakfile(lines):
    chroms = defaultdict(list)
    for line in lines:
        chrom, start, end  = line.strip().split("\t", 3)[:3]
        chroms[chrom].append((int(start), int(end)))
    for chrom, coords in chroms.items():
        coords = np.array(coords, dtype="int")
        chroms[chrom] = coords[coords[:, 0].argsort()]
    return {chrom: Regions(coords[:, 0], coords[:, 1], 1)
            for chrom, coords in chroms.items()}


def read_fragments(lines):
    cur_chrom = None
    cur_starts, cur_ends = ([], [])
    for line in lines:
        chrom, start, end = line.split("\t", 3)[:3]
        if chrom != cur_chrom:
            if cur_chrom is not None:
                yield (cur_chrom, Regions(cur_starts, cur_ends))
            cur_chrom = chrom
            cur_starts, cur_ends = ([], [])
        cur_starts.append(int(start))
        cur_ends.append(int(end))
    if cur_starts:
        yield (cur_chrom, Regions(cur_starts, cur_ends))

def read_bedgraphs_pd(file_obj, size_hint=1000000):
    cur_chrom=None
    reader = pd.read_table(file_obj, names=["chrom", "start", "end", "value"], usecols=[0, 1, 2, 3], chunksize=size_hint)
    chunks = []
    for chunk in reader:
        while chunk["chrom"].iloc[-1] != cur_chrom:
            idx = np.argmax(chunk["chrom"].values!=cur_chrom)
            chunks.append(chunk.iloc[:idx])
            if cur_chrom is not None:
                log.info("Yielding chrom %s", cur_chrom)
                yield cur_chrom, BedGraph(np.concatenate([c["start"].values for c in chunks]),
                                          np.concatenate([c["value"].values for c in chunks]),
                                          chunks[-1]["end"].values[-1])
                chunks = []
            chunk = chunk.iloc[idx:]
            cur_chrom = chunk["chrom"].iloc[0]
        chunks.append(chunk)
    yield cur_chrom, BedGraph(np.concatenate([c["start"].values for c in chunks]),
                              np.concatenate([c["value"].values for c in chunks]),
                              chunks[-1]["end"].values[-1])


def read_bedgraphs_fast(file_obj, size_hint=10000000):
    chunks = (file_obj.readlines(size_hint) for _ in repeat(None))
    chunks = takewhile(lambda x: x, chunks)
    lines =  chain.from_iterable(chunks)
    all_parted = (line.split("\t", 4)[:4] for line in lines)
    for chrom, parted in groupby(all_parted, itemgetter(0)):
        if chrom == "chrM":
            continue
        assert "alt" not in chrom, chrom
        idxs, ends, values = zip(*((int(p[1]), int(p[2]), float(p[3])) for p in parted))
        yield chrom, BedGraph(np.array(idxs, dtype="int"), np.array(values, dtype="float"), size=ends[-1])

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

def print_regions(chrom, regions):
    for start, end, _ in regions:
        print("%s\t%s\t%s" % (chrom, start, end))
