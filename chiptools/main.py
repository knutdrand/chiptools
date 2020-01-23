import sys, os
from .overlap import get_overlap, get_overlap_fraction
from .bedIO import get_chroms, print_chroms, read_bedfile, read_bedgraphs
from .signalplot import signal_plot

import numpy as np
def signal_plot():
    regions = read_bedfile(sys.argv[2])
    bedgraphs = read_bedgraphs(sys.stdin)
    signal = np.zeros(1000)
    for chrom, bedgraph in bedgraphs:
        chrom_regions = regions[chrom]
        signal_plot(bedgraph, chrom_regions, signal)
    np.save(sys.argv[3], signal)

def main():
    if sys.argv[1] == "overlap":
        data_a = get_chroms(open(sys.argv[2]))
        data_b = get_chroms(open(sys.argv[3]))
        overlap = {chrom: np.hstack((data_a[chrom], get_overlap(data_a[chrom], data_b[chrom])))
                   for chrom in data_a if chrom in data_b}
        print_chroms(overlap)
    if sys.argv[1] == "overlap_fraction":
        data_a = get_chroms(open(sys.argv[2]))
        data_b = get_chroms(open(sys.argv[3]))
        overlap = {chrom: np.hstack((data_a[chrom], get_overlap_fraction(data_a[chrom], data_b[chrom])))
                   for chrom in data_a if chrom in data_b}
        print_chroms(overlap)

    if sys.argv[1] == "signal_plot":
        regions = read_bedfile(sys.argv[2])
        bedgraphs = read_bedgraphs(sys.stdin)
        signal = np.zeros(1000)
        for chrom, bedgraph in bedgraphs:
            chrom_regions = regions[chrom]
            signal
