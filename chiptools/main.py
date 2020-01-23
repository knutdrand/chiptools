import sys, os
import numpy as np
from .overlap import get_overlap, get_overlap_fraction
from .bedIO import get_chroms, print_chroms, read_bedfile, read_bedgraphs
from .signalplot import signal_plot
from .regions import expand


def do_signalplot():
    regions = read_bedfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs(sys.stdin)
    N = 1000
    signal = np.zeros(2*N)
    for chrom, bedgraph in bedgraphs:
        print("Reading", chrom)
        chrom_regions = expand(regions[chrom], N, N)
        signal += signal_plot(bedgraph, chrom_regions, 2*N)
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

    if sys.argv[1] == "signalplot":
        do_signalplot()
