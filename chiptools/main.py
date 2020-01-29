import sys, os
import numpy as np
import matplotlib.pyplot as plt
from .overlap import get_overlap, get_overlap_fraction
from .sizehist import get_hist, get_sizes

from .bedIO import get_chroms, print_chroms, read_bedfile, read_bedgraphs, read_peakfile
from .signalplot import signal_plot
from .regions import expand


def do_tss_plot():
    regions = read_bedfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs(sys.stdin)
    N = 1000
    signal = np.zeros(2*N)
    for chrom, bedgraph in bedgraphs:
        if chrom not in regions or chrom=="chrM":
            continue
        print("Reading", chrom)
        chrom_regions = expand(regions[chrom], N, N)
        signal += signal_plot(bedgraph, chrom_regions, 2*N)
    np.save(sys.argv[3], signal)
    if len(sys.argv)>4:
        plt.plot(signal)
        plt.savefig(sys.argv[4])

def do_averageplot():
    regions = read_peakfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs(sys.stdin)
    N=2000
    signal = np.zeros(2*N)
    for chrom, bedgraph in bedgraphs:
        print("Reading", chrom)
        if chrom not in regions or chrom=="chrM":
            continue
        chrom_regions = regions[chrom]
        signal += signal_plot(bedgraph, chrom_regions, 2*N, scale_to=True)
    np.save(sys.argv[3], signal)
    if len(sys.argv)>4:
        plt.plot(signal)
        plt.savefig(sys.argv[4])

def main():
    if sys.argv[1] == "overlap":
        data_a = get_chroms(open(sys.argv[2]))
        data_b = get_chroms(open(sys.argv[3]))
        overlap = {chrom: np.hstack((data_a[chrom], get_overlap(data_a[chrom], data_b[chrom])))
                   for chrom in data_a if chrom in data_b}
        print_chroms(overlap)
    elif sys.argv[1] == "overlap_fraction":
        data_a = get_chroms(open(sys.argv[2]))
        data_b = get_chroms(open(sys.argv[3]))
        overlap = {chrom: np.hstack((data_a[chrom], get_overlap_fraction(data_a[chrom], data_b[chrom])))
                   for chrom in data_a if chrom in data_b}
        print_chroms(overlap)

    elif sys.argv[1] == "tssplot":
        do_tss_plot()

    elif sys.argv[1] == "averageplot":
        do_averageplot()

    elif sys.argv[1] == "sizehist":
        bins = get_hist(get_sizes(open(sys.argv[2])),
                        int(sys.argv[3]), int(sys.argv[4]))
        plt.xlabel("size")
        plt.ylabel("count")
        plt.savefig(sys.argv[5])
