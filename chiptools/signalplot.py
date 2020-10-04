import numpy as np

from .regions import Regions, expand

def get_tss_plot(regions, bedgraphs, N=1000, do_normalize=True):
    diffs = np.zeros(2*N)
    n_regions = 0
    coverage = 0
    for chrom, bedgraph in bedgraphs:
        if do_normalize:
            coverage += bedgraph.sum()
        if chrom not in regions:
            continue
        chrom_regions = expand(regions[chrom], N, N)
        signal_plot(bedgraph, chrom_regions, diffs)
        n_regions += chrom_regions.starts.size
    signal = np.cumsum(diffs)/n_regions
    if do_normalize:
        signal/=(coverage/1000000)
    return signal

def get_average_plot(bedgraphs, regions, N=2000, do_normalize=True):
    diffs = np.zeros(2*N)
    n_regions = 0
    coverage = 0
    for chrom, bedgraph in bedgraphs:
        if do_normalize:
            coverage += bedgraph.sum()
        if chrom not in regions or chrom=="chrM":
            continue
        chrom_regions = regions[chrom]
        sizes = chrom_regions.ends-chrom_regions.starts
        new_regions = Regions(chrom_regions.starts-sizes//2, chrom_regions.ends+sizes//2)
        signal_plot(bedgraph, new_regions, diffs, scale_to=True)
        n_regions += regions[chrom].starts.size
    signal = np.cumsum(diffs)/max(n_regions, 1)
    if do_normalize:
        signal /= (coverage/1000000)
    return signal

def signal_plot(bedgraph, regions, diffs, scale_to=False):
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    if scale_to:
        signals = signals.scale_x(diffs.size)
    signals.sum(axis=1).update_dense_diffs(diffs)

def signal_cumulative_hist(bedgraph, regions):
    H = np.zeros(1000)
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    for signal in signals:
        h = signal.hist()
        H[:h.size]+=signal.hist()
    return np.cumsum(H)
    
