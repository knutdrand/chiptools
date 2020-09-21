import pandas as pd
import numpy as np
from .bedgraph import BedGraph
read_bed = lambda filename: pd.read_csv(filename, sep="\t", usecols=(0,1,2,5), engine="c", names=("chrom", "start", "end", "direction"))

def pileup(starts, ends):
    indices = np.concatenate((starts, ends))
    args = np.argsort(indices, kind="mergesort")
    diffs = np.where(args>=starts.size, -1, 1)
    indices = indices[args]
    duplicates = np.flatnonzero(indices[1:]==indices[:-1])
    print(np.cumsum(diffs))
    return BedGraph(np.delete(indices, duplicates), np.delete(np.cumsum(diffs),duplicates), strict=False)

def find_enriched(starts, ends,  n_reads=2, ):
    starts = starts.to_numpy()
    print(pileup(starts, starts+window).threshold(n_reads))
    is_enriched = starts[n_reads:]-starts[:-n_reads] < window
    changes = np.flatnonzero(is_enriched[1:] != is_enriched[:-1])
    starts = starts[n_reads:][changes[::2]+1]
    ends = starts[n_reads:][changes[1::]-n_reads]+window
    return starts, ends
    
def fraglen(filename, window=300, n_reads=3):
    intervals = read_bed(filename)
    is_forward = intervals["direction"] == "+"
    directed_intervals = (intervals[is_forward]["start"], intervals[~is_forward]["end"])
    pos_starts = intervals[is_forward]["start"].to_numpy()
    neg_ends = intervals[~is_forward]["end"].to_numpy()
    forward_peaks = pileup(pos_starts, pos_starts+window).threshold(n_reads)
    reverse_peaks = pileup(neg_ends-window, neg_ends).threshold(n_reads)
    print(forward_peaks)
    print(reverse_peaks)
    

