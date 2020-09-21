import matplotlib.pyplot as plt
import numpy as np
bin_size = 0.05
n_bins = 250

def get_log_hist(sizes):
    bins = np.zeros(n_bins, dtype="int")
    log_sizes = (np.log(size) for size in sizes)
    size_bins = (min(s//bin_size, n_bins-1) for s in log_sizes)
    for sb in size_bins:
        bins[int(sb)] += 1
    N = 110
    ticks = np.arange(0, bins.size-N, (bins.size-N)//7)
    tick_values = np.array(np.exp((ticks+N)*bin_size), dtype="int")
    plt.plot(bins[N:])
    plt.xticks(ticks, tick_values)
    return bins

def get_size(p):
    assert len(p) >= 3, p
    return int(p[2])-int(p[1])

def get_sizes(lines):
    parts = (l.split() for l in lines)
    return (get_size(p) for p in parts)
    
def get_hist(sizes, n_bins=50, bin_size=1000):
    bins = np.zeros(n_bins, dtype="int")
    size_bins = (min(s//bin_size, n_bins-1) for s in sizes) 
    for sb in size_bins:
        bins[int(sb)] += 1
    plt.plot(np.arange(bins.size)*bin_size, bins)
    return bins


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        input_lines = open(sys.argv[2])
    else:
        input_lines = sys.stdin
    out_name = sys.argv[-1]
    sizes = get_sizes(input_lines)
    funcs[sys.argv[1]](sizes)
    plt.xlabel("size")
    plt.ylabel("count")
    plt.savefig(out_name)
