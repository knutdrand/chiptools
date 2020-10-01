import numpy as np
import logging
from collections import namedtuple

log = logging

Size = namedtuple("Size", ["x", "y"])

def vplot(bedgraph, regions, size_x, size_y, max_size):
    sizes = regions.ends-regions.starts
    diffs = np.zeros((size_y, size_x))
    mids = (regions.ends+regions.starts)//2
    signals = bedgraph.get_slices(mids-max_size//2, mids+max_size//2, regions.directions)
    Ns = np.zeros(size_y)
    for size, signal in zip(sizes, signals):
        if size>max_size:
            continue
        assert np.all(signal._values >= 0)
        graph_diffs = signal.scale_x(size_x).to_graph_diffs()
        graph_diffs.assert_positive()
        row = min(int(size/max_size*size_y), diffs.shape[0]-1)
        assert np.all(np.cumsum(diffs[row])>=0), np.where(np.cumsum(diffs[row])<0)
        graph_diffs.update_dense_array(diffs[row])
        assert np.all(np.cumsum(diffs[row])>=0), (graph_diffs, np.cumsum(diffs[row]))
        Ns[row]+=1
    tot = np.cumsum(diffs, axis=1)
    assert np.all(tot>=0)
    return tot, Ns

def get_ranks(array):
    args_tmp = np.argsort(array)
    args = np.empty_like(args_tmp)
    args[args_tmp] = np.arange(len(args))
    return args

def get_y_coords(regions, size_y):
    sizes = [r.ends-r.starts for _, r in regions.items()]
    offsets = np.cumsum([0]+[len(s) for s in sizes])
    ranks = get_ranks(np.concatenate(sizes))
    y_coords = (ranks/ranks.size*size_y).astype("int")
    return {chrom: y_coords[offsets[i]:offsets[i+1]] for i, chrom in enumerate(regions)}

def get_heatplot(regions, bedgraphs, max_size=50000, fig_size=Size(2000, 4000)):
    y_coords = get_y_coords(regions, fig_size.y)
    diffs = np.zeros((fig_size.y, fig_size.x))
    Ns = np.zeros(fig_size.y)
    for chrom, bedgraph in bedgraphs:
        if chrom not in regions or chrom=="chrM":
            continue
        assert "alt" not in chrom, chrom
        log.info("Reading %s", chrom)
        heatplot_per_chrom(bedgraph, regions[chrom], y_coords[chrom], diffs, max_size)
        for y in y_coords[chrom]:
            Ns[y] += 1

    signal = np.cumsum(diffs, axis=1)/(np.maximum(Ns, 1)[:, None])
    return signal


def heatplot_per_chrom(bedgraph, regions, y_coords, diffs, max_size):
    mids = (regions.ends+regions.starts)//2
    signals = bedgraph.get_slices(mids-max_size//2, mids+max_size//2, regions.directions)
    graph_diffs = (signal.scale_x(diffs.shape[1]).to_graph_diffs() for signal in signals)
    for y, signal in zip(y_coords, graph_diffs):
        signal.update_dense_array(diffs[y])
