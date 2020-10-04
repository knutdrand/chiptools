import numpy as np
import logging
from collections import namedtuple

log = logging

Size = namedtuple("Size", ["x", "y"])

def get_vplot(regions, bedgraphs, max_size=50000, fig_size=Size(2000, 2000)):
    diffs = np.zeros((fig_size.y, fig_size.x))
    Ns = np.zeros(fig_size.y)
    total = 0
    for chrom, bedgraph in bedgraphs:
        total += bedgraph.sum()
        if chrom not in regions:
            continue
        log.info("Processing %s", chrom)
        vplot(bedgraph, regions[chrom], diffs, max_size, Ns)

    return np.cumsum(diffs, axis=1)/(np.maximum(Ns, 1)[:, None])/total*1000000

def vplot_(bedgraph, regions, diffs, max_size, Ns):
    rows = ((regions.ends-regions.starts)/max_size*diffs.shape[0]).astype("int")
    mask = rows<diffs.shape[0]
    mids = (regions.ends[mask]+regions.starts[mask])//2
    bedgraph.get_slices(mids-max_size//2, mids+max_size//2, regions.directions[mask]).scale_x(diffs.shape[1]).update_dense_diffs(diffs, rows[mask])
    rows, counts = np.unique(rows[mask], return_counts=True)
    Ns[rows]+=counts
    # for row in rows[mask]:
    # Ns[row] += 1
    #for row, signal in zip(rows[mask], signals):
    # assert np.all(signal._values >= 0)
    # signal.update_dense_diffs(diffs[row])

def vplot(bedgraph, regions, diffs, max_size, Ns):
    rows = ((regions.ends-regions.starts)/max_size*diffs.shape[0]).astype("int")
    mask = rows<diffs.shape[0]
    print(mask.dtype, mask)
    mids = (regions.ends[mask]+regions.starts[mask])//2
    print(regions.directions)
    dirs = regions.directions[mask]
    signals = bedgraph.get_slices_normal(mids-max_size//2, mids+max_size//2, dirs)#regions.directions[mask])# .scale_x(diffs.shape[1])
    for row, signal in zip(rows[mask], signals):
        assert np.all(signal._values >= 0)
        signal.scale_x(diffs.shape[1]).update_dense_diffs(diffs[row])
        Ns[row] += 1

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
    do_normalize=True
    y_coords = get_y_coords(regions, fig_size.y)
    diffs = np.zeros((fig_size.y, fig_size.x))
    Ns = np.zeros(fig_size.y)
    total = 0
    for chrom, bedgraph in bedgraphs:
        total+=bedgraph.sum()
        if chrom not in regions:
            continue
        log.info("Reading %s", chrom)
        heatplot_per_chrom(bedgraph, regions[chrom], y_coords[chrom], diffs, max_size)
        for y in y_coords[chrom]:
            Ns[y] += 1

    signal = np.cumsum(diffs, axis=1)/(np.maximum(Ns, 1)[:, None])
    if do_normalize:
        signal/=(total/1000000)
    return signal


def heatplot_per_chrom(bedgraph, regions, y_coords, diffs, max_size):
    mids = (regions.ends+regions.starts)//2
    signals = bedgraph.get_slices(mids-max_size//2, mids+max_size//2, regions.directions)
    graph_diffs = (signal.scale_x(diffs.shape[1]).to_graph_diffs() for signal in signals)
    for y, signal in zip(y_coords, graph_diffs):
        signal.update_dense_array(diffs[y])
