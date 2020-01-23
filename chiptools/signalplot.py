import numpy as np

def _signal_plot(bedgraph, regions, size):
    diffs = np.zeros(size)
    for region in regions:
        signal = bedgraph[region.start:region.end:region.direction]
        signal.to_graph_diffs().update_dense_array(diffs)
    return np.cumsum(diffs)

def signal_plot(bedgraph, regions, size):
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    diffs = np.zeros(size)
    for signal in signals:
        signal.to_graph_diffs().update_dense_array(diffs)
    return np.cumsum(diffs)
                                        
