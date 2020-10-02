import numpy as np

def signal_plot(bedgraph, regions, diffs, scale_to=False):
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    for signal in signals:
        if scale_to:
            signal = signal.scale_x(diffs.size)
        graph_diffs = signal.to_graph_diffs()
        graph_diffs.update_dense_array(diffs)

def signal_cumulative_hist(bedgraph, regions):
    H = np.zeros(1000)
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    for signal in signals:
        h = signal.hist()
        H[:h.size]+=signal.hist()
    return np.cumsum(H)
    
