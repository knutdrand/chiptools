import numpy as np

def signal_plot(bedgraph, regions, size, scale_to=False):
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    diffs = np.zeros(size)
    
    for signal in signals:
        graph_diffs = signal.to_graph_diffs()
        if scale_to:
            graph_diffs = graph_diffs.scale_x(size)
        graph_diffs.update_dense_array(diffs)

    return np.cumsum(diffs)
