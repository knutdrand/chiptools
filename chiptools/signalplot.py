import numpy as np

def scale(indices, from_size, to_size):
    return np.rint(indices/from_size*to_size)

def signal_plot(bedgraph, regions, size, scale_to=False):
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    diffs = np.zeros(size)
    
    for signal in signals:
        graph_diffs = signal.to_graph_diffs()
        if scale_to:
            graph_diffs = graph_diffs.scale_x(size)
        graph_diffs.update_dense_array(diffs)

    return np.cumsum(diffs)
                                        
