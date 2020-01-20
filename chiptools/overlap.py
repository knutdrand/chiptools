import numpy as np

def get_overlap(x, ref):
    x = np.ravel(x)
    cum_sizes = np.insert(np.cumsum(ref[:, 1]-ref[:, 0]), 0, 0)
    idxs = np.searchsorted(np.ravel(ref), x)
    count_until = cum_sizes[idxs//2]
    mask = idxs %2 == 1
    extra = x[mask]-np.ravel(ref)[idxs[mask]-1]
    count_until[mask] += extra
    return np.diff(count_until.reshape((-1, 2)))
    # count_until[mask] += x[mask]-
    
def get_overlap_fraction(x, ref):
    overlap = get_overlap(x, ref)
    sizes = np.diff(x)
    return overlap/sizes

def _get_overlap(coords_a, coords_b):
    thr = max(coords_a[-1, -1], coords_b[-1,-1])+10
    coords_b = np.vstack((coords_b, (thr, thr+10)))
    sizes_b = coords_b[:, 1]-coords_b[:, 0]
    cum_sizes = np.insert(np.cumsum(sizes_b), 0, 0)
    start_idxs = np.searchsorted(coords_b[:, 1], coords_a[:, 0])
    end_idxs = np.searchsorted(coords_b[:, 0], coords_a[:, 1])
    start_value = cum_sizes[start_idxs]
    all_overlap = cum_sizes[end_idxs]-cum_sizes[start_idxs]
    start_overlap = np.maximum(coords_a[:, 0]-coords_b[start_idxs-1, 0], 0)
    end_overlap = np.maximum(coords_b[end_idxs-1, 1]-coords_a[:, 1], 0)
    overlap = all_overlap-start_overlap-end_overlap
    # assert np.all(overlap>=0), (coords_a[overlap<0], start_value[overlap<0])
    return (all_overlap-start_overlap-end_overlap)[:, None]
