import numpy as np
from .regions import Regions

def filterdup(regions):
    ind = np.lexsort((regions.ends, regions.starts))
    starts = regions.starts
    ends = regions.ends[ind]
    mask = starts[1:] == starts[:-1]
    mask &= ends[1:] == ends[:-1]
    return Regions(
        np.append(starts[:-1][~mask], starts[-1]),
        np.append(ends[:-1][~mask], ends[-1]))
