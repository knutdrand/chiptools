from collections import namedtuple
import numpy as np
Region = namedtuple("Region", ["start", "end", "direction"])

class Regions:
    def __init__(self, starts, ends, directions=1):
        self.starts = starts
        self.ends = ends
        self.directions = directions

    def __iter__(self):
        return (Region(*t) for t in zip(self.starts, self.ends, self.directions))

def expand(regions, upstream, downstream):
    centers = np.where(regions.directions==1, regions.starts, regions.ends-1)
    starts = centers-np.where(regions.directions==1, upstream, downstream)
    ends = centers+np.where(regions.directions==1, downstream, upstream)
    return Regions(starts, ends, regions.directions)
