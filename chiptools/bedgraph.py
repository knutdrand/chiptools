import numpy as np
import logging
from collections import Counter
from itertools import chain

from .regions import Regions

class GraphDiff:
    def __init__(self, start_value, indices, values, size=None):
        self._start_value = start_value
        self._indices = np.asanyarray(indices)
        assert np.all(np.diff(self._indices)>0), self._indices
        assert np.all(self._indices>=0), self._indices
        self._values = np.asanyarray(values)
        self._size = size

    def __eq__(self, other):
        t = self._start_value == other._start_value
        t &= np.all(self._indices == other._indices)
        return t and np.all(self._values==other._values)

    def __repr__(self):
        return "GD(%s, %s, %s)" % (self._start_value, self._indices, self._values)

    def update_dense_array(self, array):
        self.assert_positive()
        array[0] += self._start_value
        assert self._indices.size==0 or array.size>=np.max(self._indices), (array.size, np.max(self._indices))
        array[self._indices] += self._values

    def assert_positive(self):
        values = self._start_value + np.cumsum(self._values)
        assert np.all(values >= -1e-10), (values[values < -1e-10], values.dtype)

class BedGraph:
    def __init__(self, indices, values, size=None, strict=True):
        assert (not strict) or indices[0] == 0, ("Indices does not start with 0", indices[:3])
        self._indices = np.asanyarray(indices)
        assert np.all(np.diff(self._indices)>0), indices
        if size is not None:
            assert np.all(self._indices<size), (self._indices, size)

        self._values = np.asanyarray(values)
        assert np.all(self._values>=0), self._values[self._values<0]
        self._size = size
        assert size is None or size>0

    def __iter__(self):
        pairs = zip(self._indices, chain(self._indices[1:], [self._size]))
        return zip(pairs, self._values)

    def _getitem(self, index):
        idx = np.searchsorted(self._indices, index, side="right")-1
        if idx >= 0:
            return self._values[idx]
        else:
            return 0

    def sum(self):
        if self._size is not None:
            sizes = np.diff(self._indices, append=self._size)
            return np.sum(sizes*self._values)
        assert self._values[-1] == 0, self._values[-1]
        return np.diff(self._indices*self._values[:-1])

    def mean(self):
        assert self._indices[0]==0
        return self.sum()/self._size

    def hist(self):
        h = np.zeros(int(np.max(self._values))+1)
        for (start, end), value in self:
            h[int(value)]+=(end-start)
        return h

    def reverse(self):
        assert self._size is not None
        indices = self._size-self._indices[::-1]
        values = self._values[::-1]
        assert indices[0] != 0
        indices = np.insert(indices[:-1], 0, 0)
        return self.__class__(indices, values, self._size)

    def get_slices(self, starts, ends, directions):
        assert np.all(ends>starts)
        start_idxs = np.searchsorted(self._indices, starts, side="right")
        end_idxs = np.searchsorted(self._indices, ends, side="left")
        start_values = self._values[start_idxs-1]
        for start_idx, end_idx, start_value, direction, start, end in zip(start_idxs, end_idxs, start_values, directions, starts, ends):
            indices = np.insert(self._indices[start_idx:end_idx]-start, 0, 0)
            values = np.insert(self._values[start_idx:end_idx], 0, start_value)
            new_obj = self.__class__(indices, values, end-start)
            assert direction in (-1, 1)
            if direction == -1:
                yield new_obj.reverse()
            else:
                yield new_obj

    def threshold(self, value):
        over = self._values >= value
        changes = np.flatnonzero(over[1:] != over[:-1])
        starts = self._indices[changes[::2]+1]
        ends = self._indices[changes[1::2]+1]
        return Regions(starts, ends)

    def _getslice(self, slice_obj):
        assert slice_obj.step is None or slice_obj.step in (1, -1), slice_obj
        start = slice_obj.start or 0
        start_idx = np.searchsorted(self._indices, start, side="right")
        stop = slice_obj.stop
        if stop is None:
            assert self._size is not None
            stop = self._size
        assert self._size is None or stop<=self._size, (slice_obj, self)
        end_idx = np.searchsorted(self._indices, stop, side="left")
        start_value = self._values[start_idx-1]
        indices = np.insert(self._indices[start_idx:end_idx]-start, 0, 0)
        values = np.insert(self._values[start_idx:end_idx], 0, start_value)
        new_obj = self.__class__(indices, values, stop-start)
        if slice_obj.step == -1:
            return new_obj.reverse()
        return new_obj

    def get_start_idx(position):
        return np.searchsorted(self._indices, position, side="right")

    def get_end_idx(position):
        return np.searchsorted(self._indices, position, side="left")

    def scale_x(self, size):
        new_indices = (self._indices*size/self._size).astype("int")
        ds = np.concatenate((np.diff(new_indices)>0, [True]))
        return BedGraph(new_indices[ds],
                        self._values[ds], size)

    @classmethod
    def concatenate(cls, bedgraphs):
        offsets = np.cumsum([0] + [bg._size for bg in bedgraphs[:-1]])
        indices = np.concatenate([bg._indices+offset for (bg, offset) in zip(bedgraphs, offsets)])
        assert np.all(np.diff(indices)>0), (bedgraphs, indices, offsets)
        values = np.concatenate([bg._values for bg in bedgraphs])
        size = sum(bg._size for bg in bedgraphs)
        assert np.all(indices<size), (size, offsets, indices)
        return cls(indices, values, size)

    def update_dense_diffs(self, diffs):
        diffs[0]+=self._values[0]
        diffs[self._indices[1:]]+= np.diff(self._values)

    def to_graph_diffs(self):
        gd = GraphDiff(self._values[0],
                       self._indices[1:], np.diff(self._values), self._size)
        gd.assert_positive()
        return gd

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._getslice(index)
        if isinstance(index, list):
            return [self._getitem(i) for i in index]
        return self._getitem(index)

    def __eq__(self, other):
        t = np.all(self._indices==other._indices)
        return t and np.all(self._values==other._values)
    
    def __repr__(self):
        return "BG(%s, %s, %s)" % (self._indices, self._values, self._size)
