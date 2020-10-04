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
        # assert (not strict) or indices[0] == 0, ("Indices does not start with 0", indices[:3])
        self._indices = np.asanyarray(indices)
        assert np.issubdtype(self._indices.dtype, np.integer), self._indices
        # assert np.all(np.diff(self._indices)>0), indices
        if size is not None:
            assert np.all(self._indices<size), (self._indices, size)

        self._values = np.asanyarray(values)
        # assert np.all(self._values>=0), self._values[self._values<0]
        if size is not None:
            self._size = int(size)
        # assert size is None or size>0

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
        assert indices[0] != 0
        indices = np.insert(indices[:-1], 0, 0)
        values = self._values[::-1]
        return self.__class__(indices, values, self._size)

    def _get_slice_indexes(self, start_idxs, end_idxs, directions, offsets):
        index_diffs = np.zeros(offsets[-1], dtype="int")
        index_diffs[offsets[1:-1]] = np.diff(directions)
        index_diffs[0] = directions[0]
        index_diffs = np.cumsum(index_diffs)
        left_idxs = np.where(directions==1, start_idxs-1, end_idxs)#?
        right_idxs = np.where(directions==1, end_idxs-1, start_idxs)#?
        index_diffs[offsets[1:-1]] = left_idxs[1:]-right_idxs[:-1]# np.diff(left_idxs)
        index_diffs[0] = left_idxs[0]
        return np.cumsum(index_diffs)

    def get_slices(self, starts, ends, directions):
        start_idxs = np.searchsorted(self._indices, starts, side="right")
        end_idxs = np.searchsorted(self._indices, ends, side="left")
        offsets = np.insert(np.cumsum(end_idxs-start_idxs+1), 0, 0)
        slice_indexes=self._get_slice_indexes(start_idxs, end_idxs, directions, offsets)
        print(starts[:10], ends[:10], directions[:10])
        print(start_idxs[:10], end_idxs[:10], directions[:10])
        for f, e in zip(offsets[:10], offsets[1:11]):
            print(slice_indexes[f:e])
        assert np.all(slice_indexes>=0)
        assert np.all(slice_indexes<=self._values.size), (slice_indexes[slice_indexes>=self._values.size], self._values.size)
        values = np.insert(self._values, self._values.size, 0)[slice_indexes]
        start_diffs = np.zeros_like(slice_indexes)
        start_diffs[offsets[1:-1]] = np.diff(starts)
        start_diffs[0] = starts[0]
        indices = np.insert(self._indices, self._indices.size, self._size)[slice_indexes]-np.cumsum(start_diffs)
        sizes = ends-starts
        indices[offsets[:-1]] = np.where(directions==1, 0, sizes-1)
        for f, e in zip(offsets[:10], offsets[1:11]):
            print(indices[f:e])
        assert np.all(np.where(directions==1, 0, sizes-1)>=0)
        assert np.all(indices>=0)# ,  (indices, np.flatnonzero(indices<0), offsets)
        return BedGraphArray(indices, values, ends-starts, offsets)

    def get_slices_normal(self, starts, ends, directions):
        start_idxs = np.searchsorted(self._indices, starts, side="right")
        end_idxs = np.searchsorted(self._indices, ends, side="left")
        start_values = self._values[start_idxs-1]
        offsets = np.cumsum(end_idxs-start_idxs+1)
        all_idxs = np.empty(offsets[-1], dtype="int")
        all_values = np.empty(offsets[-1], dtype="float")
        all_idxs[offsets[:-1]] = 0
        all_idxs[0] = 0
        
        all_values[np.where(directions[1:]==1, offsets[:-1],offsets[1:]-1)] = start_values[1:]
        if directions[0]==1:
            all_values[0] = start_values[0]
        else:
            all_values[offsets[0]-1] = start_values[0]
        for start_idx, end_idx, direction, start, end, offset in zip(start_idxs, end_idxs, directions, starts, ends, offsets):
            S = end_idx-start_idx
            if direction == 1:
                all_idxs[offset-S:offset] = self._indices[start_idx:end_idx]-start
                all_values[offset-S:offset] = self._values[start_idx:end_idx]
            else:
                all_idxs[offset-S:offset] = end-self._indices[start_idx:end_idx][::-1]
                all_values[offset-S-1:offset-1] = self._values[start_idx:end_idx][::-1]
            yield self.__class__(all_idxs[offset-S-1:offset],
                                 all_values[offset-S-1:offset],
                                 end-start)

    def get_slices_slow(self, starts, ends, directions):
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
        assert np.issubdtype(self._indices.dtype, np.integer), self._indices
        new_indices = (self._indices*size//self._size)
        assert np.issubdtype(new_indices.dtype, np.integer), (self._indices, size, self._size)
        ds = np.concatenate((np.diff(new_indices)>0, [True]))
        return BedGraph(new_indices[ds], self._values[ds], size)
                        

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
        assert np.issubdtype(self._indices.dtype, np.integer), self._indices
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

class BedGraphArray:
    def __init__(self, indices, values, sizes, offsets):
        self._indices = indices
        assert np.all(self._indices>=0), self._indices.dtype
        self._values = values
        self._sizes = sizes
        assert offsets[-1]==self._indices.size, (offsets[-1], self._indices.size)
        self._offsets = offsets

    def scale_x(self, size):
        assert np.issubdtype(self._indices.dtype, np.integer), self._indices
        sizes = np.zeros_like(self._indices)
        sizes[0] = self._sizes[0]
        sizes[self._offsets[1:-1]] = np.diff(self._sizes)
        new_indices = (self._indices*size//np.cumsum(sizes))
        assert np.issubdtype(new_indices.dtype, np.integer), (self._indices, size, self._size)
        mask = np.concatenate((np.diff(new_indices)>0, [True]))
        mask[self._offsets[1:]-1] = True
        counts = np.cumsum(mask)
        new_offsets = np.insert(counts[self._offsets[1:]-1], 0, 0)
        return BedGraphArray(new_indices[mask], self._values[mask], size*np.ones_like(self._sizes), new_offsets)

    def _broadcast(self, values):
        assert values.size == self._offsets.size-1, (values.size, self._offsets.size)
        broadcasted = np.zeros_like(self._indices, dtype=values.dtype)
        broadcasted[self._offsets[1:-1]] = np.diff(values)
        broadcasted[0] = values[0]
        return np.cumsum(broadcasted)

    def update_dense_diffs(self, diffs, rows):
        ncols = diffs.shape[1]
        all_rows = self._broadcast(rows)
        set_a = set(rows.flatten())
        set_b = set(all_rows.flatten())
        assert set_a==set_b, (set_a-set_a, set_b-set_a)
        assert np.all(self._indices < ncols)
        assert np.all(self._indices >= 0), self._indices
        composite_indexes = all_rows*ncols + self._indices
        assert np.all(composite_indexes//ncols==all_rows), (composite_indexes//ncols, all_rows)
        set_b = set((composite_indexes//ncols).flatten())
        assert set_a==set_b, (set_a-set_a, set_b-set_a)
        args = np.argsort(composite_indexes, kind="mergesort")
        indices = composite_indexes[args]
        index_changes = np.insert(indices[:-1] != indices[1:], indices.size-1, True)
        value_diffs = np.insert(np.diff(self._values), 0, self._values[0])
        value_diffs[self._offsets[:-1]] = self._values[self._offsets[:-1]]
        value_diffs = value_diffs[args]
        totals = np.insert(np.cumsum(value_diffs)[index_changes], 0, 0)
        total_diffs = np.diff(totals)
        used_indices = indices[index_changes]
        set_a = set(rows.flatten())
        set_b = set((used_indices//ncols).flatten())
        assert set_a==set_b, (set_a-set_a, set_b-set_a)
        diffs[used_indices//ncols, used_indices % ncols] += total_diffs


    def __getitem__(self, idx):
        assert idx<self._offsets.size-1
        s = slice(self._offsets[idx], self._offsets[idx+1])
        return BedGraph(self._indices[s], self._values[s], self._sizes[idx])

    def __iter__(self):
        return (self[i] for i in range(self._sizes.size))
