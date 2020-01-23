import numpy as np

class GraphDiff:
    def __init__(self, start_value, indices, values):
        self._start_value = start_value
        self._indices = np.asanyarray(indices)
        self._values = np.asanyarray(values)

    def __eq__(self, other):
        t = self._start_value == other._start_value
        t &= np.all(self._indices == other._indices)
        return t and np.all(self._values==other._values)

    def __repr__(self):
        return "GD(%s, %s, %s)" % (self._start_value, self._indices, self._values)

    def update_dense_array(self, array):
        array[0] += self._start_value
        array[self._indices] += self._values

class BedGraph:
    def __init__(self, indices, values, size=None):
        assert indices[0] == 0, ("Indices does not start with 0", indices[:3])
        self._indices = np.asanyarray(indices)
        self._values = np.asanyarray(values)
        self._size = size

    def _getitem(self, index):
        idx = np.searchsorted(self._indices, index, side="right")-1
        if idx >= 0:
            return self._values[idx]
        else:
            return 0

    def reverse(self):
        assert self._size is not None
        indices = self._size-self._indices[::-1]
        values = self._values[::-1]
        assert indices[0] != 0
        indices = np.insert(indices[:-1], 0, 0)
        return self.__class__(indices, values, self._size)

    def get_slices(self, starts, ends, directions):
        start_idxs = np.searchsorted(self._indices, starts, side="right")
        end_idxs = np.searchsorted(self._indices, ends, side="left")
        start_values = self._values[start_idxs-1]
        for start_idx, end_idx, start_value, direction, start, end in zip(start_idxs, end_idxs, start_values, directions, starts, ends):
            indices = np.insert(self._indices[start_idx:end_idx]-start, 0, 0)
            values = np.insert(self._values[start_idx:end_idx], 0, start_value)
            new_obj = self.__class__(indices, values, end-start)
            if direction == -1:
                yield new_obj.reverse()
            else:
                yield new_obj
        

    def _getslice(self, slice_obj):
        assert slice_obj.step is None or slice_obj.step in (1, -1), slice_obj
        start_idx = np.searchsorted(self._indices, slice_obj.start, side="right")
        end_idx = np.searchsorted(self._indices, slice_obj.stop, side="left")
        start_value = self._values[start_idx-1]
        indices = np.insert(self._indices[start_idx:end_idx]-slice_obj.start, 0, 0)
        values = np.insert(self._values[start_idx:end_idx], 0, start_value)
        new_obj = self.__class__(indices, values, slice_obj.stop-slice_obj.start)
        if slice_obj.step == -1:
            return new_obj.reverse()
        return new_obj

    def get_start_idx(position):
        return np.searchsorted(self._indices, position, side="right")

    def get_end_idx(position):
        return np.searchsorted(self._indices, position, side="left")

    def to_graph_diffs(self):
        return GraphDiff(self._values[0],
                         self._indices[1:], np.diff(self._values))

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
        return "BG(%s, %s)" % (self._indices, self._values)
