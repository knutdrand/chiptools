from collections import namedtuple, defaultdict
from itertools import groupby

import numpy as np

from .regions import Regions

Anno = namedtuple("Anno", ["chrom", "strand", "txStart", "txEnd",
                           "cdsStart", "cdsEnd", "exonStarts", "exonEnds", "name", "isoform"])

def get_coding_offsets(anno):
    #assert anno.txStart<anno.cdsStart<anno.cdsEnd<anno.txEnd, anno
    offsets = np.cumsum([0]+[end-start for start, end in zip(anno.exonStarts, anno.exonEnds)])
    idxs = np.searchsorted(anno.exonStarts, [anno.cdsStart, anno.cdsEnd], side="right")-1

    for idx, pos in zip(idxs, [anno.cdsStart, anno.cdsEnd]):
        assert idx<len(anno.exonStarts), (idx, pos, anno)
        assert anno.exonStarts[idx]<=pos, (anno, idx, pos)
        assert idx==len(anno.exonStarts)-1 or pos < anno.exonStarts[idx+1], (anno, idx, pos)
    return tuple(offsets[i]+pos-anno.exonStarts[i] for i, pos in zip(idxs, [anno.cdsStart, anno.cdsEnd]))

class Annotations:
    def __init__(self, annotations):
        self._annotations = list(sorted(annotations, key=lambda x: x.name))

    def filter_largest(self):
        groups = groupby(self._annotations, key=lambda x: x.name)
        self._annotations =  [max(group, key=lambda x: sum(end-start for start, end in zip(x.exonStarts, x.exonEnds))) for _, group in groups]

    def filter_coding(self):
        self._annotations = [anno for anno in self._annotations if anno.cdsStart!=anno.cdsEnd]
        
    def position_sort(self):
        self._annotations.sort(key=lambda x: (x.chrom, x.txStart))

    def get_gene(self, name):
        return Annotations(a for a in self._annotations if a.name==name)

    def to_indexed_regions(self):
        self.position_sort()
        groups = groupby(self._annotations, lambda x: x.chrom)
        indexed_regions = {}
        for chrom, group in groups:
            annotations = list(group)
            starts = np.array([start for anno in annotations for start in anno.exonStarts])
            args = np.argsort(starts)
            ends = np.array([end for anno in annotations for end in anno.exonEnds])[args]
            indexes =  np.array([i for i, anno in enumerate(annotations) for end in anno.exonEnds])[args]
            directions =  np.array([anno.strand for anno in annotations for end in anno.exonEnds])[args]
            names = [anno.name for anno in annotations]
            regions = Regions(starts[args], ends, directions)
            indexed_regions[chrom] = IndexedRegions(names, indexes, regions)
        return indexed_regions

    def get_splice_site_dict(self):
        splice_site_dict = defaultdict(list)
        for anno in self._annotations:
            for end, start in zip(anno.exonEnds[:-1], anno.exonStarts[1:]):
                splice_site_dict[(anno.chrom, end, start)].append(anno.isoform)
        return splice_site_dict
        

class IndexedRegions:
    def __init__(self, names, indexes, regions):
        self.names = names
        self.indexes = indexes
        self.regions = regions

    
