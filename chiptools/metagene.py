from .bedgraph import BedGraph
from collections import defaultdict

def metagene(bedgraph, indexed_regions, diffs):
    regions = indexed_regions.regions
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    exons = defaultdict(list)
    directions = defaultdict(list)
    for idx, direction in zip(indexed_regions.indexes, regions.directions):
        directions[idx].append(direction)
    for index, signal in zip(indexed_regions.indexes, signals):
        exons[index].append(signal)
    joined_exons = []
    for idx, bgs in exons.items():
        dirs = directions[idx]
        assert all(d==dirs[0] for d in dirs), (idx, dirs)
        if dirs[0] == -1:
            bgs = bgs[::-1]
        joined_exons.append(BedGraph.concatenate(bgs))
    for e in joined_exons:
        e.scale_x(diffs.size).to_graph_diffs().update_dense_array(diffs)

def coding_metagene(bedgraph, indexed_regions, diffs, offset_dict):
    regions = indexed_regions.regions
    signals = bedgraph.get_slices(regions.starts, regions.ends, regions.directions)
    exons = defaultdict(list)
    directions = defaultdict(list)
    for idx, direction in zip(indexed_regions.indexes, regions.directions):
        directions[idx].append(direction)
    for index, signal in zip(indexed_regions.indexes, signals):
        exons[index].append(signal)
    joined_exons = {}
    for idx, bgs in exons.items():
        dirs = directions[idx]
        assert all(d==dirs[0] for d in dirs), (idx, dirs)
        if dirs[0] == -1:
            bgs = bgs[::-1]
        joined_exons[indexed_regions.names[idx]] = BedGraph.concatenate(bgs)
    for name, e in joined_exons.items():
        offsets = offset_dict[name]
        if offsets[0]>0:
            e[:offsets[0]].scale_x(diffs[0].size).to_graph_diffs().update_dense_array(diffs[0])
        e[offsets[0]:offsets[1]].scale_x(diffs[1].size).to_graph_diffs().update_dense_array(diffs[1])
        if offsets[1] < e._size:
            e[offsets[1]:].scale_x(diffs[2].size).to_graph_diffs().update_dense_array(diffs[2])

