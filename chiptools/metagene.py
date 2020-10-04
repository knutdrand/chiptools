from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from .bedgraph import BedGraph
from .annotation import get_coding_offsets


def metagene(bedgraph, indexed_regions, diffs):
    regions = indexed_regions.regions
    signals = bedgraph.get_slices_normal(regions.starts, regions.ends, regions.directions)
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
    signals = bedgraph.get_slices_normal(regions.starts, regions.ends, regions.directions)
    exons = defaultdict(list)
    directions = defaultdict(list)
    for idx, direction in zip(indexed_regions.indexes, regions.directions):
        directions[idx].append(direction)
    for index, signal in zip(indexed_regions.indexes, signals):
        exons[index].append(signal)
    joined_exons = {}
    exon_offsets = {}
    for idx, bgs in exons.items():
        dirs = directions[idx]
        assert all(d==dirs[0] for d in dirs), (idx, dirs)
        if dirs[0] == -1:
            bgs = bgs[::-1]
        joined_exons[indexed_regions.names[idx]] = BedGraph.concatenate(bgs)
        exon_offsets[indexed_regions.names[idx]] = np.cumsum([bg._size for bg in bgs])
    counts = np.array([0, 0, 0])
    for name, e in joined_exons.items():
        tmp = np.zeros(e._size)
        e.to_graph_diffs().update_dense_array(tmp)
        tmp = np.cumsum(tmp)
        if np.max(tmp)>=4:
            plt.plot(tmp)
            plt.title(name)
            plt.vlines(exon_offsets[name], ymin=0, ymax=tmp.max())
            plt.show()
        offsets = offset_dict[name]
        if offsets[0]>0:
            counts[0]+=1
            e[:offsets[0]].scale_x(diffs[0].size).to_graph_diffs().update_dense_array(diffs[0])
        counts[1] += 1
        e[offsets[0]:offsets[1]].scale_x(diffs[1].size).to_graph_diffs().update_dense_array(diffs[1])
        if offsets[1] < e._size:
            counts[2]+=1
            e[offsets[1]:].scale_x(diffs[2].size).to_graph_diffs().update_dense_array(diffs[2])
    return counts

def genome_metagene(anno, bedgraphs):
    names = ["Mrpl15", "Lypla1", "Tcea1", "Cops5"]
    show_many_transcripts(anno, next(bg for c, bg in bedgraphs if c=="chr1"), names)
    exit()
    anno.filter_coding()
    anno.filter_largest()
    offset_dict = {a.name: get_coding_offsets(a) for a in anno._annotations}

    cds_lens = sum(o[1]-o[0] for o in offset_dict.values())/len(offset_dict)
    utr_r_lens = sum(sum(a.exonEnds)-sum(a.exonStarts)-offset_dict[a.name][1] if a.strand==1 else offset_dict[a.name][0]
                     for a in anno._annotations)/len(anno._annotations)
    utr_l_lens = sum(sum(a.exonEnds)-sum(a.exonStarts)-offset_dict[a.name][1] if a.strand==-1 else offset_dict[a.name][0]
                     for a in anno._annotations)/len(anno._annotations)


    chroms = anno.to_indexed_regions()
    N = 1000
    diffs = [np.zeros(s) for s in (int(N*utr_l_lens/cds_lens), N, int(N*utr_r_lens/cds_lens))]
    counts = np.array([0, 0, 0])
    for chrom, bedgraph in bedgraphs:
        print("Reading", chrom)
        if chrom not in chroms or chrom=="chrM":
            continue
        indexed_regions = chroms[chrom]
        counts += coding_metagene(bedgraph, indexed_regions, diffs, offset_dict)

    return [d/count for d, count in zip(diffs, counts)]
    # t = 0
    # 
    # for d in diffs:
    #     plt.plot(t+np.arange(d.size), np.cumsum(d))
    #     t += d.size
    #     
    # plt.savefig(sys.argv[3])

def show_many_transcripts(annotation, bedgraph, names):
    for name in names:
        show_transcripts(annotation, bedgraph, name)
    
def show_transcripts(annotation, bedgraph, name):
    print("Showing transcripts for ", name)
    variants=annotation.get_gene(name)
    indexed_regions = variants.to_indexed_regions()
    chrom = list(indexed_regions.keys())[0]
    indexed_regions = indexed_regions[chrom]
    regions = indexed_regions.regions
    signals = bedgraph.get_slices_normal(regions.starts, regions.ends, regions.directions)
    exons = defaultdict(list)
    directions = defaultdict(list)
    for idx, direction in zip(indexed_regions.indexes, regions.directions):
        directions[idx].append(direction)
    for index, signal in zip(indexed_regions.indexes, signals):
        exons[index].append(signal)
    joined_exons = {}
    exon_offsets = {}
    for idx, bgs in exons.items():
        dirs = directions[idx]
        assert all(d==dirs[0] for d in dirs), (idx, dirs)
        if dirs[0] == -1:
            bgs = bgs[::-1]
        joined_exons[idx] = BedGraph.concatenate(bgs)
        exon_offsets[idx] = np.cumsum([bg._size for bg in bgs])
    fig, axs = plt.subplots(len(joined_exons))
    xmax = max(e._size for e in joined_exons.values())
    for idx, e in joined_exons.items():
        tmp = np.zeros(e._size)
        e.to_graph_diffs().update_dense_array(tmp)
        tmp = np.cumsum(tmp)
        axs[idx].plot(tmp)
        axs[idx].vlines(exon_offsets[idx], ymin=0, ymax=tmp.max())
        axs[idx].set_xlim(0, xmax)
    plt.title(name)
    plt.show()
