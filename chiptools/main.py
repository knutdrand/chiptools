import sys, os
import cProfile
import pstats
import numpy as np
import trackhub
import matplotlib.pyplot as plt
import logging
import gzip
from collections import defaultdict

from .annotation import get_coding_offsets
from .metagene import metagene, coding_metagene, genome_metagene, show_transcripts
from .regions import Regions, get_holes
from .alignscores import alignscore
from .overlap import get_overlap, get_overlap_fraction
from .sizehist import get_hist, get_sizes
from .filterdup import filterdup
from .chainfile import parse_lines
from .bedIO import get_chroms, print_chroms, read_bedfile, read_bedgraphs, read_peakfile, read_fragments, print_regions, read_bedgraphs_fast
from .bedgraph import BedGraph
from .refSeqIO import parse_refseq_file
from .signalplot import signal_plot, signal_cumulative_hist
from .vplot import get_vplot, get_heatplot
from .regions import expand
from .genomebrowser import get_color, histone_track, coverage_track
from .rna_bam_to_bed import rna_bam_to_bed, get_splices
from .clipbed import clip_bed
from .fraglen import fraglen

log = logging

def do_tss_plot():
    do_normalize=True
    regions = read_bedfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs_fast(open(sys.argv[3]))
    N = 1000
    diffs = np.zeros(2*N)
    total = 0
    for chrom, bedgraph in bedgraphs:
        if chrom not in regions:
            continue
        print("Reading", chrom)
        chrom_regions = expand(regions[chrom], N, N)
        signal_plot(bedgraph, chrom_regions, diffs)
        total += bedgraph.sum()
    signal = np.cumsum(diffs)
    if do_normalize:
        signal/=(total/1000000)
    np.save(sys.argv[4], signal)
    if len(sys.argv)>5:
        plt.plot(np.arange(-N, N), signal)
        plt.savefig(sys.argv[5])

def do_signalhist():
    regions = read_bedfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs(sys.stdin)
    H = np.zeros(1000)
    for chrom, bedgraph in bedgraphs:
        if chrom not in regions or chrom=="chrM":
            continue
        assert "alt" not in chrom, chrom
        print("Reading", chrom)
        chrom_regions = regions[chrom]
        H += signal_cumulative_hist(bedgraph, chrom_regions)
    np.save(sys.argv[3], H)

def do_vplot():
    regions = read_bedfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs_fast(open(sys.argv[3]))
    signal = get_vplot(regions, bedgraphs)
    np.save(sys.argv[4], signal)
    if len(sys.argv)>5:
        plt.imshow(signal, cmap='gray_r') # , interpolation='nearest')
        plt.savefig(sys.argv[5])

def do_heatplot():
    regions = read_peakfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs_fast(open(sys.argv[3]))
    signal = get_heatplot(regions, bedgraphs)
    # singal = signal/(np.maximum(Ns, 1)[:, None])
    np.save(sys.argv[4], signal)
    if len(sys.argv) > 5:
        plt.figure(figsize=(15, 30))
        plt.imshow(signal, cmap='gray_r')
        plt.savefig(sys.argv[5])

def do_metagene():
    anno = parse_refseq_file(gzip.open(sys.argv[2], "rt"))
    bedgraphs = read_bedgraphs(sys.stdin)
    # anno.filter_coding()
    # anno.filter_largest()
    # offset_dict = {a.name: get_coding_offsets(a) for a in anno._annotations}
    # #utr_l_lens = sum(o[0] for o in offset_dict.values())/len(offset_dict)
    # cds_lens = sum(o[1]-o[0] for o in offset_dict.values())/len(offset_dict)
    # utr_r_lens = sum(sum(a.exonEnds)-sum(a.exonStarts)-offset_dict[a.name][1] if a.strand==1 else offset_dict[a.name][0]
    #                  for a in anno._annotations)/len(anno._annotations)
    # utr_l_lens = sum(sum(a.exonEnds)-sum(a.exonStarts)-offset_dict[a.name][1] if a.strand==-1 else offset_dict[a.name][0]
    #                  for a in anno._annotations)/len(anno._annotations)
    # 
    # print(utr_l_lens, cds_lens, utr_r_lens)
    # chroms = anno.to_indexed_regions()
    # 
    # N = 1000
    # diffs = [np.zeros(s) for s in (int(N*utr_l_lens/cds_lens), N, int(N*utr_r_lens/cds_lens))]
    # for chrom, bedgraph in bedgraphs:
    #     print("Reading", chrom)
    #     if chrom not in chroms or chrom=="chrM":
    #         continue
    #     indexed_regions = chroms[chrom]
    #     coding_metagene(bedgraph, indexed_regions, diffs, offset_dict)
    diffs = genome_metagene(anno, bedgraphs)
    t = 0
    for d in diffs:
        print(t)
        plt.plot(t+np.arange(d.size), np.cumsum(d))
        t += d.size
    plt.savefig(sys.argv[3])
    np.save(sys.argv[4], np.concatenate([np.cumsum(d) for d in diffs]))

def do_averageplot(gene=False):
    if gene:
        regions = read_bedfile(open(sys.argv[2]))
    else:
        regions = read_peakfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs_fast(open(sys.argv[3]))
    N=2000
    diffs = np.zeros(2*N)
    total = 0
    for chrom, bedgraph in bedgraphs:
        print("Reading", chrom)
        if chrom not in regions or chrom=="chrM":
            continue
        chrom_regions = regions[chrom]
        sizes = chrom_regions.ends-chrom_regions.starts
        new_regions = Regions(chrom_regions.starts-sizes/2, chrom_regions.ends+sizes/2)
        signal_plot(bedgraph, new_regions, diffs, scale_to=True)
        total+=regions[chrom].starts.size
    signal = np.cumsum(diffs)/total
    np.save(sys.argv[4], signal)
    if len(sys.argv)>5:
        plt.plot(signal)
        plt.savefig(sys.argv[5])

def do_regionmeans(args):
    regions = read_peakfile(open(sys.argv[2]))
    if "-holes" in args:
        regions = {chrom: get_holes(r) for chrom, r in regions.items()}
    bedgraphs = read_bedgraphs(sys.stdin)
    for chrom, bedgraph in bedgraphs:
        chrom_regions = regions[chrom]
        means = (s.mean() for s in bedgraph.get_slices(chrom_regions.starts, chrom_regions.ends, chrom_regions.directions))
        for start, end, mean in zip(chrom_regions.starts, chrom_regions.ends, means):
            print(f"{chrom}\t{start}\t{end}\t{mean}")

def get_input_lines(arg):
    if arg=="-":
        return sys.stdin
    return open(arg)

def main():
    if sys.argv[1] == "overlap":
        data_a = get_chroms(open(sys.argv[2]))
        data_b = get_chroms(open(sys.argv[3]))
        overlap = {chrom: np.hstack((data_a[chrom], get_overlap(data_a[chrom], data_b[chrom])))
                   for chrom in data_a if chrom in data_b}
        print_chroms(overlap)
    elif sys.argv[1] == "overlap_fraction":
        data_a = get_chroms(open(sys.argv[2]))
        data_b = get_chroms(open(sys.argv[3]))
        overlap = {chrom: np.hstack((data_a[chrom], get_overlap_fraction(data_a[chrom], data_b[chrom])))
                   for chrom in data_a if chrom in data_b}
        print_chroms(overlap)

    elif sys.argv[1] == "tssplot":
        do_tss_plot()

    elif sys.argv[1] in ("geneaverageplot", "averageplot"):
        do_averageplot(sys.argv[1]=="geneaverageplot")

    elif sys.argv[1] == "vplot":
        cProfile.runctx("do_vplot()", globals(), locals(), "profiling")
        stats = pstats.Stats("profiling")
        stats.sort_stats("cumulative")
        stats.print_stats()

    elif sys.argv[1] == "signalhist":
        do_signalhist()

    elif sys.argv[1] == "heatplot":
        do_heatplot()

    elif sys.argv[1] == "hist":
        sizes = [float(l.strip()) for l in get_input_lines(sys.argv[2])]
        plt.hist(sizes, range=(0, float(sys.argv[3])), bins=int(sys.argv[4]))
        if len(sys.argv)>5:
            plt.savefig(sys.argv[5])
        else:
            plt.show()

    elif sys.argv[1] == "sizehist":
        bins = get_hist(get_sizes(get_input_lines(sys.argv[2])),
                        int(sys.argv[3]), int(sys.argv[4]))
        plt.xlabel("size")
        plt.ylabel("count")
        plt.savefig(sys.argv[5])
        if len(sys.argv) > 6:
            np.save(sys.argv[6], bins)

    elif sys.argv[1] == "filterdup":
        chroms = read_fragments(open(sys.argv[2]))
        for chrom, data in chroms:
            print_regions(chrom, filterdup(data))

    elif sys.argv[1] == "chain":
        input_lines = get_input_lines(sys.argv[2])
        output_lines = parse_lines(input_lines)
        for output_line in output_lines:
            print("\t".join(str(c) for c in output_line))

    elif sys.argv[1] == "trackdb":
        names = [os.path.basename(n).split(".")[0] for n in sys.argv[2:]]
        f = histone_track
        if names[0] == "single":
            f = coverage_track
            names = names[1:]
        colors = [get_color(i, len(names)) for i in range(len(names))]
        hub, genomes_file, genome, trackdb = trackhub.default_hub(
            hub_name="testing",
            genome="hg38",
            email="knutdrand@gmail.com")
        #db = trackhub.TrackDb()
        trackdb.add_tracks([f(*pair)
                       for pair in zip(names, colors)])
        print(trackdb)

    elif sys.argv[1] == "plot":
        names = sys.argv[2:-1]
        lines = [plt.plot(np.load(name))[0] for name in names]
        plt.legend(lines, [name.split(".")[0] for name in names])
        plt.savefig(sys.argv[-1])

    elif sys.argv[1] == "rnabam2bed":
        lines = (f"{chrom}\t{start}\t{end}\t.\t.\t{direction}\t{is_end}" 
                 for chrom, start, end, direction, is_end in rna_bam_to_bed(sys.stdin))
        for line in lines:
            print(line)

    elif sys.argv[1] == "getsplices":
        lines = (f"{chrom}\t{start}\t{end}\t.\t.\t{direction}" 
                 for chrom, start, end, direction in get_splices(sys.stdin))
        for line in lines:
            print(line)

    elif sys.argv[1] == "clipbed":
        parts = [line.strip().split() for line in open(sys.argv[3])]
        chrom_sizes = {p[0]: int(p[1]) for p in parts if p}
        for line in clip_bed(open(sys.argv[2]), chrom_sizes):
            sys.stdout.write(line)

    elif sys.argv[1] == "alignscores":
        read_length = int(sys.argv[2])
        align_scores = alignscore(sys.stdin, read_length)
        np.save(sys.argv[3], align_scores)

    elif sys.argv[1] == "fraglen":
        fraglen(sys.argv[2])

    elif sys.argv[1] == "pairbed":
        from .pair_intervals import pair_lines
        for new_entry in pair_lines(sys.stdin):
            print(new_entry)

    elif sys.argv[1] == "regionmeans":
        do_regionmeans(sys.argv)

    elif sys.argv[1] == "metagene":
        do_metagene()
