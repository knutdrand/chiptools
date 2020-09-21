import sys, os
import numpy as np
import trackhub
import matplotlib.pyplot as plt
from .regions import Regions, get_holes
from .alignscores import alignscore
from .overlap import get_overlap, get_overlap_fraction
from .sizehist import get_hist, get_sizes
from .filterdup import filterdup
from .chainfile import parse_lines
from .bedIO import get_chroms, print_chroms, read_bedfile, read_bedgraphs, read_peakfile, read_fragments, print_regions
from .signalplot import signal_plot, signal_cumulative_hist
from .vplot import vplot, get_heatplot
from .regions import expand
from .genomebrowser import get_color, histone_track, coverage_track
from .rna_bam_to_bed import rna_bam_to_bed
from .clipbed import clip_bed
from .fraglen import fraglen

def do_tss_plot():
    regions = read_bedfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs(sys.stdin)
    N = 1000
    signal = np.zeros(2*N)
    for chrom, bedgraph in bedgraphs:
        if chrom not in regions or chrom=="chrM":
            continue
        print("Reading", chrom)
        chrom_regions = expand(regions[chrom], N, N)
        signal += signal_plot(bedgraph, chrom_regions, 2*N)
    np.save(sys.argv[3], signal)
    if len(sys.argv)>4:
        plt.plot(np.arange(-N, N), signal)
        plt.savefig(sys.argv[4])

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
    bedgraphs = read_bedgraphs(sys.stdin)
    size_x=10000
    size_y=5000
    max_size=10000
    signal = np.zeros((size_y, size_x))
    Ns = np.zeros(size_y)
    for chrom, bedgraph in bedgraphs:
        if chrom not in regions or chrom=="chrM":
            continue
        assert "alt" not in chrom, chrom
        print("Reading", chrom)
        chrom_regions = regions[chrom]
        s, Ns_tmp = vplot(bedgraph, chrom_regions, size_x, size_y, max_size)
        signal += s
        Ns += Ns_tmp
    singal = signal/(np.maximum(Ns, 1)[:, None])
    np.save(sys.argv[3], signal)
    if len(sys.argv)>4:
        plt.imshow(signal, cmap='hot', interpolation='nearest')
        plt.savefig(sys.argv[4])

def do_heatplot():
    regions = read_peakfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs(sys.stdin)
    signal = get_heatplot(regions, bedgraphs)
    # singal = signal/(np.maximum(Ns, 1)[:, None])
    np.save(sys.argv[3], signal)
    if len(sys.argv) > 4:
        plt.figure(figsize=(15, 30))
        plt.imshow(signal, cmap='gray_r')
        plt.savefig(sys.argv[4])

def do_averageplot():
    regions = read_peakfile(open(sys.argv[2]))
    bedgraphs = read_bedgraphs(sys.stdin)
    N=2000
    signal = np.zeros(2*N)
    total = 0
    for chrom, bedgraph in bedgraphs:
        print("Reading", chrom)
        if chrom not in regions or chrom=="chrM":
            continue
        chrom_regions = regions[chrom]
        sizes = chrom_regions.ends-chrom_regions.starts
        new_regions = Regions(chrom_regions.starts-sizes/2, chrom_regions.ends+sizes/2)
        signal += signal_plot(bedgraph, new_regions, 2*N, scale_to=True)
        total+=regions[chrom].starts.size
    signal /= total
    np.save(sys.argv[3], signal)
    if len(sys.argv)>4:
        plt.plot(signal)
        plt.savefig(sys.argv[4])

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

    elif sys.argv[1] == "averageplot":
        do_averageplot()

    elif sys.argv[1] == "vplot":
        do_vplot()

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
        lines = (f"{chrom}\t{start}\t{end}\t.\t.\t{direction}" 
                 for chrom, start, end, direction in rna_bam_to_bed(sys.stdin))
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
