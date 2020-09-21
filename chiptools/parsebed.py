from collections import defaultdict
dtype= [('start', np.unicode_, 16), ('grades', np.float64, (2,))]

def get_chroms(lines,  buffer_size=10000):
    chroms = defaultdict(list)
    for line in lines:
        chrom, start, end = line.split("\t", 3)[:3]
        chroms[chrom].append((int(start), int(end), attrs))
    for chrom, coords in chroms.items():
        chroms[chrom] = np.array(coords, dtype="int")
    return {chrom: coords[coords[:, 0].argsort()] for chrom, coords in chroms.items()}

def print_chroms(chroms):
    for chrom, lines in chroms.items():
        for line in lines:
            print(chrom+"\t"+"\t".join(line))


def parse_fasta(filename):
    


def parse_fastq(filename):
    lines = open(filename)
    args = [iter(lines)] * 4
    entries = zip_longest(*args, fillvalue=fillvalue)
    return (FastqEntry(name, seq, qual) for name, seq, _ qual in entries)

def parse_bed(filename, n=6):
    assert n <= 12
    header = ("chrom", "start", "end", "name", "score", "strand", 
              "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
    BedEntry = namedtuple(f"Bed{n}Entry", header[:n])
    types = (str, int, int, str, float, str, int, int, str, int, int, int)
    split_lines = (lines.split(None, n)[:n] for line in open(filename))
    return (BedEntry(*(f(part) for f, part in zip(types, line)))
            for line in split_lines)

    
