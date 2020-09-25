from .annotation import Anno, Annotations

def parse_refseq_line(line):
    parts = line.split()
    _, _, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, _, exonStarts, exonEnds, _, name, _, _, _ = line.split()
    txStart, txEnd = (int(txStart), int(txEnd))
    cdsStart, cdsEnd = (int(cdsStart), int(cdsEnd))
    exonStarts = [int(s) for s in exonStarts.split(",") if s]
    exonEnds = [int(s) for s in exonEnds.split(",") if s]
    exonStarts = [int(s) for s in exonStarts]
    exonEnds = [int(s) for s in exonEnds]
    assert strand in ("+", "-"), (strand, line)
    strand = 1 if strand == "+" else -1
    return Anno(chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonStarts, exonEnds, name)

def parse_refseq_file(f):
    return Annotations(parse_refseq_line(line) for line in f)
