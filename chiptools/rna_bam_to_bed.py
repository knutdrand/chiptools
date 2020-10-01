from cigar import Cigar
consuming = {"M", "D", "N", "=", "X"}

def get_splices(lines):
    parted = (line.split() for line in lines)
    spliced = (p for p in parted if "N" in p[5])
    for parts in spliced:
        chrom = parts[2]
        direction = "-" if (16 & int(parts[1])) else "+"
        cigar = list(Cigar(parts[5]).items())
        cur_pos = int(parts[3])-1
        for d, code in cigar:
            if code == "N":
                yield chrom, cur_pos, cur_pos+d, direction
            if code in consuming:
                cur_pos += d

def rna_bam_to_bed(lines):
    for line in lines:
        parts = line.split()
        direction = "-" if 16 & int(parts[1]) else "+"
        cigar = list(Cigar(parts[5]).items())
        start_pos = int(parts[3])-1
        cur_pos = start_pos
        for d, code in cigar:
            if code == "N":
                yield (parts[2], start_pos, cur_pos, direction, 1)
                start_pos = cur_pos+d
            if code in consuming:
                cur_pos += d

        yield (parts[2], start_pos, cur_pos, direction, 0)
