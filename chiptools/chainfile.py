def parse_chain_line(line):
    parts = line.split()
    chrom = parts[2]
    strand = parts[4]
    start = int(parts[5])
    end = int(parts[6])
    return (chrom, start, end)

def parse_lines(lines):
    cur_chrom = None
    cur_idx = None
    cur_end = None
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        if line.startswith("chain"):
            assert cur_end is None or cur_idx == cur_end, (cur_idx, cur_end)
            cur_chrom, cur_idx, cur_end = parse_chain_line(line)
            continue
        tmp = [int(c) for c in line.split()]
        if len(tmp) != 3:
            size=tmp[0]
            dt = 0
        else:
            size, dt, _ = tmp

        yield cur_chrom, cur_idx, cur_idx+size
        cur_idx += size+dt
