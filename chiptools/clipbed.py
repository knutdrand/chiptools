def clip_bed(bed_lines, chrom_sizes):
    parts = (line.split(None, 3) for line in bed_lines)
    for p in parts:
        if not p[0] in chrom_sizes:
            continue
        size = chrom_sizes[p[0]]
        start = int(p[1])
        end = int(p[2])
        if  (start >= size) or (end<=0) or start>=end:
            continue
        start = min(max(0, start), size-1)
        end = min(max(1, end), size)
        p[1] = str(start)
        p[2] = str(end)
        yield "\t".join(p)
