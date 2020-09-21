from collections import namedtuple
Entry = namedtuple("Entry", ["chrom", "start", "end", "name", "qual", "strand"])

def pair_pair(first, second):
    assert first.chrom==second.chrom, (first, second)
    assert first.name[:-1]==second.name[:-1], (first, second)
    start = min(int(first.start), int(second.start))
    end = max(int(first.end), int(second.end))
    return (first.chrom, str(start), str(end))

def pair_lines(bed_lines):
    entries = (Entry(*line.split()) for line in bed_lines)
    
    while True:
        try:
            first = next(entries)
        except StopIteration:
            break
        try:
            second = next(entries)
        except StopIteration:
            yield "\t".join(pair_pair(first, first))
            break
        while first.name[:-1]!=second.name[:-1]:
            yield "\t".join(pair_pair(first, first))
            try:
                first, second = (second, next(entries))
            except StopIteration:
                yield "\t".join(pair_pair(second, second))
                break
        if first.chrom == second.chrom:
            yield "\t".join(pair_pair(first, second))
#     
# 
#     for first in entries:
#         second = next(entries)
#         if first.name[:-1]==second.name[:-1], (first, second)
#         r = next(entries)
# ...     print (x, next(it))
#     for line in lines:
        


    return ("\t".join(pair_pair(*pair)) for pair in zip(*[iter(entries)]*2))
            
