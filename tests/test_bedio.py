from chiptools.bedIO import read_bedgraphs, read_bedfile
from chiptools import BedGraph

def test_read_bedgraphs():
    lines = ["chr1\t0\t10\t0",
             "chr1\t10\t25\t1",
             "chr1\t25\t35\t10",
             "chr2\t0\t5\t0",
             "chr2\t5\t10\t2"]
    bedgraphs = list(read_bedgraphs(lines))
    assert bedgraphs == [("chr1", BedGraph([0, 10, 25], [0, 1, 10])),
                         ("chr2", BedGraph([0, 5], [0, 2]))]
