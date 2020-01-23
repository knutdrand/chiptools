from chiptools.bedgraph import BedGraph, GraphDiff
import pytest
@pytest.fixture
def bedgraph():
    return BedGraph([0, 10, 15, 25, 40], [0, 1, 2, 3, 4], size=50)

@pytest.fixture
def graphdiff():
    return GraphDiff(0, [10, 15, 25, 40], [1, 1, 1, 1], 50)

def test_scale_x(graphdiff):
    assert graphdiff.scale_x(100) == GraphDiff(0, [20, 30, 50, 80], [1, 1, 1, 1], 100)
    assert graphdiff.scale_x(75) == GraphDiff(0, [15, 22, 37, 60], [1, 1, 1, 1], 75)

def test_getitem(bedgraph):
    values = bedgraph[[10, 11, 14, 15, 16]]
    assert values == [1, 1, 1, 2, 2]

def test_getitem_slice(bedgraph):
    assert bedgraph[9:25] == BedGraph([0, 1, 6], [0, 1, 2])
    assert bedgraph[10:26] == BedGraph([0, 5, 15], [1, 2, 3])
    assert bedgraph[0:9] == BedGraph([0], [0])
    assert bedgraph[1:9] == BedGraph([0], [0])
    assert bedgraph[15:42] == BedGraph([0, 10, 25], [2, 3, 4])
    assert bedgraph[16:41] == BedGraph([0, 9, 24], [2, 3, 4])

def test_to_graph_diffs(bedgraph, graphdiff):
    assert bedgraph.to_graph_diffs() == graphdiff
    
def test_reverse(bedgraph):
    assert bedgraph.reverse() == BedGraph([0, 10, 25, 35, 40], [4, 3, 2, 1, 0], size=50)
