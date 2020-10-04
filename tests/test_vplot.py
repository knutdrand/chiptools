import pytest
import numpy as np
from chiptools.vplot import vplot
from chiptools import Regions, BedGraph

@pytest.fixture
def regions():
    return Regions([100, 200, 300, 400, 500], [101, 202, 303, 404, 505], [1, 1, 1, 1, 1])

@pytest.fixture
def graph():
    indices = [0, 100, 101, 200, 201, 300, 301, 400, 401, 500, 501]
    return BedGraph(indices,
                    [0, 1]*(len(indices)//2) + [0], size=600)

def test_vplot(regions, graph):
    diffs = np.zeros((6, 6))
    Ns = np.zeros(6)
    p = vplot(graph, regions, diffs, 6, Ns)
    print(p)
    # assert False
