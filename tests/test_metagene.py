import pytest
import numpy as np
from chiptools.metagene import metagene
from chiptools.annotation import IndexedRegions
from chiptools import Regions, BedGraph
from fixtures import bedgraph

@pytest.fixture
def indexed_regions():
    regions = Regions([0,20,40], [10, 30, 50])
    return IndexedRegions(None, [0, 0,0], regions)

@pytest.fixture
def adv_indexed_regions():
    regions = Regions([0, 12, 24, 36], [10, 22, 34, 46], [1, -1, 1, -1])
    return IndexedRegions(None, [0, 1, 0, 1], regions)

def test_metagene(indexed_regions, bedgraph):
    diffs = np.zeros(30)
    metagene(bedgraph, indexed_regions, diffs)
    true = np.zeros_like(diffs)
    true[[10,15,20]] =  [2,1,1]
    assert np.all(diffs==true)

# regions = Regions([0, 12, 24, 36], [10, 22, 34, 46], [1, -1, 1, -1])
def test_adv_metagene(adv_indexed_regions, bedgraph):
    diffs = np.zeros(20)
    metagene(bedgraph, adv_indexed_regions, diffs)
    true = np.zeros_like(diffs)
    true[[10,11]] +=  [2, 1]
    true[[0, 6, 10, 17]] +=  [4, -1, -1, -1]
    assert np.all(diffs==true)


