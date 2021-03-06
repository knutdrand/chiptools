from fixtures import bedgraph, regions_10b
from chiptools.signalplot import signal_plot
import numpy as np

def test_signalplot(bedgraph, regions_10b):
    diffs = np.zeros(10)
    signal_plot(bedgraph, regions_10b, diffs)
    signal = np.cumsum(diffs)
    true = np.sum([[0,0,0,0,0,0,0,0,1,1],
                   [2,2,2,2,2,2,2,2,1,1],
                   [2,2,2,2,2,2,2,2,3,3]], axis=0)
    assert np.all(signal==true)
