import numpy as np
from chiptools.overlap import get_overlap


def _test_count_until():
    a = np.array([10, 20, 30, 70, 110, 130, 140, 150, 160, 170, 180, 190, 200, 210]).reshape((-1, 2))
    b = np.array([3, 9, 11, 18, 35, 75, 105, 133, 140, 150, 170, 180]).reshape((-1, 2))
    correct = [6, 13, 13, 48, 110, 130, 140, 150, 160, 170, 180, 190, 200, 210]
    assert get_count_until(b, np.ravel(a))


def test_overlap():
    a = np.array([10, 20, 30, 70, 110, 130, 140, 150, 160, 170, 180, 190, 200, 210]).reshape((-1, 2))
    b = np.array([3, 9, 11, 18, 35, 75, 105, 133, 140, 150, 170, 180]).reshape((-1, 2))
    overlap = get_overlap(a, b)
    assert np.all(overlap.ravel()==[7, 35, 20, 10, 0, 0, 0])
