from ..array import SortedArray
from ..table import Table
import pytest
import numpy as np

@pytest.fixture
def array():
    # composite index
    col0 = [x % 2 for x in range(1, 11)]
    col1 = [x for x in range(1, 11)]
    t = Table([col0, col1, [x for x in range(1, 11)]])
    return SortedArray(t[t.argsort()])

def test_array_find(array):
    for i in range(1, 11):
        print("Searching for {0}".format(i))
        assert array.find((i % 2, i)) == [i]
    assert array.find((1, 4)) == []

def test_array_range(array):
    assert np.all(array.range((0, 8), (1, 3), (True, True)) == [8, 10, 1, 3])
    assert np.all(array.range((0, 8), (1, 3), (False, True)) == [10, 1, 3])
    assert np.all(array.range((0, 8), (1, 3), (True, False)) == [8, 10, 1])
