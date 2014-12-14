import numpy as np

from ....tests.helper import pytest
from ....io import fits

try:
    import matplotlib
    HAS_MATPLOTLIB = True
    from ..fits2bitmap import fits2bitmap
except:
    HAS_MATPLOTLIB = False


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_fits2bitmap_function(tmpdir):
    filename = tmpdir.join('test.fits').strpath
    fits.writeto(filename, np.ones((128, 128)))
    fits2bitmap(filename)
