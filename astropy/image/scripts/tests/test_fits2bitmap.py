import numpy as np
from astropy.io import fits

from ..fits2bitmap import fits2bitmap

def test_fits2bitmap_function(tmpdir):
    filename = tmpdir.join('test.fits').strpath
    fits.writeto(filename, np.ones((128, 128)))
    fits2bitmap(filename)
