import pytest
from astropy.io import fits
import numpy as np
from pathlib import Path

def test_fitsdiff_openfile(tmpdir):
    """Make sure that failing FITSDiff doesn't leave open files"""
    path1 = str(tmpdir.join("file1.fits"))
    path2 = str(tmpdir.join("file2.fits"))

    hdulist = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=np.zeros((10)))])
    hdulist.writeto(path1)
    hdulist[1].data[0] = 1
    hdulist.writeto(path2)

    diff = fits.FITSDiff(path1, path2)
    assert diff.identical, diff.report()
