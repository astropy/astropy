# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Can `astropy.io.fits.open` access (remote) data using the fsspec package?
"""
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.compat.optional_deps import HAS_FSSPEC, HAS_S3FS  # noqa
from astropy.utils.data import get_pkg_data_filename

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import pytest


@pytest.mark.skipif("not HAS_FSSPEC")
def test_fsspec_local():
    """Can we use fsspec to read a local file?"""
    fn = get_pkg_data_filename('data/test0.fits')
    hdulist_classic = fits.open(fn, use_fsspec=False)
    hdulist_fsspec = fits.open(fn, use_fsspec=True)
    assert_array_equal(hdulist_classic[2].data, hdulist_fsspec[2].data)
    assert_array_equal(hdulist_classic[2].section[3:5], hdulist_fsspec[2].section[3:5])
    assert "partially read" not in repr(hdulist_classic)
    assert "partially read" in repr(hdulist_fsspec)
    hdulist_classic.close()
    hdulist_fsspec.close()


@pytest.mark.skipif("not HAS_FSSPEC")
def test_fsspec_local_write(tmpdir):
    """Can we write to a local file that was opened using fsspec?"""
    fn = get_pkg_data_filename('data/test0.fits')
    fn_tmp = tmpdir / "tmp.fits"
    with fits.open(fn, use_fsspec=True) as hdul:
        # writing to a section is never allowed
        with pytest.raises(TypeError):
            hdul[1].section[0, 0] = -999
        # however writing to .data should work
        hdul[1].data[2, 3] = -999
        assert hdul[1].data[2,3] == -999
        hdul.writeto(fn_tmp)

    # Is the new value present when we re-open the file?
    with fits.open(fn_tmp) as hdul:
        assert hdul[1].data[2,3] == -999

    # Does fsspec support `mode="update"`?
    with fits.open(str(fn_tmp), use_fsspec=True, mode="update") as hdul:
        hdul[1].data[2,3] = 42
    with fits.open(fn_tmp) as hdul:
        assert hdul[1].data[2,3] == 42


@pytest.mark.remote_data
@pytest.mark.skipif("not HAS_FSSPEC")
def test_fsspec_http():
    """Can we use fsspec to open a remote FITS file via http?"""
    uri = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HST/product/j8pu0y010_drc.fits"
    # Expected array was obtained by downloading the file locally and executing:
    # with fits.open(local_path) as hdul:
    #     hdul[1].data[1000:1002, 2000:2003]
    expected = np.array([[0.00545289, 0.0051066, -0.00034149],
                         [0.00120684, 0.00782754, 0.00546404]])
    with fits.open(uri, use_fsspec=True) as hdul:
        # Do we retrieve the expected array?
        assert_allclose(hdul[1].section[1000:1002, 2000:2003], expected, atol=1e-7)
        # The file has multiple extensions which are not yet downloaded;
        # the repr and string representation should reflect this.
        assert "partially read" in repr(hdul)
        assert "partially read" in str(hdul)


@pytest.mark.remote_data
@pytest.mark.skipif("not HAS_S3FS")
def test_fsspec_s3():
    """Can we use fsspec to open a FITS file in a public Amazon S3 bucket?"""
    uri = f"s3://stpubdata/hst/public/j8pu/j8pu0y010/j8pu0y010_drc.fits"
    # Expected array was obtained by downloading the file locally and executing:
    # with fits.open(local_path) as hdul:
    #     hdul[1].data[1000:1002, 2000:2003]
    expected = np.array([[0.00545289, 0.0051066, -0.00034149],
                         [0.00120684, 0.00782754, 0.00546404]])
    with fits.open(uri) as hdul:  # s3:// paths should default to use_fsspec=True
        # Do we retrieve the expected array?
        assert_allclose(hdul[1].section[1000:1002, 2000:2003], expected, atol=1e-7)
        # The file has multiple extensions which are not yet downloaded;
        # the repr and string representation should reflect this.
        assert "partially read" in repr(hdul)
        assert "partially read" in str(hdul)


@pytest.mark.skipif("not HAS_FSSPEC")
def test_fsspec_cutout2d():
    """Does Cutout2D work with data loaded lazily using fsspec and .section?"""
    fn = get_pkg_data_filename('data/test0.fits')
    with fits.open(fn, use_fsspec=True) as hdul:
        position = (10, 20)
        size = (2, 3)
        cutout1 = Cutout2D(hdul[1].data, position, size)
        cutout2 = Cutout2D(hdul[1].section, position, size)
        assert_allclose(cutout1.data, cutout2.data)


@pytest.mark.skipif("not HAS_FSSPEC")
def test_fsspec_compressed():
    """Does fsspec support compressed data correctly?"""
    # comp.fits[1] is a compressed image with shape (440, 300)
    fn = get_pkg_data_filename('data/comp.fits')
    with fits.open(fn, use_fsspec=True) as hdul:
        # The .data attribute should work as normal
        assert hdul[1].data[0,0] == 7
        # However the .section attribute does not support compressed data
        with pytest.raises(AttributeError) as excinfo:
            hdul[1].section[1,2]
        assert "'CompImageHDU' object has no attribute 'section'" in str(excinfo.value)
