# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Can `astropy.io.fits.open` access (remote) data using the fsspec package?
"""
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.compat.optional_deps import HAS_FSSPEC, HAS_S3FS
from astropy.utils.data import get_pkg_data_filename

if HAS_FSSPEC:
    import fsspec


@pytest.mark.skipif(not HAS_FSSPEC, reason="requires fsspec")
def test_fsspec_local():
    """Can we use fsspec to read a local file?"""
    fn = get_pkg_data_filename("data/test0.fits")
    with fits.open(fn) as hdulist_classic:
        with fits.open(fn, use_fsspec=True) as hdulist_fsspec:
            assert_array_equal(hdulist_classic[2].data, hdulist_fsspec[2].data)
            assert_array_equal(
                hdulist_classic[2].section[3:5], hdulist_fsspec[2].section[3:5]
            )


@pytest.mark.skipif(not HAS_FSSPEC, reason="requires fsspec")
def test_fsspec_local_write(tmp_path):
    """Can we write to a local file that was opened using fsspec?"""
    fn = get_pkg_data_filename("data/test0.fits")
    fn_tmp = tmp_path / "tmp.fits"
    with fits.open(fn, use_fsspec=True) as hdul:
        # writing to a section is never allowed
        with pytest.raises(TypeError):
            hdul[1].section[0, 0] = -999
        # however writing to .data should work
        hdul[1].data[2, 3] = -999
        assert hdul[1].data[2, 3] == -999
        hdul.writeto(fn_tmp)

    # Is the new value present when we re-open the file?
    with fits.open(fn_tmp) as hdul:
        assert hdul[1].data[2, 3] == -999


@pytest.mark.skipif(not HAS_FSSPEC, reason="requires fsspec")
def test_fsspec_cutout2d():
    """Does Cutout2D work with data loaded lazily using fsspec and .section?"""
    fn = get_pkg_data_filename("data/test0.fits")
    with fits.open(fn, use_fsspec=True) as hdul:
        position = (10, 20)
        size = (2, 3)
        cutout1 = Cutout2D(hdul[1].data, position, size)
        cutout2 = Cutout2D(hdul[1].section, position, size)
        assert_allclose(cutout1.data, cutout2.data)


@pytest.mark.skipif(not HAS_FSSPEC, reason="requires fsspec")
def test_fsspec_compressed():
    """Does fsspec support compressed data correctly?"""
    # comp.fits[1] is a compressed image with shape (440, 300)
    fn = get_pkg_data_filename("data/comp.fits")
    with fits.open(fn, use_fsspec=True) as hdul:
        # The .data attribute should work as normal
        assert hdul[1].data[0, 0] == 7
        # And the .section attribute should work too
        assert hdul[1].section[0, 0] == 7


@pytest.mark.remote_data
class TestFsspecRemote:
    """Test obtaining cutouts from FITS files via HTTP (from MAST) and S3 (from Amazon)."""

    def setup_class(self):
        # The test file (ibxl50020_jif.fits) is a Hubble jitter FITS file (*.jif)
        # rather than a real image, because jitter files are less likely to
        # change due to reprocessing.
        self.http_url = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HST/product/ibxl50020_jif.fits"
        self.s3_uri = "s3://stpubdata/hst/public/ibxl/ibxl50020/ibxl50020_jif.fits"
        # Random slice was selected for testing:
        self.slice = (slice(31, 33), slice(27, 30))
        # The expected cutout array below was obtained by downloading the URIs
        # listed above to a local path and and executing:
        # with fits.open(local_path) as hdul:
        #     expected_cutout = hdul[1].data[31:33, 27:30]
        self.expected_cutout = np.array([[24, 88, 228], [35, 132, 305]], dtype=np.int32)

    @pytest.mark.skipif(not HAS_FSSPEC, reason="requires fsspec")
    def test_fsspec_http(self):
        """Can we use fsspec to open a remote FITS file via http?"""
        with fits.open(self.http_url, use_fsspec=True) as hdul:
            # Do we retrieve the expected array?
            assert_array_equal(hdul[1].section[self.slice], self.expected_cutout)
            # The file has multiple extensions which are not yet downloaded;
            # the repr and string representation should reflect this.
            assert "partially read" in repr(hdul)
            assert "partially read" in str(hdul)

        # Can the user also pass an fsspec file object directly to fits open?
        with fsspec.open(self.http_url) as fileobj:
            with fits.open(fileobj) as hdul2:
                assert_array_equal(hdul2[1].section[self.slice], self.expected_cutout)
                assert "partially read" in repr(hdul)
                assert "partially read" in str(hdul)

    @pytest.mark.skipif(not HAS_S3FS, reason="requires s3fs")
    def test_fsspec_s3(self):
        """Can we use fsspec to open a FITS file in a public Amazon S3 bucket?"""
        with fits.open(
            self.s3_uri, fsspec_kwargs={"anon": True}
        ) as hdul:  # s3:// paths should default to use_fsspec=True
            # Do we retrieve the expected array?
            assert_array_equal(hdul[1].section[self.slice], self.expected_cutout)
            # The file has multiple extensions which are not yet downloaded;
            # the repr and string representation should reflect this.
            assert "partially read" in repr(hdul)
            assert "partially read" in str(hdul)

        # Can the user also pass an fsspec file object directly to fits open?
        with fsspec.open(self.s3_uri, anon=True) as fileobj:
            with fits.open(fileobj) as hdul2:
                assert_array_equal(hdul2[1].section[self.slice], self.expected_cutout)
                assert "partially read" in repr(hdul)
                assert "partially read" in str(hdul)
