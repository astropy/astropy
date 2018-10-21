# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import tempfile
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal

from ...tests.helper import assert_quantity_allclose
from ..utils import (extract_array, add_array, subpixel_indices,
                     block_reduce, block_replicate,
                     overlap_slices, NoOverlapError, PartialOverlapError,
                     Cutout2D, make_cutouts, cutouts_from_fits)
from ...wcs import WCS, Sip
from ...wcs.utils import proj_plane_pixel_area
from ...coordinates import SkyCoord
from ... import units as u
from ...io import fits
from ...io.fits.tests import FitsTestCase
from ...table import Table
from ..nddata import NDData

try:
    import skimage  # pylint: disable=W0611
    HAS_SKIMAGE = True
except ImportError:
    HAS_SKIMAGE = False


test_positions = [(10.52, 3.12), (5.62, 12.97), (31.33, 31.77),
                  (0.46, 0.94), (20.45, 12.12), (42.24, 24.42)]

test_position_indices = [(0, 3), (0, 2), (4, 1),
                         (4, 2), (4, 3), (3, 4)]

test_slices = [slice(10.52, 3.12), slice(5.62, 12.97),
               slice(31.33, 31.77), slice(0.46, 0.94),
               slice(20.45, 12.12), slice(42.24, 24.42)]

subsampling = 5

test_pos_bad = [(-1, -4), (-1, 0), (6, 2), (6, 6)]


def test_slices_different_dim():
    '''Overlap from arrays with different number of dim is undefined.'''
    with pytest.raises(ValueError) as e:
        overlap_slices((4, 5, 6), (1, 2), (0, 0))
    assert "the same number of dimensions" in str(e.value)


def test_slices_pos_different_dim():
    '''Position must have same dim as arrays.'''
    with pytest.raises(ValueError) as e:
        overlap_slices((4, 5), (1, 2), (0, 0, 3))
    assert "the same number of dimensions" in str(e.value)


@pytest.mark.parametrize('pos', test_pos_bad)
def test_slices_no_overlap(pos):
    '''If there is no overlap between arrays, an error should be raised.'''
    with pytest.raises(NoOverlapError):
        overlap_slices((5, 5), (2, 2), pos)


def test_slices_partial_overlap():
    '''Compute a slice for partially overlapping arrays.'''
    temp = overlap_slices((5,), (3,), (0,))
    assert temp == ((slice(0, 2, None),), (slice(1, 3, None),))

    temp = overlap_slices((5,), (3,), (0,), mode='partial')
    assert temp == ((slice(0, 2, None),), (slice(1, 3, None),))

    for pos in [0, 4]:
        with pytest.raises(PartialOverlapError) as e:
            temp = overlap_slices((5,), (3,), (pos,), mode='strict')
        assert 'Arrays overlap only partially.' in str(e.value)


def test_slices_overlap_wrong_mode():
    '''Call overlap_slices with non-existing mode.'''
    with pytest.raises(ValueError) as e:
        overlap_slices((5,), (3,), (0,), mode='full')
    assert "Mode can be only" in str(e.value)


def test_extract_array_wrong_mode():
    '''Call extract_array with non-existing mode.'''
    with pytest.raises(ValueError) as e:
        extract_array(np.arange(4), (2, ), (0, ), mode='full')
    assert "Valid modes are 'partial', 'trim', and 'strict'." == str(e.value)


def test_extract_array_1d_even():
    '''Extract 1 d arrays.

    All dimensions are treated the same, so we can test in 1 dim.
    '''
    assert np.all(extract_array(np.arange(4), (2, ), (0, ), fill_value=-99) == np.array([-99, 0]))
    for i in [1, 2, 3]:
        assert np.all(extract_array(np.arange(4), (2, ), (i, )) == np.array([i - 1, i]))
    assert np.all(extract_array(np.arange(4.), (2, ), (4, ), fill_value=np.inf) == np.array([3, np.inf]))


def test_extract_array_1d_odd():
    '''Extract 1 d arrays.

    All dimensions are treated the same, so we can test in 1 dim.
    The first few lines test the most error-prone part: Extraction of an
    array on the boundaries.
    Additional tests (e.g. dtype of return array) are done for the last
    case only.
    '''
    assert np.all(extract_array(np.arange(4), (3,), (-1, ), fill_value=-99) == np.array([-99, -99, 0]))
    assert np.all(extract_array(np.arange(4), (3,), (0, ), fill_value=-99) == np.array([-99, 0, 1]))
    for i in [1, 2]:
        assert np.all(extract_array(np.arange(4), (3,), (i, )) == np.array([i-1, i, i+1]))
    assert np.all(extract_array(np.arange(4), (3,), (3, ), fill_value=-99) == np.array([2, 3, -99]))
    arrayin = np.arange(4.)
    extracted = extract_array(arrayin, (3,), (4, ))
    assert extracted[0] == 3
    assert np.isnan(extracted[1])  # since I cannot use `==` to test for nan
    assert extracted.dtype == arrayin.dtype


def test_extract_array_1d():
    """In 1d, shape can be int instead of tuple"""
    assert np.all(extract_array(np.arange(4), 3, (-1, ), fill_value=-99) == np.array([-99, -99, 0]))
    assert np.all(extract_array(np.arange(4), 3, -1, fill_value=-99) == np.array([-99, -99, 0]))


def test_extract_Array_float():
    """integer is at bin center"""
    for a in np.arange(2.51, 3.49, 0.1):
        assert np.all(extract_array(np.arange(5), 3, a) == np.array([2, 3, 4]))


def test_extract_array_1d_trim():
    '''Extract 1 d arrays.

    All dimensions are treated the same, so we can test in 1 dim.
    '''
    assert np.all(extract_array(np.arange(4), (2, ), (0, ), mode='trim') == np.array([0]))
    for i in [1, 2, 3]:
        assert np.all(extract_array(np.arange(4), (2, ), (i, ), mode='trim') == np.array([i - 1, i]))
    assert np.all(extract_array(np.arange(4.), (2, ), (4, ), mode='trim') == np.array([3]))


@pytest.mark.parametrize('mode', ['partial', 'trim', 'strict'])
def test_extract_array_easy(mode):
    """
    Test extract_array utility function.

    Test by extracting an array of ones out of an array of zeros.
    """
    large_test_array = np.zeros((11, 11))
    small_test_array = np.ones((5, 5))
    large_test_array[3:8, 3:8] = small_test_array
    extracted_array = extract_array(large_test_array, (5, 5), (5, 5), mode=mode)
    assert np.all(extracted_array == small_test_array)


def test_extract_array_return_pos():
    '''Check that the return position is calculated correctly.

    The result will differ by mode. All test here are done in 1d because it's
    easier to construct correct test cases.
    '''
    large_test_array = np.arange(5)
    for i in np.arange(-1, 6):
        extracted, new_pos = extract_array(large_test_array, 3, i,
                                           mode='partial', return_position=True)
        assert new_pos == (1, )
    # Now check an array with an even number
    for i, expected in zip([1.49, 1.51, 3], [1.49, 0.51, 1]):
        extracted, new_pos = extract_array(large_test_array, (2,), (i,),
                                           mode='strict', return_position=True)
        assert new_pos == (expected, )
    # For mode='trim' the answer actually depends
    for i, expected in zip(np.arange(-1, 6), (-1, 0, 1, 1, 1, 1, 1)):
        extracted, new_pos = extract_array(large_test_array, (3,), (i,),
                                           mode='trim', return_position=True)
        assert new_pos == (expected, )


def test_add_array_odd_shape():
    """
    Test add_array utility function.

    Test by adding an array of ones out of an array of zeros.
    """
    large_test_array = np.zeros((11, 11))
    small_test_array = np.ones((5, 5))
    large_test_array_ref = large_test_array.copy()
    large_test_array_ref[3:8, 3:8] += small_test_array

    added_array = add_array(large_test_array, small_test_array, (5, 5))
    assert np.all(added_array == large_test_array_ref)


def test_add_array_even_shape():
    """
    Test add_array_2D utility function.

    Test by adding an array of ones out of an array of zeros.
    """
    large_test_array = np.zeros((11, 11))
    small_test_array = np.ones((4, 4))
    large_test_array_ref = large_test_array.copy()
    large_test_array_ref[0:2, 0:2] += small_test_array[2:4, 2:4]

    added_array = add_array(large_test_array, small_test_array, (0, 0))
    assert np.all(added_array == large_test_array_ref)


@pytest.mark.parametrize(('position', 'subpixel_index'),
                         zip(test_positions, test_position_indices))
def test_subpixel_indices(position, subpixel_index):
    """
    Test subpixel_indices utility function.

    Test by asserting that the function returns correct results for
    given test values.
    """
    assert np.all(subpixel_indices(position, subsampling) == subpixel_index)


@pytest.mark.skipif('not HAS_SKIMAGE')
class TestBlockReduce:
    def test_1d(self):
        """Test 1D array."""
        data = np.arange(4)
        expected = np.array([1, 5])
        result = block_reduce(data, 2)
        assert np.all(result == expected)

    def test_1d_mean(self):
        """Test 1D array with func=np.mean."""
        data = np.arange(4)
        block_size = 2.
        expected = block_reduce(data, block_size, func=np.sum) / block_size
        result_mean = block_reduce(data, block_size, func=np.mean)
        assert np.all(result_mean == expected)

    def test_2d(self):
        """Test 2D array."""
        data = np.arange(4).reshape(2, 2)
        expected = np.array([[6]])
        result = block_reduce(data, 2)
        assert np.all(result == expected)

    def test_2d_mean(self):
        """Test 2D array with func=np.mean."""
        data = np.arange(4).reshape(2, 2)
        block_size = 2.
        expected = (block_reduce(data, block_size, func=np.sum) /
                    block_size**2)
        result = block_reduce(data, block_size, func=np.mean)
        assert np.all(result == expected)

    def test_2d_trim(self):
        """
        Test trimming of 2D array when size is not perfectly divisible
        by block_size.
        """

        data1 = np.arange(15).reshape(5, 3)
        result1 = block_reduce(data1, 2)
        data2 = data1[0:4, 0:2]
        result2 = block_reduce(data2, 2)
        assert np.all(result1 == result2)

    def test_block_size_broadcasting(self):
        """Test scalar block_size broadcasting."""
        data = np.arange(16).reshape(4, 4)
        result1 = block_reduce(data, 2)
        result2 = block_reduce(data, (2, 2))
        assert np.all(result1 == result2)

    def test_block_size_len(self):
        """Test block_size length."""
        data = np.ones((2, 2))
        with pytest.raises(ValueError):
            block_reduce(data, (2, 2, 2))


@pytest.mark.skipif('not HAS_SKIMAGE')
class TestBlockReplicate:
    def test_1d(self):
        """Test 1D array."""
        data = np.arange(2)
        expected = np.array([0, 0, 0.5, 0.5])
        result = block_replicate(data, 2)
        assert np.all(result == expected)

    def test_1d_conserve_sum(self):
        """Test 1D array with conserve_sum=False."""
        data = np.arange(2)
        block_size = 2.
        expected = block_replicate(data, block_size) * block_size
        result = block_replicate(data, block_size, conserve_sum=False)
        assert np.all(result == expected)

    def test_2d(self):
        """Test 2D array."""
        data = np.arange(2).reshape(2, 1)
        expected = np.array([[0, 0], [0, 0], [0.25, 0.25], [0.25, 0.25]])
        result = block_replicate(data, 2)
        assert np.all(result == expected)

    def test_2d_conserve_sum(self):
        """Test 2D array with conserve_sum=False."""
        data = np.arange(6).reshape(2, 3)
        block_size = 2.
        expected = block_replicate(data, block_size) * block_size**2
        result = block_replicate(data, block_size, conserve_sum=False)
        assert np.all(result == expected)

    def test_block_size_broadcasting(self):
        """Test scalar block_size broadcasting."""
        data = np.arange(4).reshape(2, 2)
        result1 = block_replicate(data, 2)
        result2 = block_replicate(data, (2, 2))
        assert np.all(result1 == result2)

    def test_block_size_len(self):
        """Test block_size length."""
        data = np.arange(5)
        with pytest.raises(ValueError):
            block_replicate(data, (2, 2))


class TestCutout2D:
    def setup_class(self):
        self.data = np.arange(20.).reshape(5, 4)
        self.position = SkyCoord('13h11m29.96s -01d19m18.7s', frame='icrs')
        wcs = WCS(naxis=2)
        rho = np.pi / 3.
        scale = 0.05 / 3600.
        wcs.wcs.cd = [[scale*np.cos(rho), -scale*np.sin(rho)],
                      [scale*np.sin(rho), scale*np.cos(rho)]]
        wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        wcs.wcs.crval = [self.position.ra.to_value(u.deg),
                         self.position.dec.to_value(u.deg)]
        wcs.wcs.crpix = [3, 3]
        self.wcs = wcs

        # add SIP
        sipwcs = wcs.deepcopy()
        sipwcs.wcs.ctype = ['RA---TAN-SIP', 'DEC--TAN-SIP']
        a = np.array(
            [[0, 0, 5.33092692e-08, 3.73753773e-11, -2.02111473e-13],
             [0, 2.44084308e-05, 2.81394789e-11, 5.17856895e-13, 0.0],
             [-2.41334657e-07, 1.29289255e-10, 2.35753629e-14, 0.0, 0.0],
             [-2.37162007e-10, 5.43714947e-13, 0.0, 0.0, 0.0],
             [ -2.81029767e-13, 0.0, 0.0, 0.0, 0.0]]
        )
        b = np.array(
            [[0, 0, 2.99270374e-05, -2.38136074e-10, 7.23205168e-13],
             [0, -1.71073858e-07, 6.31243431e-11, -5.16744347e-14, 0.0],
             [6.95458963e-06, -3.08278961e-10, -1.75800917e-13, 0.0, 0.0],
             [3.51974159e-11, 5.60993016e-14, 0.0, 0.0, 0.0],
             [-5.92438525e-13, 0.0, 0.0, 0.0, 0.0]]
        )
        sipwcs.sip = Sip(a, b, None, None, wcs.wcs.crpix)
        sipwcs.wcs.set()
        self.sipwcs = sipwcs

    def test_cutout(self):
        sizes = [3, 3*u.pixel, (3, 3), (3*u.pixel, 3*u.pix), (3., 3*u.pixel),
                 (2.9, 3.3)]
        for size in sizes:
            position = (2.1, 1.9)
            c = Cutout2D(self.data, position, size)
            assert c.data.shape == (3, 3)
            assert c.data[1, 1] == 10
            assert c.origin_original == (1, 1)
            assert c.origin_cutout == (0, 0)
            assert c.input_position_original == position
            assert_allclose(c.input_position_cutout, (1.1, 0.9))
            assert c.position_original == (2., 2.)
            assert c.position_cutout == (1., 1.)
            assert c.center_original == (2., 2.)
            assert c.center_cutout == (1., 1.)
            assert c.bbox_original == ((1, 3), (1, 3))
            assert c.bbox_cutout == ((0, 2), (0, 2))
            assert c.slices_original == (slice(1, 4), slice(1, 4))
            assert c.slices_cutout == (slice(0, 3), slice(0, 3))

    def test_size_length(self):
        with pytest.raises(ValueError):
            Cutout2D(self.data, (2, 2), (1, 1, 1))

    def test_size_units(self):
        for size in [3 * u.cm, (3, 3 * u.K)]:
            with pytest.raises(ValueError):
                Cutout2D(self.data, (2, 2), size)

    def test_size_pixel(self):
        """
        Check size in derived pixel units.
        """
        size = 0.3*u.arcsec / (0.1*u.arcsec/u.pixel)
        c = Cutout2D(self.data, (2, 2), size)
        assert c.data.shape == (3, 3)
        assert c.data[0, 0] == 5
        assert c.slices_original == (slice(1, 4), slice(1, 4))
        assert c.slices_cutout == (slice(0, 3), slice(0, 3))

    def test_size_angle(self):
        c = Cutout2D(self.data, (2, 2), (0.1*u.arcsec), wcs=self.wcs)
        assert c.data.shape == (2, 2)
        assert c.data[0, 0] == 5
        assert c.slices_original == (slice(1, 3), slice(1, 3))
        assert c.slices_cutout == (slice(0, 2), slice(0, 2))

    def test_size_angle_without_wcs(self):
        with pytest.raises(ValueError):
            Cutout2D(self.data, (2, 2), (3, 3 * u.arcsec))

    def test_cutout_trim_overlap(self):
        c = Cutout2D(self.data, (0, 0), (3, 3), mode='trim')
        assert c.data.shape == (2, 2)
        assert c.data[0, 0] == 0
        assert c.slices_original == (slice(0, 2), slice(0, 2))
        assert c.slices_cutout == (slice(0, 2), slice(0, 2))

    def test_cutout_partial_overlap(self):
        c = Cutout2D(self.data, (0, 0), (3, 3), mode='partial')
        assert c.data.shape == (3, 3)
        assert c.data[1, 1] == 0
        assert c.slices_original == (slice(0, 2), slice(0, 2))
        assert c.slices_cutout == (slice(1, 3), slice(1, 3))

    def test_cutout_partial_overlap_fill_value(self):
        fill_value = -99
        c = Cutout2D(self.data, (0, 0), (3, 3), mode='partial',
                     fill_value=fill_value)
        assert c.data.shape == (3, 3)
        assert c.data[1, 1] == 0
        assert c.data[0, 0] == fill_value

    def test_copy(self):
        data = np.copy(self.data)
        c = Cutout2D(data, (2, 3), (3, 3))
        xy = (0, 0)
        value = 100.
        c.data[xy] = value
        xy_orig = c.to_original_position(xy)
        yx = xy_orig[::-1]
        assert data[yx] == value

        data = np.copy(self.data)
        c2 = Cutout2D(self.data, (2, 3), (3, 3), copy=True)
        c2.data[xy] = value
        assert data[yx] != value

    def test_to_from_large(self):
        position = (2, 2)
        c = Cutout2D(self.data, position, (3, 3))
        xy = (0, 0)
        result = c.to_cutout_position(c.to_original_position(xy))
        assert_allclose(result, xy)

    def test_skycoord_without_wcs(self):
        with pytest.raises(ValueError):
            Cutout2D(self.data, self.position, (3, 3))

    def test_skycoord(self):
        c = Cutout2D(self.data, self.position, (3, 3), wcs=self.wcs)
        skycoord_original = self.position.from_pixel(c.center_original[1],
                                                     c.center_original[0],
                                                     self.wcs)
        skycoord_cutout = self.position.from_pixel(c.center_cutout[1],
                                                   c.center_cutout[0], c.wcs)
        assert_quantity_allclose(skycoord_original.ra, skycoord_cutout.ra)
        assert_quantity_allclose(skycoord_original.dec, skycoord_cutout.dec)

    def test_skycoord_partial(self):
        c = Cutout2D(self.data, self.position, (3, 3), wcs=self.wcs,
                     mode='partial')
        skycoord_original = self.position.from_pixel(c.center_original[1],
                                                     c.center_original[0],
                                                     self.wcs)
        skycoord_cutout = self.position.from_pixel(c.center_cutout[1],
                                                   c.center_cutout[0], c.wcs)
        assert_quantity_allclose(skycoord_original.ra, skycoord_cutout.ra)
        assert_quantity_allclose(skycoord_original.dec, skycoord_cutout.dec)

    def test_naxis_update(self):
        xsize = 2
        ysize = 3
        c = Cutout2D(self.data, self.position, (ysize, xsize), wcs=self.wcs)
        assert c.wcs._naxis[0] == xsize
        assert c.wcs._naxis[1] == ysize

    def test_crpix_maps_to_crval(self):
        w = Cutout2D(self.data, (0, 0), (3, 3), wcs=self.sipwcs,
                     mode='partial').wcs
        pscale = np.sqrt(proj_plane_pixel_area(w))
        assert_allclose(
            w.wcs_pix2world(*w.wcs.crpix, 1), w.wcs.crval,
            rtol=0.0, atol=1e-6 * pscale
        )
        assert_allclose(
            w.all_pix2world(*w.wcs.crpix, 1), w.wcs.crval,
            rtol=0.0, atol=1e-6 * pscale
        )


class TestCutoutTools(FitsTestCase):
    def setup(self):
        self.temp_dir = tempfile.mkdtemp(prefix='fits-test-')
        fits.conf.enable_record_valued_keyword_cards = True
        fits.conf.extension_name_case_sensitive = False
        fits.conf.strip_header_whitespace = True
        fits.conf.use_memmap = True

    @staticmethod
    def construct_test_image():
        # Construct data: 9 X 9 array where
        # the value at each pixel (x, y) is x*10 + y
        data = np.array([i * 9 + np.arange(9) for i in range(9)])

        # Construct a WCS for test image
        # N.B: Pixel scale should be 1 deg/pix
        w = WCS()
        w.wcs.crpix = [5, 5]
        w.wcs.crval = [0, 45]
        w.wcs.cdelt = [-1, 1]
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        # Construct fits header
        h = w.to_header()

        # Construct Image and return:
        return fits.PrimaryHDU(data, h)

    def test_make_cutouts_inputs(self):
        # Construct image:
        image_hdu = self.construct_test_image()

        # Construct catalog
        ra = [0, 1] * u.deg
        dec = [45, 46] * u.deg
        ids = ["Target_1", "Target_2"]
        cutout_width = cutout_height = [4.0, 4.0] * u.pix

        catalog = Table(
            data=[ids, ra, dec, cutout_width, cutout_height],
            names=['id', 'ra', 'dec', 'cutout_width', 'cutout_height'])

        # From NumPy Array and WCS:
        array = image_hdu.data
        w = WCS(image_hdu.header)
        assert None not in make_cutouts(array, catalog, wcs=w)

        # From NDData Array Only:
        w = WCS(image_hdu.header)
        array = NDData(data=image_hdu.data, wcs=w)
        assert None not in make_cutouts(array, catalog)

    def test_cutouts_from_fits_inputs(self):
        # Construct image:
        image_hdu = self.construct_test_image()

        # Construct catalog
        ra = [0, 1] * u.deg
        dec = [45, 46] * u.deg
        ids = ["Target_1", "Target_2"]
        cutout_width = cutout_height = [4.0, 4.0] * u.pix

        catalog = Table(
            data=[ids, ra, dec, cutout_width, cutout_height],
            names=['id', 'ra', 'dec', 'cutout_width', 'cutout_height'])

        # Test with no rotation

        # - From memory:
        assert None not in cutouts_from_fits(image_hdu, catalog)

        # - From fits file:
        fits_file = self.temp("input_image.fits")
        image_hdu.writeto(fits_file)
        assert None not in cutouts_from_fits(fits_file, catalog)

        # - From HDUList:
        hdu_list = fits.HDUList(image_hdu)
        assert None not in cutouts_from_fits(hdu_list, catalog)

        # - From ECSV file
        ecsv_file = self.temp("input_catalog.ecsv")
        catalog.write(ecsv_file, format="ascii.ecsv")
        assert None not in cutouts_from_fits(image_hdu, ecsv_file)

        # Test with rotation column:
        pa = [90, 45] * u.deg
        catalog.add_column(pa, name="cutout_pa")
        assert None not in cutouts_from_fits(image_hdu, catalog)

    def test_cutout_tool_correctness(self):
        # Construct image:
        image_hdu = self.construct_test_image()

        # Construct catalog
        ra = [0] * u.deg  # Center pixel
        dec = [45] * u.deg  # Center pixel
        ids = ["target_1"]

        # Cutout should be 3 by 3:
        cutout_width = cutout_height = [3.0] * u.pix

        catalog = Table(
            data=[ids, ra, dec, cutout_width, cutout_height],
            names=['id', 'ra', 'dec', 'cutout_width', 'cutout_height'])

        # Call cutout tool and extract the first cutout:
        cutout = cutouts_from_fits(image_hdu, catalog)[0]

        # check if values are correct:
        w_orig = WCS(image_hdu.header)
        w_new = WCS(cutout.header)

        for x_new, x_orig in enumerate(range(3, 6)):
            for y_new, y_orig in enumerate(range(3, 6)):
                coords_orig = SkyCoord.from_pixel(x_orig, y_orig, w_orig, origin=0)
                coords_new = SkyCoord.from_pixel(x_new, y_new, w_new, origin=0)

                assert_almost_equal(image_hdu.data[x_orig][y_orig], cutout.data[x_new][y_new])
                assert_almost_equal(coords_orig.ra.value, coords_new.ra.value)
                assert_almost_equal(coords_orig.dec.value, coords_new.dec.value)

        # Test for rotation:
        pa = [180] * u.deg
        catalog.add_column(pa, name="cutout_pa")

        # Call cutout tool:
        cutout = cutouts_from_fits(image_hdu, catalog)[0]

        # check if values are correct:
        w_orig = WCS(image_hdu.header)
        w_new = WCS(cutout.header)

        for x_new, x_orig in enumerate(range(5, 2, -1)):
            for y_new, y_orig in enumerate(range(5, 2, -1)):
                coords_orig = SkyCoord.from_pixel(x_orig, y_orig, w_orig, origin=0)
                coords_new = SkyCoord.from_pixel(x_new, y_new, w_new, origin=0)

                assert_almost_equal(image_hdu.data[x_orig][y_orig], cutout.data[x_new][y_new])
                assert_almost_equal(coords_orig.ra.value, coords_new.ra.value)
                assert_almost_equal(coords_orig.dec.value, coords_new.dec.value)
