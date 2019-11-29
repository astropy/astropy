# Licensed under a 3-clause BSD style license - see LICENSE.rst

import textwrap

import numpy as np
import pytest

from astropy.io import fits
from astropy.nddata.nduncertainty import (
    StdDevUncertainty, MissingDataAssociationException, VarianceUncertainty,
    InverseVariance)
from astropy import units as u
from astropy import log
from astropy.wcs import WCS, FITSFixedWarning
from astropy.tests.helper import catch_warnings
from astropy.utils import NumpyRNGContext
from astropy.utils.data import (get_pkg_data_filename, get_pkg_data_filenames,
                                get_pkg_data_contents)
from astropy.utils.exceptions import AstropyWarning

from astropy.nddata.ccddata import CCDData
from astropy.nddata import _testing as nd_testing
from astropy.table import Table

DEFAULT_DATA_SIZE = 100

with NumpyRNGContext(123):
    _random_array = np.random.normal(size=[DEFAULT_DATA_SIZE, DEFAULT_DATA_SIZE])


def create_ccd_data():
    """
    Return a CCDData object of size DEFAULT_DATA_SIZE x DEFAULT_DATA_SIZE
    with units of ADU.
    """
    data = _random_array.copy()
    fake_meta = {'my_key': 42, 'your_key': 'not 42'}
    ccd = CCDData(data, unit=u.adu)
    ccd.header = fake_meta
    return ccd


def test_ccddata_empty():
    with pytest.raises(TypeError):
        CCDData()  # empty initializer should fail


def test_ccddata_must_have_unit():
    with pytest.raises(ValueError):
        CCDData(np.zeros([2, 2]))


def test_ccddata_unit_cannot_be_set_to_none():
    ccd_data = create_ccd_data()
    with pytest.raises(TypeError):
        ccd_data.unit = None


def test_ccddata_meta_header_conflict():
    with pytest.raises(ValueError) as exc:
        CCDData([1, 2, 3], unit='', meta={1: 1}, header={2: 2})
        assert "can't have both header and meta." in str(exc.value)


def test_ccddata_simple():
    ccd_data = create_ccd_data()
    assert ccd_data.shape == (DEFAULT_DATA_SIZE, DEFAULT_DATA_SIZE)
    assert ccd_data.size == DEFAULT_DATA_SIZE * DEFAULT_DATA_SIZE
    assert ccd_data.dtype == np.dtype(float)


def test_ccddata_init_with_string_electron_unit():
    ccd = CCDData(np.zeros([2, 2]), unit="electron")
    assert ccd.unit is u.electron


def test_initialize_from_FITS(tmpdir):
    ccd_data = create_ccd_data()
    hdu = fits.PrimaryHDU(ccd_data)
    hdulist = fits.HDUList([hdu])
    filename = tmpdir.join('afile.fits').strpath
    hdulist.writeto(filename)
    cd = CCDData.read(filename, unit=u.electron)
    assert cd.shape == (DEFAULT_DATA_SIZE, DEFAULT_DATA_SIZE)
    assert cd.size == DEFAULT_DATA_SIZE * DEFAULT_DATA_SIZE
    assert np.issubdtype(cd.data.dtype, np.floating)
    for k, v in hdu.header.items():
        assert cd.meta[k] == v


def test_initialize_from_fits_with_unit_in_header(tmpdir):
    fake_img = np.zeros([2, 2])
    hdu = fits.PrimaryHDU(fake_img)
    hdu.header['bunit'] = u.adu.to_string()
    filename = tmpdir.join('afile.fits').strpath
    hdu.writeto(filename)
    ccd = CCDData.read(filename)
    # ccd should pick up the unit adu from the fits header...did it?
    assert ccd.unit is u.adu

    # An explicit unit in the read overrides any unit in the FITS file
    ccd2 = CCDData.read(filename, unit="photon")
    assert ccd2.unit is u.photon


def test_initialize_from_fits_with_ADU_in_header(tmpdir):
    fake_img = np.zeros([2, 2])
    hdu = fits.PrimaryHDU(fake_img)
    hdu.header['bunit'] = 'ADU'
    filename = tmpdir.join('afile.fits').strpath
    hdu.writeto(filename)
    ccd = CCDData.read(filename)
    # ccd should pick up the unit adu from the fits header...did it?
    assert ccd.unit is u.adu


def test_initialize_from_fits_with_invalid_unit_in_header(tmpdir):
    hdu = fits.PrimaryHDU(np.ones((2, 2)))
    hdu.header['bunit'] = 'definetely-not-a-unit'
    filename = tmpdir.join('afile.fits').strpath
    hdu.writeto(filename)
    with pytest.raises(ValueError):
        CCDData.read(filename)


def test_initialize_from_fits_with_technically_invalid_but_not_really(tmpdir):
    hdu = fits.PrimaryHDU(np.ones((2, 2)))
    hdu.header['bunit'] = 'ELECTRONS/S'
    filename = tmpdir.join('afile.fits').strpath
    hdu.writeto(filename)
    ccd = CCDData.read(filename)
    assert ccd.unit == u.electron/u.s


def test_initialize_from_fits_with_data_in_different_extension(tmpdir):
    fake_img = np.arange(4).reshape(2, 2)
    hdu1 = fits.PrimaryHDU()
    hdu2 = fits.ImageHDU(fake_img)
    hdus = fits.HDUList([hdu1, hdu2])
    filename = tmpdir.join('afile.fits').strpath
    hdus.writeto(filename)
    with catch_warnings(FITSFixedWarning) as w:
        ccd = CCDData.read(filename, unit='adu')
    assert len(w) == 0
    # ccd should pick up the unit adu from the fits header...did it?
    np.testing.assert_array_equal(ccd.data, fake_img)
    # check that the header is the combined header
    assert hdu2.header + hdu1.header == ccd.header


def test_initialize_from_fits_with_extension(tmpdir):
    fake_img1 = np.zeros([2, 2])
    fake_img2 = np.arange(4).reshape(2, 2)
    hdu0 = fits.PrimaryHDU()
    hdu1 = fits.ImageHDU(fake_img1)
    hdu2 = fits.ImageHDU(fake_img2)
    hdus = fits.HDUList([hdu0, hdu1, hdu2])
    filename = tmpdir.join('afile.fits').strpath
    hdus.writeto(filename)
    ccd = CCDData.read(filename, hdu=2, unit='adu')
    # ccd should pick up the unit adu from the fits header...did it?
    np.testing.assert_array_equal(ccd.data, fake_img2)


def test_write_unit_to_hdu():
    ccd_data = create_ccd_data()
    ccd_unit = ccd_data.unit
    hdulist = ccd_data.to_hdu()
    assert 'bunit' in hdulist[0].header
    assert hdulist[0].header['bunit'] == ccd_unit.to_string()


def test_initialize_from_FITS_bad_keyword_raises_error(tmpdir):
    # There are two fits.open keywords that are not permitted in ccdproc:
    #     do_not_scale_image_data and scale_back
    ccd_data = create_ccd_data()
    filename = tmpdir.join('test.fits').strpath
    ccd_data.write(filename)

    with pytest.raises(TypeError):
        CCDData.read(filename, unit=ccd_data.unit,
                     do_not_scale_image_data=True)
    with pytest.raises(TypeError):
        CCDData.read(filename, unit=ccd_data.unit, scale_back=True)


def test_ccddata_writer(tmpdir):
    ccd_data = create_ccd_data()
    filename = tmpdir.join('test.fits').strpath
    ccd_data.write(filename)

    ccd_disk = CCDData.read(filename, unit=ccd_data.unit)
    np.testing.assert_array_equal(ccd_data.data, ccd_disk.data)


def test_ccddata_meta_is_case_sensitive():
    ccd_data = create_ccd_data()
    key = 'SoMeKEY'
    ccd_data.meta[key] = 10
    assert key.lower() not in ccd_data.meta
    assert key.upper() not in ccd_data.meta
    assert key in ccd_data.meta


def test_ccddata_meta_is_not_fits_header():
    ccd_data = create_ccd_data()
    ccd_data.meta = {'OBSERVER': 'Edwin Hubble'}
    assert not isinstance(ccd_data.meta, fits.Header)


def test_fromMEF(tmpdir):
    ccd_data = create_ccd_data()
    hdu = fits.PrimaryHDU(ccd_data)
    hdu2 = fits.PrimaryHDU(2 * ccd_data.data)
    hdulist = fits.HDUList(hdu)
    hdulist.append(hdu2)
    filename = tmpdir.join('afile.fits').strpath
    hdulist.writeto(filename)
    # by default, we reading from the first extension
    cd = CCDData.read(filename, unit=u.electron)
    np.testing.assert_array_equal(cd.data, ccd_data.data)
    # but reading from the second should work too
    cd = CCDData.read(filename, hdu=1, unit=u.electron)
    np.testing.assert_array_equal(cd.data, 2 * ccd_data.data)


def test_metafromheader():
    hdr = fits.header.Header()
    hdr['observer'] = 'Edwin Hubble'
    hdr['exptime'] = '3600'

    d1 = CCDData(np.ones((5, 5)), meta=hdr, unit=u.electron)
    assert d1.meta['OBSERVER'] == 'Edwin Hubble'
    assert d1.header['OBSERVER'] == 'Edwin Hubble'


def test_metafromdict():
    dic = {'OBSERVER': 'Edwin Hubble', 'EXPTIME': 3600}
    d1 = CCDData(np.ones((5, 5)), meta=dic, unit=u.electron)
    assert d1.meta['OBSERVER'] == 'Edwin Hubble'


def test_header2meta():
    hdr = fits.header.Header()
    hdr['observer'] = 'Edwin Hubble'
    hdr['exptime'] = '3600'

    d1 = CCDData(np.ones((5, 5)), unit=u.electron)
    d1.header = hdr
    assert d1.meta['OBSERVER'] == 'Edwin Hubble'
    assert d1.header['OBSERVER'] == 'Edwin Hubble'


def test_metafromstring_fail():
    hdr = 'this is not a valid header'
    with pytest.raises(TypeError):
        CCDData(np.ones((5, 5)), meta=hdr, unit=u.adu)


def test_setting_bad_uncertainty_raises_error():
    ccd_data = create_ccd_data()
    with pytest.raises(TypeError):
        # Uncertainty is supposed to be an instance of NDUncertainty
        ccd_data.uncertainty = 10


def test_setting_uncertainty_with_array():
    ccd_data = create_ccd_data()
    ccd_data.uncertainty = None
    fake_uncertainty = np.sqrt(np.abs(ccd_data.data))
    ccd_data.uncertainty = fake_uncertainty.copy()
    np.testing.assert_array_equal(ccd_data.uncertainty.array, fake_uncertainty)


def test_setting_uncertainty_wrong_shape_raises_error():
    ccd_data = create_ccd_data()
    with pytest.raises(ValueError):
        ccd_data.uncertainty = np.zeros([3, 4])


def test_to_hdu():
    ccd_data = create_ccd_data()
    ccd_data.meta = {'observer': 'Edwin Hubble'}
    fits_hdulist = ccd_data.to_hdu()
    assert isinstance(fits_hdulist, fits.HDUList)
    for k, v in ccd_data.meta.items():
        assert fits_hdulist[0].header[k] == v
    np.testing.assert_array_equal(fits_hdulist[0].data, ccd_data.data)


def test_copy():
    ccd_data = create_ccd_data()
    ccd_copy = ccd_data.copy()
    np.testing.assert_array_equal(ccd_copy.data, ccd_data.data)
    assert ccd_copy.unit == ccd_data.unit
    assert ccd_copy.meta == ccd_data.meta


@pytest.mark.parametrize('operation,affects_uncertainty', [
                         ("multiply", True),
                         ("divide", True),
                         ])
@pytest.mark.parametrize('operand', [
                         2.0,
                         2 * u.dimensionless_unscaled,
                         2 * u.photon / u.adu,
                         ])
@pytest.mark.parametrize('with_uncertainty', [
                         True,
                         False])
def test_mult_div_overload(operand, with_uncertainty,
                           operation, affects_uncertainty):
    ccd_data = create_ccd_data()
    if with_uncertainty:
        ccd_data.uncertainty = StdDevUncertainty(np.ones_like(ccd_data))
    method = getattr(ccd_data, operation)
    np_method = getattr(np, operation)
    result = method(operand)
    assert result is not ccd_data
    assert isinstance(result, CCDData)
    assert (result.uncertainty is None or
            isinstance(result.uncertainty, StdDevUncertainty))
    try:
        op_value = operand.value
    except AttributeError:
        op_value = operand

    np.testing.assert_array_equal(result.data,
                                  np_method(ccd_data.data, op_value))
    if with_uncertainty:
        if affects_uncertainty:
            np.testing.assert_array_equal(result.uncertainty.array,
                                          np_method(ccd_data.uncertainty.array,
                                                    op_value))
        else:
            np.testing.assert_array_equal(result.uncertainty.array,
                                          ccd_data.uncertainty.array)
    else:
        assert result.uncertainty is None

    if isinstance(operand, u.Quantity):
        # Need the "1 *" below to force arguments to be Quantity to work around
        # astropy/astropy#2377
        expected_unit = np_method(1 * ccd_data.unit, 1 * operand.unit).unit
        assert result.unit == expected_unit
    else:
        assert result.unit == ccd_data.unit


@pytest.mark.parametrize('operation,affects_uncertainty', [
                         ("add", False),
                         ("subtract", False),
                         ])
@pytest.mark.parametrize('operand,expect_failure', [
                         (2.0, u.UnitsError),  # fail--units don't match image
                         (2 * u.dimensionless_unscaled, u.UnitsError),  # same
                         (2 * u.adu, False),
                         ])
@pytest.mark.parametrize('with_uncertainty', [
                         True,
                         False])
def test_add_sub_overload(operand, expect_failure, with_uncertainty,
                          operation, affects_uncertainty):
    ccd_data = create_ccd_data()
    if with_uncertainty:
        ccd_data.uncertainty = StdDevUncertainty(np.ones_like(ccd_data))
    method = getattr(ccd_data, operation)
    np_method = getattr(np, operation)
    if expect_failure:
        with pytest.raises(expect_failure):
            result = method(operand)
        return
    else:
        result = method(operand)
    assert result is not ccd_data
    assert isinstance(result, CCDData)
    assert (result.uncertainty is None or
            isinstance(result.uncertainty, StdDevUncertainty))
    try:
        op_value = operand.value
    except AttributeError:
        op_value = operand

    np.testing.assert_array_equal(result.data,
                                  np_method(ccd_data.data, op_value))
    if with_uncertainty:
        if affects_uncertainty:
            np.testing.assert_array_equal(result.uncertainty.array,
                                          np_method(ccd_data.uncertainty.array,
                                                    op_value))
        else:
            np.testing.assert_array_equal(result.uncertainty.array,
                                          ccd_data.uncertainty.array)
    else:
        assert result.uncertainty is None

    if isinstance(operand, u.Quantity):
        assert (result.unit == ccd_data.unit and result.unit == operand.unit)
    else:
        assert result.unit == ccd_data.unit


def test_arithmetic_overload_fails():
    ccd_data = create_ccd_data()
    with pytest.raises(TypeError):
        ccd_data.multiply("five")

    with pytest.raises(TypeError):
        ccd_data.divide("five")

    with pytest.raises(TypeError):
        ccd_data.add("five")

    with pytest.raises(TypeError):
        ccd_data.subtract("five")


def test_arithmetic_no_wcs_compare():
    ccd = CCDData(np.ones((10, 10)), unit='')
    assert ccd.add(ccd, compare_wcs=None).wcs is None
    assert ccd.subtract(ccd, compare_wcs=None).wcs is None
    assert ccd.multiply(ccd, compare_wcs=None).wcs is None
    assert ccd.divide(ccd, compare_wcs=None).wcs is None


def test_arithmetic_with_wcs_compare():
    def return_true(_, __):
        return True

    wcs1, wcs2 = nd_testing.create_two_equal_wcs(naxis=2)
    ccd1 = CCDData(np.ones((10, 10)), unit='', wcs=wcs1)
    ccd2 = CCDData(np.ones((10, 10)), unit='', wcs=wcs2)
    nd_testing.assert_wcs_seem_equal(
        ccd1.add(ccd2, compare_wcs=return_true).wcs,
        wcs1)
    nd_testing.assert_wcs_seem_equal(
        ccd1.subtract(ccd2, compare_wcs=return_true).wcs,
        wcs1)
    nd_testing.assert_wcs_seem_equal(
        ccd1.multiply(ccd2, compare_wcs=return_true).wcs,
        wcs1)
    nd_testing.assert_wcs_seem_equal(
        ccd1.divide(ccd2, compare_wcs=return_true).wcs,
        wcs1)


def test_arithmetic_with_wcs_compare_fail():
    def return_false(_, __):
        return False

    ccd1 = CCDData(np.ones((10, 10)), unit='', wcs=WCS())
    ccd2 = CCDData(np.ones((10, 10)), unit='', wcs=WCS())
    with pytest.raises(ValueError):
        ccd1.add(ccd2, compare_wcs=return_false)
    with pytest.raises(ValueError):
        ccd1.subtract(ccd2, compare_wcs=return_false)
    with pytest.raises(ValueError):
        ccd1.multiply(ccd2, compare_wcs=return_false)
    with pytest.raises(ValueError):
        ccd1.divide(ccd2, compare_wcs=return_false)


def test_arithmetic_overload_ccddata_operand():
    ccd_data = create_ccd_data()
    ccd_data.uncertainty = StdDevUncertainty(np.ones_like(ccd_data))
    operand = ccd_data.copy()
    result = ccd_data.add(operand)
    assert len(result.meta) == 0
    np.testing.assert_array_equal(result.data,
                                  2 * ccd_data.data)
    np.testing.assert_array_almost_equal_nulp(
        result.uncertainty.array,
        np.sqrt(2) * ccd_data.uncertainty.array
    )

    result = ccd_data.subtract(operand)
    assert len(result.meta) == 0
    np.testing.assert_array_equal(result.data,
                                  0 * ccd_data.data)
    np.testing.assert_array_almost_equal_nulp(
        result.uncertainty.array,
        np.sqrt(2) * ccd_data.uncertainty.array
    )

    result = ccd_data.multiply(operand)
    assert len(result.meta) == 0
    np.testing.assert_array_equal(result.data,
                                  ccd_data.data ** 2)
    expected_uncertainty = (np.sqrt(2) * np.abs(ccd_data.data) *
                            ccd_data.uncertainty.array)
    np.testing.assert_allclose(result.uncertainty.array,
                               expected_uncertainty)

    result = ccd_data.divide(operand)
    assert len(result.meta) == 0
    np.testing.assert_array_equal(result.data,
                                  np.ones_like(ccd_data.data))
    expected_uncertainty = (np.sqrt(2) / np.abs(ccd_data.data) *
                            ccd_data.uncertainty.array)
    np.testing.assert_allclose(result.uncertainty.array,
                               expected_uncertainty)


def test_arithmetic_overload_differing_units():
    a = np.array([1, 2, 3]) * u.m
    b = np.array([1, 2, 3]) * u.cm
    ccddata = CCDData(a)

    # TODO: Could also be parametrized.
    res = ccddata.add(b)
    np.testing.assert_array_almost_equal(res.data, np.add(a, b).value)
    assert res.unit == np.add(a, b).unit

    res = ccddata.subtract(b)
    np.testing.assert_array_almost_equal(res.data, np.subtract(a, b).value)
    assert res.unit == np.subtract(a, b).unit

    res = ccddata.multiply(b)
    np.testing.assert_array_almost_equal(res.data, np.multiply(a, b).value)
    assert res.unit == np.multiply(a, b).unit

    res = ccddata.divide(b)
    np.testing.assert_array_almost_equal(res.data, np.divide(a, b).value)
    assert res.unit == np.divide(a, b).unit


def test_arithmetic_add_with_array():
    ccd = CCDData(np.ones((3, 3)), unit='')
    res = ccd.add(np.arange(3))
    np.testing.assert_array_equal(res.data, [[1, 2, 3]] * 3)

    ccd = CCDData(np.ones((3, 3)), unit='adu')
    with pytest.raises(ValueError):
        ccd.add(np.arange(3))


def test_arithmetic_subtract_with_array():
    ccd = CCDData(np.ones((3, 3)), unit='')
    res = ccd.subtract(np.arange(3))
    np.testing.assert_array_equal(res.data, [[1, 0, -1]] * 3)

    ccd = CCDData(np.ones((3, 3)), unit='adu')
    with pytest.raises(ValueError):
        ccd.subtract(np.arange(3))


def test_arithmetic_multiply_with_array():
    ccd = CCDData(np.ones((3, 3)) * 3, unit=u.m)
    res = ccd.multiply(np.ones((3, 3)) * 2)
    np.testing.assert_array_equal(res.data, [[6, 6, 6]] * 3)
    assert res.unit == ccd.unit


def test_arithmetic_divide_with_array():
    ccd = CCDData(np.ones((3, 3)), unit=u.m)
    res = ccd.divide(np.ones((3, 3)) * 2)
    np.testing.assert_array_equal(res.data, [[0.5, 0.5, 0.5]] * 3)
    assert res.unit == ccd.unit


def test_history_preserved_if_metadata_is_fits_header(tmpdir):
    fake_img = np.zeros([2, 2])
    hdu = fits.PrimaryHDU(fake_img)
    hdu.header['history'] = 'one'
    hdu.header['history'] = 'two'
    hdu.header['history'] = 'three'
    assert len(hdu.header['history']) == 3
    tmp_file = tmpdir.join('temp.fits').strpath
    hdu.writeto(tmp_file)

    ccd_read = CCDData.read(tmp_file, unit="adu")
    assert ccd_read.header['history'] == hdu.header['history']


def test_infol_logged_if_unit_in_fits_header(tmpdir):
    ccd_data = create_ccd_data()
    tmpfile = tmpdir.join('temp.fits')
    ccd_data.write(tmpfile.strpath)
    log.setLevel('INFO')
    explicit_unit_name = "photon"
    with log.log_to_list() as log_list:
        ccd_from_disk = CCDData.read(tmpfile.strpath, unit=explicit_unit_name)
        assert explicit_unit_name in log_list[0].message


def test_wcs_attribute(tmpdir):
    """
    Check that WCS attribute gets added to header, and that if a CCDData
    object is created from a FITS file with a header, and the WCS attribute
    is modified, then the CCDData object is turned back into an hdu, the
    WCS object overwrites the old WCS information in the header.
    """
    ccd_data = create_ccd_data()
    tmpfile = tmpdir.join('temp.fits')
    # This wcs example is taken from the astropy.wcs docs.
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = np.array(ccd_data.shape) / 2
    wcs.wcs.cdelt = np.array([-0.066667, 0.066667])
    wcs.wcs.crval = [0, -90]
    wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    wcs.wcs.set_pv([(2, 1, 45.0)])
    ccd_data.header = ccd_data.to_hdu()[0].header
    ccd_data.header.extend(wcs.to_header(), useblanks=False)
    ccd_data.write(tmpfile.strpath)

    # Get the header length after it has been extended by the WCS keywords
    original_header_length = len(ccd_data.header)

    ccd_new = CCDData.read(tmpfile.strpath)
    # WCS attribute should be set for ccd_new
    assert ccd_new.wcs is not None
    # WCS attribute should be equal to wcs above.
    assert ccd_new.wcs.wcs == wcs.wcs

    # Converting CCDData object with wcs to an hdu shouldn't
    # create duplicate wcs-related entries in the header.
    ccd_new_hdu = ccd_new.to_hdu()[0]
    assert len(ccd_new_hdu.header) == original_header_length

    # Making a CCDData with WCS (but not WCS in the header) should lead to
    # WCS information in the header when it is converted to an HDU.
    ccd_wcs_not_in_header = CCDData(ccd_data.data, wcs=wcs, unit="adu")
    hdu = ccd_wcs_not_in_header.to_hdu()[0]
    wcs_header = wcs.to_header()
    for k in wcs_header.keys():
        # Skip these keywords if they are in the WCS header because they are
        # not WCS-specific.
        if k in ['', 'COMMENT', 'HISTORY']:
            continue
        # No keyword from the WCS should be in the header.
        assert k not in ccd_wcs_not_in_header.header
        # Every keyword in the WCS should be in the header of the HDU
        assert hdu.header[k] == wcs_header[k]

    # Now check that if WCS of a CCDData is modified, then the CCDData is
    # converted to an HDU, the WCS keywords in the header are overwritten
    # with the appropriate keywords from the header.
    #
    # ccd_new has a WCS and WCS keywords in the header, so try modifying
    # the WCS.
    ccd_new.wcs.wcs.cdelt *= 2
    ccd_new_hdu_mod_wcs = ccd_new.to_hdu()[0]
    assert ccd_new_hdu_mod_wcs.header['CDELT1'] == ccd_new.wcs.wcs.cdelt[0]
    assert ccd_new_hdu_mod_wcs.header['CDELT2'] == ccd_new.wcs.wcs.cdelt[1]


def test_wcs_keywords_removed_from_header():
    """
    Test, for the file included with the nddata tests, that WCS keywords are
    properly removed from header.
    """
    from astropy.nddata.ccddata import _KEEP_THESE_KEYWORDS_IN_HEADER
    keepers = set(_KEEP_THESE_KEYWORDS_IN_HEADER)
    data_file = get_pkg_data_filename('data/sip-wcs.fits')
    ccd = CCDData.read(data_file)
    with pytest.warns(AstropyWarning,
                      match=r'Some non-standard WCS keywords were excluded'):
        wcs_header = ccd.wcs.to_header()
    assert not (set(wcs_header) & set(ccd.meta) - keepers)

    # Make sure that exceptions are not raised when trying to remove missing
    # keywords. o4sp040b0_raw.fits of io.fits is missing keyword 'PC1_1'.
    data_file1 = get_pkg_data_filename('../../io/fits/tests/data/o4sp040b0_raw.fits')
    with pytest.warns(FITSFixedWarning, match=r"'unitfix' made the change"):
        ccd = CCDData.read(data_file1, unit='count')


def test_wcs_SIP_coefficient_keywords_removed():
    # If SIP polynomials are present, check that no more polynomial
    # coefficients remain in the header. See #8598

    # The SIP paper is ambiguous as to whether keywords like
    # A_0_0 can appear in the header for a 2nd order or higher
    # polynomial. The paper clearly says that the corrections
    # are only for quadratic or higher order, so A_0_0 and the like
    # should be zero if they are present, but they apparently can be
    # there (or at least astrometry.net produces them).

    # astropy WCS does not write those coefficients, so they were
    # not being removed from the header even though they are WCS-related.

    data_file = get_pkg_data_filename('data/sip-wcs.fits')
    test_keys = ['A_0_0', 'B_0_1']

    # Make sure the keywords added to this file for testing are there
    with fits.open(data_file) as hdu:
        for key in test_keys:
            assert key in hdu[0].header

    ccd = CCDData.read(data_file)

    # Now the test...the two keywords above should have been removed.
    for key in test_keys:
        assert key not in ccd.header


@pytest.mark.filterwarnings('ignore')
def test_wcs_keyword_removal_for_wcs_test_files():
    """
    Test, for the WCS test files, that keyword removal works as
    expected. Those cover a much broader range of WCS types than
    test_wcs_keywords_removed_from_header.

    Includes regression test for #8597
    """
    from astropy.nddata.ccddata import _generate_wcs_and_update_header
    from astropy.nddata.ccddata import (_KEEP_THESE_KEYWORDS_IN_HEADER,
                                        _CDs, _PCs)

    keepers = set(_KEEP_THESE_KEYWORDS_IN_HEADER)
    wcs_headers = get_pkg_data_filenames('../../wcs/tests/data',
                                         pattern='*.hdr')

    for hdr in wcs_headers:
        # Skip the files that are expected to be bad...
        if 'invalid' in hdr or 'nonstandard' in hdr or 'segfault' in hdr:
            continue
        header_string = get_pkg_data_contents(hdr)
        header = fits.Header.fromstring(header_string)

        wcs = WCS(header_string)
        header_from_wcs = wcs.to_header(relax=True)

        new_header, new_wcs = _generate_wcs_and_update_header(header)
        new_wcs_header = new_wcs.to_header(relax=True)

        # Make sure all of the WCS-related keywords generated by astropy
        # have been removed.
        assert not (set(new_header) &
                    set(new_wcs_header) -
                    keepers)

        # Check that new_header contains no remaining WCS information.
        # Specifically, check that
        # 1. The combination of new_header and new_wcs does not contain
        #    both PCi_j and CDi_j keywords. See #8597.

        # Check for 1
        final_header = new_header + new_wcs_header
        final_header_set = set(final_header)

        if _PCs & final_header_set:
            assert not (_CDs & final_header_set)
        elif _CDs & final_header_set:
            assert not (_PCs & final_header_set)

        # Check that the new wcs is the same as the old.
        for k, v in new_wcs_header.items():
            if isinstance(v, str):
                assert header_from_wcs[k] == v
            else:
                np.testing.assert_almost_equal(header_from_wcs[k], v)


def test_read_wcs_not_creatable(tmpdir):
    # The following Header can't be converted to a WCS object. See also #6499.
    hdr_txt_example_WCS = textwrap.dedent('''
    SIMPLE  =                    T / Fits standard
    BITPIX  =                   16 / Bits per pixel
    NAXIS   =                    2 / Number of axes
    NAXIS1  =                 1104 / Axis length
    NAXIS2  =                 4241 / Axis length
    CRVAL1  =         164.98110962 / Physical value of the reference pixel X
    CRVAL2  =          44.34089279 / Physical value of the reference pixel Y
    CRPIX1  =                -34.0 / Reference pixel in X (pixel)
    CRPIX2  =               2041.0 / Reference pixel in Y (pixel)
    CDELT1  =           0.10380000 / X Scale projected on detector (#/pix)
    CDELT2  =           0.10380000 / Y Scale projected on detector (#/pix)
    CTYPE1  = 'RA---TAN'           / Pixel coordinate system
    CTYPE2  = 'WAVELENGTH'         / Pixel coordinate system
    CUNIT1  = 'degree  '           / Units used in both CRVAL1 and CDELT1
    CUNIT2  = 'nm      '           / Units used in both CRVAL2 and CDELT2
    CD1_1   =           0.20760000 / Pixel Coordinate translation matrix
    CD1_2   =           0.00000000 / Pixel Coordinate translation matrix
    CD2_1   =           0.00000000 / Pixel Coordinate translation matrix
    CD2_2   =           0.10380000 / Pixel Coordinate translation matrix
    C2YPE1  = 'RA---TAN'           / Pixel coordinate system
    C2YPE2  = 'DEC--TAN'           / Pixel coordinate system
    C2NIT1  = 'degree  '           / Units used in both C2VAL1 and C2ELT1
    C2NIT2  = 'degree  '           / Units used in both C2VAL2 and C2ELT2
    RADECSYS= 'FK5     '           / The equatorial coordinate system
    ''')
    with catch_warnings(FITSFixedWarning):
        hdr = fits.Header.fromstring(hdr_txt_example_WCS, sep='\n')
    hdul = fits.HDUList([fits.PrimaryHDU(np.ones((4241, 1104)), header=hdr)])
    filename = tmpdir.join('afile.fits').strpath
    hdul.writeto(filename)
    # The hdr cannot be converted to a WCS object because of an
    # InconsistentAxisTypesError but it should still open the file
    ccd = CCDData.read(filename, unit='adu')
    assert ccd.wcs is None


def test_header():
    ccd_data = create_ccd_data()
    a = {'Observer': 'Hubble'}
    ccd = CCDData(ccd_data, header=a)
    assert ccd.meta == a


def test_wcs_arithmetic():
    ccd_data = create_ccd_data()
    wcs = WCS(naxis=2)
    ccd_data.wcs = wcs
    result = ccd_data.multiply(1.0)
    nd_testing.assert_wcs_seem_equal(result.wcs, wcs)


@pytest.mark.parametrize('operation',
                         ['multiply', 'divide', 'add', 'subtract'])
def test_wcs_arithmetic_ccd(operation):
    ccd_data = create_ccd_data()
    ccd_data2 = ccd_data.copy()
    ccd_data.wcs = WCS(naxis=2)
    method = getattr(ccd_data, operation)
    result = method(ccd_data2)
    nd_testing.assert_wcs_seem_equal(result.wcs, ccd_data.wcs)
    assert ccd_data2.wcs is None


def test_wcs_sip_handling():
    """
    Check whether the ctypes RA---TAN-SIP and DEC--TAN-SIP survive
    a roundtrip unchanged.
    """
    data_file = get_pkg_data_filename('data/sip-wcs.fits')

    def check_wcs_ctypes(header):
        expected_wcs_ctypes = {
            'CTYPE1': 'RA---TAN-SIP',
            'CTYPE2': 'DEC--TAN-SIP'
        }

        return [header[k] == v for k, v in expected_wcs_ctypes.items()]

    ccd_original = CCDData.read(data_file)
    # After initialization the keywords should be in the WCS, not in the
    # meta.
    with fits.open(data_file) as raw:
        good_ctype = check_wcs_ctypes(raw[0].header)
    assert all(good_ctype)

    ccd_new = ccd_original.to_hdu()
    good_ctype = check_wcs_ctypes(ccd_new[0].header)
    assert all(good_ctype)

    # Try converting to header with wcs_relax=False and
    # the header should contain the CTYPE keywords without
    # the -SIP

    ccd_no_relax = ccd_original.to_hdu(wcs_relax=False)
    good_ctype = check_wcs_ctypes(ccd_no_relax[0].header)
    assert not any(good_ctype)
    assert ccd_no_relax[0].header['CTYPE1'] == 'RA---TAN'
    assert ccd_no_relax[0].header['CTYPE2'] == 'DEC--TAN'


@pytest.mark.parametrize('operation',
                         ['multiply', 'divide', 'add', 'subtract'])
def test_mask_arithmetic_ccd(operation):
    ccd_data = create_ccd_data()
    ccd_data2 = ccd_data.copy()
    ccd_data.mask = (ccd_data.data > 0)
    method = getattr(ccd_data, operation)
    result = method(ccd_data2)
    np.testing.assert_equal(result.mask, ccd_data.mask)


def test_write_read_multiextensionfits_mask_default(tmpdir):
    # Test that if a mask is present the mask is saved and loaded by default.
    ccd_data = create_ccd_data()
    ccd_data.mask = ccd_data.data > 10
    filename = tmpdir.join('afile.fits').strpath
    ccd_data.write(filename)
    ccd_after = CCDData.read(filename)
    assert ccd_after.mask is not None
    np.testing.assert_array_equal(ccd_data.mask, ccd_after.mask)


@pytest.mark.parametrize(
    'uncertainty_type',
    [StdDevUncertainty, VarianceUncertainty, InverseVariance])
def test_write_read_multiextensionfits_uncertainty_default(
        tmpdir, uncertainty_type):
    # Test that if a uncertainty is present it is saved and loaded by default.
    ccd_data = create_ccd_data()
    ccd_data.uncertainty = uncertainty_type(ccd_data.data * 10)
    filename = tmpdir.join('afile.fits').strpath
    ccd_data.write(filename)
    ccd_after = CCDData.read(filename)
    assert ccd_after.uncertainty is not None
    assert type(ccd_after.uncertainty) is uncertainty_type
    np.testing.assert_array_equal(ccd_data.uncertainty.array,
                                  ccd_after.uncertainty.array)


@pytest.mark.parametrize(
    'uncertainty_type',
    [StdDevUncertainty, VarianceUncertainty, InverseVariance])
def test_write_read_multiextensionfits_uncertainty_different_uncertainty_key(
        tmpdir, uncertainty_type):
    # Test that if a uncertainty is present it is saved and loaded by default.
    ccd_data = create_ccd_data()
    ccd_data.uncertainty = uncertainty_type(ccd_data.data * 10)
    filename = tmpdir.join('afile.fits').strpath
    ccd_data.write(filename, key_uncertainty_type='Blah')
    ccd_after = CCDData.read(filename, key_uncertainty_type='Blah')
    assert ccd_after.uncertainty is not None
    assert type(ccd_after.uncertainty) is uncertainty_type
    np.testing.assert_array_equal(ccd_data.uncertainty.array,
                                  ccd_after.uncertainty.array)


def test_write_read_multiextensionfits_not(tmpdir):
    # Test that writing mask and uncertainty can be disabled
    ccd_data = create_ccd_data()
    ccd_data.mask = ccd_data.data > 10
    ccd_data.uncertainty = StdDevUncertainty(ccd_data.data * 10)
    filename = tmpdir.join('afile.fits').strpath
    ccd_data.write(filename, hdu_mask=None, hdu_uncertainty=None)
    ccd_after = CCDData.read(filename)
    assert ccd_after.uncertainty is None
    assert ccd_after.mask is None


def test_write_read_multiextensionfits_custom_ext_names(tmpdir):
    # Test writing mask, uncertainty in another extension than default
    ccd_data = create_ccd_data()
    ccd_data.mask = ccd_data.data > 10
    ccd_data.uncertainty = StdDevUncertainty(ccd_data.data * 10)
    filename = tmpdir.join('afile.fits').strpath
    ccd_data.write(filename, hdu_mask='Fun', hdu_uncertainty='NoFun')

    # Try reading with defaults extension names
    ccd_after = CCDData.read(filename)
    assert ccd_after.uncertainty is None
    assert ccd_after.mask is None

    # Try reading with custom extension names
    ccd_after = CCDData.read(filename, hdu_mask='Fun', hdu_uncertainty='NoFun')
    assert ccd_after.uncertainty is not None
    assert ccd_after.mask is not None
    np.testing.assert_array_equal(ccd_data.mask, ccd_after.mask)
    np.testing.assert_array_equal(ccd_data.uncertainty.array,
                                  ccd_after.uncertainty.array)


def test_read_old_style_multiextensionfits(tmpdir):
    # Regression test for https://github.com/astropy/ccdproc/issues/664
    #
    # Prior to astropy 3.1 there was no uncertainty type saved
    # in the multiextension fits files generated by CCDData
    # because the uncertainty had to be StandardDevUncertainty.
    #
    # Current version should be able to read those in.
    #
    size = 4
    # Value of the variables below are not important to the test.
    data = np.zeros([size, size])
    mask = data > 0.9
    uncert = np.sqrt(data)

    ccd = CCDData(data=data, mask=mask, uncertainty=uncert, unit='adu')
    # We'll create the file manually to ensure we have the
    # right extension names and no uncertainty type.
    hdulist = ccd.to_hdu()
    del hdulist[2].header['UTYPE']
    file_name = tmpdir.join('old_ccddata_mef.fits').strpath
    hdulist.writeto(file_name)

    ccd = CCDData.read(file_name)

    assert isinstance(ccd.uncertainty, StdDevUncertainty)


def test_wcs():
    ccd_data = create_ccd_data()
    wcs = WCS(naxis=2)
    ccd_data.wcs = wcs
    assert ccd_data.wcs is wcs


def test_recognized_fits_formats_for_read_write(tmpdir):
    # These are the extensions that are supposed to be supported.
    ccd_data = create_ccd_data()
    supported_extensions = ['fit', 'fits', 'fts']

    for ext in supported_extensions:
        path = tmpdir.join(f"test.{ext}")
        ccd_data.write(path.strpath)
        from_disk = CCDData.read(path.strpath)
        assert (ccd_data.data == from_disk.data).all()


def test_stddevuncertainty_compat_descriptor_no_parent():
    with pytest.raises(MissingDataAssociationException):
        StdDevUncertainty(np.ones((10, 10))).parent_nddata


def test_stddevuncertainty_compat_descriptor_no_weakref():
    # TODO: Remove this test if astropy 1.0 isn't supported anymore
    # This test might create a Memoryleak on purpose, so the last lines after
    # the assert are IMPORTANT cleanup.
    ccd = CCDData(np.ones((10, 10)), unit='')
    uncert = StdDevUncertainty(np.ones((10, 10)))
    uncert._parent_nddata = ccd
    assert uncert.parent_nddata is ccd
    uncert._parent_nddata = None


# https://github.com/astropy/astropy/issues/7595
def test_read_returns_image(tmpdir):
    # Test if CCData.read returns a image when reading a fits file containing
    # a table and image, in that order.
    tbl = Table(np.ones(10).reshape(5, 2))
    img = np.ones((5, 5))
    hdul = fits.HDUList(hdus=[fits.PrimaryHDU(), fits.TableHDU(tbl.as_array()),
                              fits.ImageHDU(img)])
    filename = tmpdir.join('table_image.fits').strpath
    hdul.writeto(filename)
    ccd = CCDData.read(filename, unit='adu')
    # Expecting to get (5, 5), the size of the image
    assert ccd.data.shape == (5, 5)


# https://github.com/astropy/astropy/issues/9664
def test_sliced_ccdata_to_hdu():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = 10, 10
    ccd = CCDData(np.ones((10, 10)), wcs=wcs, unit='pixel')
    trimmed = ccd[2:-2, 2:-2]
    hdul = trimmed.to_hdu()
    assert isinstance(hdul, fits.HDUList)
    assert hdul[0].header['CRPIX1'] == 8
    assert hdul[0].header['CRPIX2'] == 8
