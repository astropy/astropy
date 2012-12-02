import glob
import os
import sys
import warnings

import numpy as np
from numpy.testing import (
    assert_array_almost_equal, assert_array_almost_equal_nulp)

from ...tests.helper import raises
from ... import wcs
from ...utils.data import (
    get_pkg_data_filenames, get_pkg_data_contents, get_pkg_data_filename)


# test_maps() is a generator
def test_maps():

    # test_map() is the function that is called to perform the generated test
    def test_map(filename):

        # the test parameter is the base name of the file to use; find
        # the file in the installed wcs test directory
        header = get_pkg_data_contents(
            os.path.join("maps", filename), encoding='binary')
        wcsobj = wcs.WCS(header)

        world = wcsobj.wcs_pix2world([[97, 97]], 1)

        assert_array_almost_equal(world, [[285.0, -66.25]], decimal=1)

        pix = wcsobj.wcs_world2pix([[285.0, -66.25]], 1)

        assert_array_almost_equal(pix, [[97, 97]], decimal=0)

    # get the list of the hdr files that we want to test
    hdr_file_list = list(get_pkg_data_filenames("maps", "*.hdr"))

    # actually perform a test for each one
    for filename in hdr_file_list:

        # use the base name of the file, because everything we yield
        # will show up in the test name in the pandokia report
        filename = os.path.basename(filename)

        # yield a function name and parameters to make a generated test
        yield test_map, filename

    # AFTER we tested with every file that we found, check to see that we
    # actually have the list we expect.  If N=0, we will not have performed
    # any tests at all.  If N < n_data_files, we are missing some files,
    # so we will have skipped some tests.  Without this check, both cases
    # happen silently!

    # how many do we expect to see?
    n_data_files = 28

    if len(hdr_file_list) != n_data_files:
        assert False, (
            "test_maps has wrong number data files: found %d, expected "
            " %d" % (len(hdr_file_list), n_data_files))
        # b.t.w.  If this assert happens, py.test reports one more test
        # than it would have otherwise.


# test_spectra() is a generator
def test_spectra():

    # test_spectrum() is the function that is called to perform the
    # generated test
    def test_spectrum(filename):

        # the test parameter is the base name of the file to use; find
        # the file in the installed wcs test directory
        header = get_pkg_data_contents(
            os.path.join("spectra", filename), encoding='binary')

        wcsobj = wcs.WCS(header)

        all = wcs.find_all_wcs(header)
        assert len(all) == 9

    # get the list of the hdr files that we want to test
    hdr_file_list = list(get_pkg_data_filenames("spectra", "*.hdr"))

    # actually perform a test for each one
    for filename in hdr_file_list:

        # use the base name of the file, because everything we yield
        # will show up in the test name in the pandokia report
        filename = os.path.basename(filename)

        # yield a function name and parameters to make a generated test
        yield test_spectrum, filename

    # AFTER we tested with every file that we found, check to see that we
    # actually have the list we expect.  If N=0, we will not have performed
    # any tests at all.  If N < n_data_files, we are missing some files,
    # so we will have skipped some tests.  Without this check, both cases
    # happen silently!

    # how many do we expect to see?
    n_data_files = 6

    if len(hdr_file_list) != n_data_files:
        assert False, (
            "test_spectra has wrong number data files: found %d, expected "
            " %d" % (len(hdr_file_list), n_data_files))
        # b.t.w.  If this assert happens, py.test reports one more test
        # than it would have otherwise.


def test_units():
    u = wcs.UnitConverter("log(MHz)", "ln(Hz)")
    print(u.convert([1, 2, 3, 4]))

basic_units = "m s g rad sr K A mol cd".split()
derived_units = "Hz J W V N Pa C Ohm ohm S F Wb T H lm lx".split()
add_all_units = "eV Jy R G barn".split()
add_sup_units = "a yr pc bit byte Byte".split()
add_sub_units = "mag".split()
general_units = (
    "deg arcmin arcsec mas d h min erg Ry u D DEGREE DEGREES".split())
astro_units = "Angstrom angstrom AU lyr beam solRad solMass solLum Sun".split()
device_units = "adu bin chan count ct photon ph pixel pix voxel".split()

sub_prefixes = "y z a f p n u m c d".split()
sup_prefixes = "da h k M G T P E Z Y".split()


def test_all_units():
    def test_self(x):
        # x appears in the test name. If we would have had an ambiguous
        # test name, we had -xxx added to the unit name.  Remove it if
        # necessary.
        if '-' in x:
            x = x.split('-')[0]

        # here is the test:
        try:
            u = wcs.UnitConverter(x, x)
        except ValueError:
            e = sys.exc_info()[1]
            if str(e).startswith("ERROR 12 in wcsutrne") and \
                    x in ("S", "H", "D"):
                return
            else:
                raise
        assert u.scale == 1.0
        assert u.offset == 0.0
        assert u.power == 1.0

    # list of all the units to test
    all = sorted(basic_units + derived_units + add_all_units + add_sup_units
                 + add_sub_units + general_units + astro_units + device_units)

    # Pandokia has non-case-sensitve test names; since the unit name is
    # showing up in the test name, we want to disambiguate any name collisions.
    # Here is a list of all the lower-cased unit name names.
    all_lower = [x.lower() for x in all]

    # here are serial numbers to use to disambiguate
    unique_tags = {}

    for unit in all:
        # disambiguate the test name, if necessary
        l_unit = unit.lower()
        if unit != l_unit and l_unit in all_lower:
            n = unique_tags.get(l_unit, 1)
            unique_tags[n] = n + 1

            # the test will tear off the part after the '-'
            unit = '%s-%d' % (unit, n)

        # perform the test
        yield test_self, unit


def test_unit_prefixes():
    def test_self(x, p):
        unit = p + x
        try:
            u = wcs.UnitConverter(unit, unit)
        except ValueError:
            e = sys.exc_info()[1]
            if str(e) == "Potentially unsafe translation" and \
                    x in ("S", "H", "D"):
                return
            else:
                raise
        assert u.scale == 1.0
        assert u.offset == 0.0
        assert u.power == 1.0

    for unit in (basic_units + derived_units + add_all_units):
        for prefix in (sub_prefixes + sup_prefixes):
            yield test_self, unit, prefix

    for unit in add_sup_units:
        for prefix in sup_prefixes:
            yield test_self, unit, prefix

    for unit in add_sub_units:
        for prefix in sub_prefixes:
            yield test_self, unit, prefix


def test_fixes():
    """
    From github issue #36
    """
    def run():
        header = get_pkg_data_contents(
            'data/nonstandard_units.hdr', encoding='binary')
        w = wcs.WCS(header)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        run()
        assert len(w) == 2
        for item in w:
            assert issubclass(item.category, wcs.FITSFixedWarning)
            if 'unitfix' in str(item.message):
                assert 'Hz' in str(item.message)


def test_outside_sky():
    """
    From github issue #107
    """
    header = get_pkg_data_contents(
        'data/outside_sky.hdr', encoding='binary')
    w = wcs.WCS(header)

    assert np.all(np.isnan(w.wcs_pix2world([[100., 500.]], 0)))  # outside sky
    assert np.all(np.isnan(w.wcs_pix2world([[200., 200.]], 0)))  # outside sky
    assert not np.any(np.isnan(w.wcs_pix2world([[1000., 1000.]], 0)))


def test_load_fits_path():
    fits = get_pkg_data_filename('data/sip.fits')
    w = wcs.WCS(fits)


def test_backward_compatible():
    fits = get_pkg_data_filename('data/sip.fits')
    w = wcs.WCS(fits)

    data = np.random.rand(100, 2)
    assert np.all(w.wcs_pix2world(data, 0) == w.wcs_pix2sky(data, 0))
    assert np.all(w.wcs_world2pix(data, 0) == w.wcs_sky2pix(data, 0))


def test_dict_init():
    """
    Test that WCS can be initialized with a dict-like object
    """

    # Dictionary with no actual WCS, returns identity transform
    w = wcs.WCS({})

    xp, yp = w.wcs_world2pix(41., 2., 1)

    assert_array_almost_equal_nulp(xp[0], 41., 10)
    assert_array_almost_equal_nulp(yp[0], 2., 10)

    # Valid WCS
    w = wcs.WCS({'CTYPE1': 'GLON-CAR',
                 'CTYPE2': 'GLAT-CAR',
                 'CUNIT1': 'deg',
                 'CUNIT2': 'deg',
                 'CRPIX1': 1,
                 'CRPIX2': 1,
                 'CRVAL1': 40.,
                 'CRVAL2': 0.,
                 'CDELT1': -0.1,
                 'CDELT2': 0.1})

    xp, yp = w.wcs_world2pix(41., 2., 0)

    assert_array_almost_equal_nulp(xp[0], -10., 10)
    assert_array_almost_equal_nulp(yp[0], 20., 10)


@raises(TypeError)
def test_extra_kwarg():
    """
    Issue #444
    """
    w = wcs.WCS()
    data = np.random.rand(100, 2)
    w.wcs_pix2sky(data, origin=1)


def test_3d_shapes():
    """
    Issue #444
    """
    w = wcs.WCS(naxis=3)
    data = np.random.rand(100, 3)
    result = w.wcs_pix2sky(data, 1)
    assert result.shape == (100, 3)
    result = w.wcs_pix2sky(
        data[..., 0], data[..., 1], data[..., 2], 1)
    assert len(result) == 3


def test_preserve_shape():

    w = wcs.WCS(naxis=2)

    x = np.random.random((2,3,4))
    y = np.random.random((2,3,4))

    xw, yw = w.wcs_pix2world(x, y, 1)

    assert xw.shape == (2,3,4)
    assert yw.shape == (2,3,4)

    xp, yp = w.wcs_world2pix(x, y, 1)

    assert xp.shape == (2,3,4)
    assert yp.shape == (2,3,4)


def test_shape_mismatch():

    w = wcs.WCS(naxis=2)

    x = np.random.random((2,3,4))
    y = np.random.random((3,2,4))

    with pytest.raises(ValueError) as exc:
        xw, yw = w.wcs_pix2world(x, y, 1)
    assert exc.value.args[0] == "coordinate arrays are not the same shape"

    with pytest.raises(ValueError) as exc:
        xp, yp = w.wcs_world2pix(x, y, 1)
    assert exc.value.args[0] == "coordinate arrays are not the same shape"
