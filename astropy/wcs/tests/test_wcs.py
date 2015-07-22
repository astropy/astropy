# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

from __future__ import absolute_import, division, print_function, unicode_literals

from ...extern import six

import io
import os
import warnings
from datetime import datetime

import numpy as np
from numpy.testing import (
    assert_allclose, assert_array_almost_equal, assert_array_almost_equal_nulp,
    assert_array_equal)

from ...tests.helper import raises, catch_warnings, pytest
from ... import wcs
from .. import _wcs
from ...utils.data import (
    get_pkg_data_filenames, get_pkg_data_contents, get_pkg_data_filename)
from ...utils.misc import NumpyRNGContext
from ...io import fits


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

        all_wcs = wcs.find_all_wcs(header)
        assert len(all_wcs) == 9

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


def test_fixes():
    """
    From github issue #36
    """
    def run():
        header = get_pkg_data_contents(
            'data/nonstandard_units.hdr', encoding='binary')
        try:
            w = wcs.WCS(header, translate_units='dhs')
        except wcs.InvalidTransformError:
            pass
        else:
            assert False, "Expected InvalidTransformError"

    with catch_warnings(wcs.FITSFixedWarning) as w:
        run()

    assert len(w) == 2
    for item in w:
        if 'unitfix' in str(item.message):
            assert 'Hz' in str(item.message)
            assert 'M/S' in str(item.message)
            assert 'm/s' in str(item.message)


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


def test_pix2world():
    """
    From github issue #1463
    """
    # TODO: write this to test the expected output behavior of pix2world,
    # currently this just makes sure it doesn't error out in unexpected ways
    filename = get_pkg_data_filename('data/sip2.fits')
    with catch_warnings(wcs.wcs.FITSFixedWarning) as caught_warnings:
        # this raises a warning unimportant for this testing the pix2world
        #   FITSFixedWarning(u'The WCS transformation has more axes (2) than the
        #        image it is associated with (0)')
        ww = wcs.WCS(filename)

        # might as well monitor for changing behavior
        assert len(caught_warnings) == 1

    n = 3
    pixels = (np.arange(n)*np.ones((2, n))).T
    result = ww.wcs_pix2world(pixels, 0, ra_dec_order=True)

    # Catch #2791
    ww.wcs_pix2world(pixels[..., 0], pixels[..., 1], 0, ra_dec_order=True)

    close_enough = 1e-8
    # assuming that the data of sip2.fits doesn't change
    answer = np.array([[0.00024976, 0.00023018],
                       [0.00023043, -0.00024997]])

    assert np.all(np.abs(ww.wcs.pc-answer) < close_enough)

    answer = np.array([[ 202.39265216,   47.17756518],
                       [ 202.39335826,   47.17754619],
                       [ 202.39406436,   47.1775272 ]])

    assert  np.all(np.abs(result-answer) < close_enough)


def test_load_fits_path():
    fits_name = get_pkg_data_filename('data/sip.fits')
    w = wcs.WCS(fits_name)


def test_dict_init():
    """
    Test that WCS can be initialized with a dict-like object
    """

    # Dictionary with no actual WCS, returns identity transform
    w = wcs.WCS({})

    xp, yp = w.wcs_world2pix(41., 2., 1)

    assert_array_almost_equal_nulp(xp, 41., 10)
    assert_array_almost_equal_nulp(yp, 2., 10)

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

    assert_array_almost_equal_nulp(xp, -10., 10)
    assert_array_almost_equal_nulp(yp, 20., 10)


@raises(TypeError)
def test_extra_kwarg():
    """
    Issue #444
    """
    w = wcs.WCS()
    with NumpyRNGContext(123456789):
        data = np.random.rand(100, 2)
        w.wcs_pix2world(data, origin=1)


def test_3d_shapes():
    """
    Issue #444
    """
    w = wcs.WCS(naxis=3)
    with NumpyRNGContext(123456789):
        data = np.random.rand(100, 3)
        result = w.wcs_pix2world(data, 1)
        assert result.shape == (100, 3)
        result = w.wcs_pix2world(
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


def test_broadcasting():
    w = wcs.WCS(naxis=2)

    x = np.random.random((2,3,4))
    y = 1

    xp, yp = w.wcs_world2pix(x, y, 1)

    assert xp.shape == (2,3,4)
    assert yp.shape == (2,3,4)


def test_shape_mismatch():
    w = wcs.WCS(naxis=2)

    x = np.random.random((2,3,4))
    y = np.random.random((3,2,4))

    with pytest.raises(ValueError) as exc:
        xw, yw = w.wcs_pix2world(x, y, 1)
    assert exc.value.args[0] == "Coordinate arrays are not broadcastable to each other"

    with pytest.raises(ValueError) as exc:
        xp, yp = w.wcs_world2pix(x, y, 1)
    assert exc.value.args[0] == "Coordinate arrays are not broadcastable to each other"

    # There are some ambiguities that need to be worked around when
    # naxis == 1
    w = wcs.WCS(naxis=1)

    x = np.random.random((42, 1))
    xw = w.wcs_pix2world(x, 1)
    assert xw.shape == (42, 1)

    x = np.random.random((42,))
    xw, = w.wcs_pix2world(x, 1)
    assert xw.shape == (42,)


def test_invalid_shape():
    # Issue #1395
    w = wcs.WCS(naxis=2)

    xy = np.random.random((2, 3))
    with pytest.raises(ValueError) as exc:
        xy2 = w.wcs_pix2world(xy, 1)
    assert exc.value.args[0] == 'When providing two arguments, the array must be of shape (N, 2)'

    xy = np.random.random((2, 1))
    with pytest.raises(ValueError) as exc:
        xy2 = w.wcs_pix2world(xy, 1)
    assert exc.value.args[0] == 'When providing two arguments, the array must be of shape (N, 2)'


def test_warning_about_defunct_keywords():
    def run():
        header = get_pkg_data_contents(
            'data/defunct_keywords.hdr', encoding='binary')
        w = wcs.WCS(header)

    with catch_warnings(wcs.FITSFixedWarning) as w:
        run()

    assert len(w) == 4
    for item in w:
        assert 'PCi_ja' in str(item.message)


@raises(wcs.FITSFixedWarning)
def test_warning_about_defunct_keywords_exception():
    def run():
        header = get_pkg_data_contents(
            'data/defunct_keywords.hdr', encoding='binary')
        w = wcs.WCS(header)

    with catch_warnings(wcs.FITSFixedWarning) as w:
        warnings.simplefilter("error", wcs.FITSFixedWarning)
        run()


def test_to_header_string():
    header_string = """
    WCSAXES =                    2 / Number of coordinate axes                      CRPIX1  =                  0.0 / Pixel coordinate of reference point            CRPIX2  =                  0.0 / Pixel coordinate of reference point            CDELT1  =                  1.0 / Coordinate increment at reference point        CDELT2  =                  1.0 / Coordinate increment at reference point        CRVAL1  =                  0.0 / Coordinate value at reference point            CRVAL2  =                  0.0 / Coordinate value at reference point            LATPOLE =                 90.0 / [deg] Native latitude of celestial pole        END"""

    w = wcs.WCS()
    h0 = fits.Header.fromstring(w.to_header_string().strip())
    if 'COMMENT' in h0:
        del h0['COMMENT']
    if '' in h0:
        del h0['']
    h1 = fits.Header.fromstring(header_string.strip())
    assert dict(h0) == dict(h1)


def test_to_fits():
    w = wcs.WCS()
    header_string = w.to_header()
    wfits = w.to_fits()
    assert isinstance(wfits, fits.HDUList)
    assert isinstance(wfits[0], fits.PrimaryHDU)
    assert header_string == wfits[0].header[-8:]


def test_to_header_warning():
    fits_name = get_pkg_data_filename('data/sip.fits')
    x = wcs.WCS(fits_name)
    with catch_warnings() as w:
        x.to_header()
    assert len(w) == 1
    assert 'A_ORDER' in str(w[0])


def test_no_comments_in_header():
    w = wcs.WCS()
    header = w.to_header()
    assert w.wcs.alt not in header
    assert 'COMMENT' + w.wcs.alt.strip() not in header
    assert 'COMMENT' not in header
    wkey = 'P'
    header = w.to_header(key=wkey)
    assert wkey not in header
    assert 'COMMENT' not in header
    assert 'COMMENT' + w.wcs.alt.strip() not in header


@raises(wcs.InvalidTransformError)
def test_find_all_wcs_crash():
    """
    Causes a double free without a recent fix in wcslib_wrap.C
    """
    with open(get_pkg_data_filename("data/too_many_pv.hdr")) as fd:
        header = fd.read()
    # We have to set fix=False here, because one of the fixing tasks is to
    # remove redundant SCAMP distortion parameters when SIP distortion
    # parameters are also present.
    wcses = wcs.find_all_wcs(header, fix=False)


def test_validate():
    with catch_warnings():
        results = wcs.validate(get_pkg_data_filename("data/validate.fits"))
        results_txt = repr(results)
        if wcs._wcs.__version__[0] == '5':
            filename = 'data/validate.5.0.txt'
        else:
            filename = 'data/validate.txt'
        with open(get_pkg_data_filename(filename), "r") as fd:
            lines = fd.readlines()
            assert set([x.strip() for x in lines]) == set([
                x.strip() for x in results_txt.splitlines()])


def test_validate_with_2_wcses():
    # From Issue #2053
    results = wcs.validate(get_pkg_data_filename("data/2wcses.hdr"))

    assert "WCS key 'A':" in six.text_type(results)


def test_all_world2pix(fname=None, ext=0,
                       tolerance=1.0e-4, origin=0,
                       random_npts=25000,
                       adaptive=False, maxiter=20,
                       detect_divergence=True):
    """Test all_world2pix, iterative inverse of all_pix2world"""

    # Open test FITS file:
    if fname is None:
        fname = get_pkg_data_filename('data/j94f05bgq_flt.fits')
        ext = ('SCI',1)
    if not os.path.isfile(fname):
        raise IOError("Input file '{:s}' to 'test_all_world2pix' not found."
                      .format(fname))
    h = fits.open(fname)
    w = wcs.WCS(h[ext].header, h)
    h.close()
    del h

    crpix = w.wcs.crpix
    ncoord = crpix.shape[0]

    # Assume that CRPIX is at the center of the image and that the image has
    # a power-of-2 number of pixels along each axis. Only use the central
    # 1/64 for this testing purpose:
    naxesi_l = list((7./16*crpix).astype(np.int))
    naxesi_u = list((9./16*crpix).astype(np.int))

    # Generate integer indices of pixels (image grid):
    img_pix = np.dstack([i.flatten() for i in
                         np.meshgrid(*map(range, naxesi_l, naxesi_u))])[0]

    # Generage random data (in image coordinates):
    with NumpyRNGContext(123456789):
        rnd_pix = np.random.rand(random_npts, ncoord)

    # Scale random data to cover the central part of the image
    mwidth = 2 * (crpix * 1./8)
    rnd_pix = crpix - 0.5*mwidth + (mwidth-1) * rnd_pix

    # Reference pixel coordinates in image coordinate system (CS):
    test_pix = np.append(img_pix, rnd_pix, axis=0)
    # Reference pixel coordinates in sky CS using forward transformation:
    all_world = w.all_pix2world(test_pix, origin)

    try:
        runtime_begin = datetime.now()
        # Apply the inverse iterative process to pixels in world coordinates
        # to recover the pixel coordinates in image space.
        all_pix = w.all_world2pix(
            all_world, origin, tolerance=tolerance, adaptive=adaptive,
            maxiter=maxiter, detect_divergence=detect_divergence)
        runtime_end = datetime.now()
    except wcs.wcs.NoConvergence as e:
        runtime_end = datetime.now()
        ndiv = 0
        if e.divergent is not None:
            ndiv = e.divergent.shape[0]
            print("There are {} diverging solutions.".format(ndiv))
            print("Indices of diverging solutions:\n{}"
                  .format(e.divergent))
            print("Diverging solutions:\n{}\n"
                  .format(e.best_solution[e.divergent]))
            print("Mean radius of the diverging solutions: {}"
                  .format(np.mean(
                      np.linalg.norm(e.best_solution[e.divergent], axis=1))))
            print("Mean accuracy of the diverging solutions: {}\n"
                  .format(np.mean(
                      np.linalg.norm(e.accuracy[e.divergent], axis=1))))
        else:
            print("There are no diverging solutions.")

        nslow = 0
        if e.slow_conv is not None:
            nslow = e.slow_conv.shape[0]
            print("There are {} slowly converging solutions."
                  .format(nslow))
            print("Indices of slowly converging solutions:\n{}"
                  .format(e.slow_conv))
            print("Slowly converging solutions:\n{}\n"
                  .format(e.best_solution[e.slow_conv]))
        else:
            print("There are no slowly converging solutions.\n")

        print("There are {} converged solutions."
              .format(e.best_solution.shape[0] - ndiv - nslow))
        print("Best solutions (all points):\n{}"
              .format(e.best_solution))
        print("Accuracy:\n{}\n".format(e.accuracy))
        print("\nFinished running 'test_all_world2pix' with errors.\n"
              "ERROR: {}\nRun time: {}\n"
              .format(e.args[0], runtime_end - runtime_begin))
        raise e

    # Compute differences between reference pixel coordinates and
    # pixel coordinates (in image space) recovered from reference
    # pixels in world coordinates:
    errors = np.sqrt(np.sum(np.power(all_pix - test_pix, 2), axis=1))
    meanerr = np.mean(errors)
    maxerr = np.amax(errors)
    print("\nFinished running 'test_all_world2pix'.\n"
          "Mean error = {0:e}  (Max error = {1:e})\n"
          "Run time: {2}\n"
          .format(meanerr, maxerr, runtime_end - runtime_begin))

    assert(maxerr < 2.0*tolerance)


def test_scamp_sip_distortion_parameters():
    """
    Test parsing of WCS parameters with redundant SIP and SCAMP distortion
    parameters.
    """
    header = get_pkg_data_contents('data/validate.fits', encoding='binary')
    w = wcs.WCS(header)
    # Just check that this doesn't raise an exception.
    w.all_pix2world(0, 0, 0)


def test_fixes2():
    """
    From github issue #1854
    """
    header = get_pkg_data_contents(
        'data/nonstandard_units.hdr', encoding='binary')
    with pytest.raises(wcs.InvalidTransformError):
        w = wcs.WCS(header, fix=False)


def test_unit_normalization():
    """
    From github issue #1918
    """
    header = get_pkg_data_contents(
        'data/unit.hdr', encoding='binary')
    w = wcs.WCS(header)
    assert w.wcs.cunit[2] == 'm/s'


def test_footprint_to_file(tmpdir):
    """
    From github issue #1912
    """
    # Arbitrary keywords from real data
    w = wcs.WCS({'CTYPE1': 'RA---ZPN', 'CRUNIT1': 'deg',
                 'CRPIX1': -3.3495999e+02, 'CRVAL1': 3.185790700000e+02,
                 'CTYPE2': 'DEC--ZPN', 'CRUNIT2': 'deg',
                 'CRPIX2': 3.0453999e+03, 'CRVAL2': 4.388538000000e+01,
                 'PV2_1': 1., 'PV2_3': 220.})
    # Just check that this doesn't raise an exception:
    w.footprint_to_file(str(tmpdir.join('test.txt')))


def test_validate_faulty_wcs():
    """
    From github issue #2053
    """
    h = fits.Header()
    # Illegal WCS:
    h['RADESYSA'] = 'ICRS'
    h['PV2_1'] = 1.0
    hdu = fits.PrimaryHDU([[0]], header=h)
    hdulist = fits.HDUList([hdu])
    # Check that this doesn't raise a NameError exception:
    wcs.validate(hdulist)


def test_error_message():
    header = get_pkg_data_contents(
        'data/invalid_header.hdr', encoding='binary')

    with pytest.raises(wcs.InvalidTransformError):
        # Both lines are in here, because 0.4 calls .set within WCS.__init__,
        # whereas 0.3 and earlier did not.
        w = wcs.WCS(header, _do_set=False)
        c = w.all_pix2world([[536.0, 894.0]], 0)


def test_out_of_bounds():
    # See #2107
    header = get_pkg_data_contents('data/zpn-hole.hdr', encoding='binary')
    w = wcs.WCS(header)

    ra, dec = w.wcs_pix2world(110, 110, 0)

    assert np.isnan(ra)
    assert np.isnan(dec)

    ra, dec = w.wcs_pix2world(0, 0, 0)

    assert not np.isnan(ra)
    assert not np.isnan(dec)


def test_calc_footprint_1():
    fits = get_pkg_data_filename('data/sip.fits')
    w = wcs.WCS(fits)

    axes = (1000, 1051)
    ref = np.array([[ 202.39314493,   47.17753352],
                    [ 202.71885939,   46.94630488],
                    [ 202.94631893,   47.15855022],
                    [ 202.72053428,   47.37893142]])
    footprint = w.calc_footprint(axes=axes)
    assert_allclose(footprint, ref)


def test_calc_footprint_2():
    """ Test calc_footprint without distortion. """
    fits = get_pkg_data_filename('data/sip.fits')
    w = wcs.WCS(fits)

    axes = (1000, 1051)
    ref = np.array([[ 202.39265216,   47.17756518],
                    [ 202.7469062 ,   46.91483312],
                    [ 203.11487481,   47.14359319],
                    [ 202.76092671,   47.40745948]])
    footprint = w.calc_footprint(axes=axes, undistort=False)
    assert_allclose(footprint, ref)


def test_calc_footprint_3():
    """ Test calc_footprint with corner of the pixel."""
    w = wcs.WCS()
    w.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]
    w.wcs.crpix = [1.5, 5.5]
    w.wcs.cdelt = [-0.1, 0.1]
    axes = (2, 10)
    ref = np.array([[0.1, -0.5],
                    [0.1, 0.5],
                    [359.9, 0.5],
                    [359.9, -0.5]])

    footprint = w.calc_footprint(axes=axes, undistort=False, center=False)
    assert_allclose(footprint, ref)


def test_sip():
    # See #2107
    header = get_pkg_data_contents('data/irac_sip.hdr', encoding='binary')
    w = wcs.WCS(header)

    x0, y0 = w.sip_pix2foc(200, 200, 0)

    assert_allclose(72, x0, 1e-3)
    assert_allclose(72, y0, 1e-3)

    x1, y1 = w.sip_foc2pix(x0, y0, 0)

    assert_allclose(200, x1, 1e-3)
    assert_allclose(200, y1, 1e-3)


def test_printwcs():
    """
    Just make sure that it runs
    """
    h = get_pkg_data_contents('spectra/orion-freq-1.hdr', encoding='binary')
    w = wcs.WCS(h)
    w.printwcs()
    h = get_pkg_data_contents('data/3d_cd.hdr', encoding='binary')
    w = wcs.WCS(h)
    w.printwcs()


def test_invalid_spherical():
    header = six.text_type("""
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                    8 / array data type
WCSAXES =                    2 / no comment
CTYPE1  = 'RA---TAN' / TAN (gnomic) projection
CTYPE2  = 'DEC--TAN' / TAN (gnomic) projection
EQUINOX =               2000.0 / Equatorial coordinates definition (yr)
LONPOLE =                180.0 / no comment
LATPOLE =                  0.0 / no comment
CRVAL1  =        16.0531567459 / RA  of reference point
CRVAL2  =        23.1148929108 / DEC of reference point
CRPIX1  =                 2129 / X reference pixel
CRPIX2  =                 1417 / Y reference pixel
CUNIT1  = 'deg     ' / X pixel scale units
CUNIT2  = 'deg     ' / Y pixel scale units
CD1_1   =    -0.00912247310646 / Transformation matrix
CD1_2   =    -0.00250608809647 / no comment
CD2_1   =     0.00250608809647 / no comment
CD2_2   =    -0.00912247310646 / no comment
IMAGEW  =                 4256 / Image width,  in pixels.
IMAGEH  =                 2832 / Image height, in pixels.
""")

    f = io.StringIO(header)
    header = fits.Header.fromtextfile(f)

    w = wcs.WCS(header)
    x, y = w.wcs_world2pix(211, -26, 0)
    assert np.isnan(x) and np.isnan(y)


def test_no_iteration():

    # Regression test for #3066

    w = wcs.WCS(naxis=2)

    with pytest.raises(TypeError) as exc:
        iter(w)
    assert exc.value.args[0] == "'WCS' object is not iterable"

    class NewWCS(wcs.WCS):
        pass

    w = NewWCS(naxis=2)

    with pytest.raises(TypeError) as exc:
        iter(w)
    assert exc.value.args[0] == "'NewWCS' object is not iterable"


@pytest.mark.skipif('_wcs.__version__[0] < "5"',
                    reason="TPV only works with wcslib 5.x or later")
def test_sip_tpv_agreement():
    sip_header = get_pkg_data_contents(
        os.path.join("data", "siponly.hdr"), encoding='binary')
    tpv_header = get_pkg_data_contents(
        os.path.join("data", "tpvonly.hdr"), encoding='binary')

    w_sip = wcs.WCS(sip_header)
    w_tpv = wcs.WCS(tpv_header)

    assert_array_almost_equal(
        w_sip.all_pix2world([w_sip.wcs.crpix], 1),
        w_tpv.all_pix2world([w_tpv.wcs.crpix], 1))

    w_sip2 = wcs.WCS(w_sip.to_header())
    w_tpv2 = wcs.WCS(w_tpv.to_header())

    assert_array_almost_equal(
        w_sip.all_pix2world([w_sip.wcs.crpix], 1),
        w_sip2.all_pix2world([w_sip.wcs.crpix], 1))
    assert_array_almost_equal(
        w_tpv.all_pix2world([w_sip.wcs.crpix], 1),
        w_tpv2.all_pix2world([w_sip.wcs.crpix], 1))
    assert_array_almost_equal(
        w_sip2.all_pix2world([w_sip.wcs.crpix], 1),
        w_tpv2.all_pix2world([w_tpv.wcs.crpix], 1))


@pytest.mark.skipif('_wcs.__version__[0] < "5"',
                    reason="TPV only works with wcslib 5.x or later")
def test_tpv_copy():
    # See #3904

    tpv_header = get_pkg_data_contents(
        os.path.join("data", "tpvonly.hdr"), encoding='binary')

    w_tpv = wcs.WCS(tpv_header)

    ra, dec = w_tpv.wcs_pix2world([0, 100, 200], [0, -100, 200], 0)
    assert ra[0] != ra[1] and ra[1] != ra[2]
    assert dec[0] != dec[1] and dec[1] != dec[2]


def test_hst_wcs():
    path = get_pkg_data_filename("data/dist_lookup.fits.gz")

    hdulist = fits.open(path)
    # wcslib will complain about the distortion parameters if they
    # weren't correctly deleted from the header
    wcs.WCS(hdulist[1].header, hdulist)
