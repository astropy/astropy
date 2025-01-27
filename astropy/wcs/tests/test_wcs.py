# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
from contextlib import nullcontext
from datetime import datetime

import numpy as np
import pytest
from numpy.testing import (
    assert_allclose,
    assert_array_almost_equal,
    assert_array_almost_equal_nulp,
    assert_array_equal,
)
from packaging.version import Version

from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.nddata import Cutout2D
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import (
    get_pkg_data_contents,
    get_pkg_data_filename,
    get_pkg_data_filenames,
)
from astropy.utils.exceptions import (
    AstropyDeprecationWarning,
    AstropyUserWarning,
    AstropyWarning,
)
from astropy.utils.misc import NumpyRNGContext
from astropy.wcs import _wcs

_WCSLIB_VER = Version(_wcs.__version__)


# NOTE: User can choose to use system wcslib instead of bundled.
def ctx_for_v71_dateref_warnings():
    if _WCSLIB_VER >= Version("7.1") and _WCSLIB_VER < Version("7.3"):
        ctx = pytest.warns(
            wcs.FITSFixedWarning,
            match=(
                r"'datfix' made the change 'Set DATE-REF to '1858-11-17' from"
                r" MJD-REF'\."
            ),
        )
    else:
        ctx = nullcontext()
    return ctx


class TestMaps:
    def setup_method(self):
        # get the list of the hdr files that we want to test
        self._file_list = list(get_pkg_data_filenames("data/maps", pattern="*.hdr"))

    def test_consistency(self):
        # Check to see that we actually have the list we expect, so that we
        # do not get in a situation where the list is empty or incomplete and
        # the tests still seem to pass correctly.

        # how many do we expect to see?
        n_data_files = 28

        assert len(self._file_list) == n_data_files, (
            f"test_spectra has wrong number data files: found {len(self._file_list)},"
            f" expected  {n_data_files}"
        )

    def test_maps(self):
        for filename in self._file_list:
            # use the base name of the file, so we get more useful messages
            # for failing tests.
            filename = os.path.basename(filename)
            # Now find the associated file in the installed wcs test directory.
            header = get_pkg_data_contents(
                os.path.join("data", "maps", filename), encoding="binary"
            )
            # finally run the test.
            wcsobj = wcs.WCS(header)
            world = wcsobj.wcs_pix2world([[97, 97]], 1)
            assert_array_almost_equal(world, [[285.0, -66.25]], decimal=1)
            pix = wcsobj.wcs_world2pix([[285.0, -66.25]], 1)
            assert_array_almost_equal(pix, [[97, 97]], decimal=0)


class TestSpectra:
    def setup_method(self):
        self._file_list = list(get_pkg_data_filenames("data/spectra", pattern="*.hdr"))

    def test_consistency(self):
        # Check to see that we actually have the list we expect, so that we
        # do not get in a situation where the list is empty or incomplete and
        # the tests still seem to pass correctly.

        # how many do we expect to see?
        n_data_files = 6

        assert len(self._file_list) == n_data_files, (
            f"test_spectra has wrong number data files: found {len(self._file_list)},"
            f" expected  {n_data_files}"
        )

    def test_spectra(self):
        for filename in self._file_list:
            # use the base name of the file, so we get more useful messages
            # for failing tests.
            filename = os.path.basename(filename)
            # Now find the associated file in the installed wcs test directory.
            header = get_pkg_data_contents(
                os.path.join("data", "spectra", filename), encoding="binary"
            )
            # finally run the test.
            if _WCSLIB_VER >= Version("7.4"):
                ctx = pytest.warns(
                    wcs.FITSFixedWarning,
                    match=(
                        r"'datfix' made the change 'Set MJD-OBS to 53925\.853472 from"
                        r" DATE-OBS'\."
                    ),
                )
            else:
                ctx = nullcontext()
            with ctx:
                all_wcs = wcs.find_all_wcs(header)

            assert len(all_wcs) == 9


def test_fixes():
    """
    From github issue #36
    """
    header = get_pkg_data_contents("data/nonstandard_units.hdr", encoding="binary")

    with (
        pytest.raises(wcs.InvalidTransformError),
        pytest.warns(wcs.FITSFixedWarning) as w,
    ):
        wcs.WCS(header, translate_units="dhs")

    if Version("7.4") <= _WCSLIB_VER < Version("7.6"):
        assert len(w) == 3
        assert "'datfix' made the change 'Success'." in str(w.pop().message)
    else:
        assert len(w) == 2

    first_wmsg = str(w[0].message)
    assert "unitfix" in first_wmsg and "Hz" in first_wmsg and "M/S" in first_wmsg
    assert "plane angle" in str(w[1].message) and "m/s" in str(w[1].message)


# Ignore "PV2_2 = 0.209028857410973 invalid keyvalue" warning seen on Windows.
@pytest.mark.filterwarnings(r"ignore:PV2_2")
def test_outside_sky():
    """
    From github issue #107
    """
    header = get_pkg_data_contents("data/outside_sky.hdr", encoding="binary")
    w = wcs.WCS(header)

    assert np.all(np.isnan(w.wcs_pix2world([[100.0, 500.0]], 0)))  # outside sky
    assert np.all(np.isnan(w.wcs_pix2world([[200.0, 200.0]], 0)))  # outside sky
    assert not np.any(np.isnan(w.wcs_pix2world([[1000.0, 1000.0]], 0)))


def test_pix2world():
    """
    From github issue #1463
    """
    # TODO: write this to test the expected output behavior of pix2world,
    # currently this just makes sure it doesn't error out in unexpected ways
    # (and compares `wcs.pc` and `result` values?)
    filename = get_pkg_data_filename("data/sip2.fits")
    with pytest.warns(wcs.FITSFixedWarning) as caught_warnings:
        # this raises a warning unimportant for this testing the pix2world
        #   FITSFixedWarning(u'The WCS transformation has more axes (2) than
        #        the image it is associated with (0)')
        ww = wcs.WCS(filename)

    # might as well monitor for changing behavior
    if Version("7.4") <= _WCSLIB_VER < Version("7.6"):
        assert len(caught_warnings) == 2
    else:
        assert len(caught_warnings) == 1

    n = 3
    pixels = (np.arange(n) * np.ones((2, n))).T
    result = ww.wcs_pix2world(pixels, 0, ra_dec_order=True)

    # Catch #2791
    ww.wcs_pix2world(pixels[..., 0], pixels[..., 1], 0, ra_dec_order=True)

    # assuming that the data of sip2.fits doesn't change
    answer = np.array([[0.00024976, 0.00023018], [0.00023043, -0.00024997]])

    assert np.allclose(ww.wcs.pc, answer, atol=1.0e-8)

    answer = np.array(
        [
            [202.39265216, 47.17756518],
            [202.39335826, 47.17754619],
            [202.39406436, 47.1775272],
        ]
    )

    assert np.allclose(result, answer, atol=1.0e-8, rtol=1.0e-10)


def test_load_fits_path():
    fits_name = get_pkg_data_filename("data/sip.fits")
    with pytest.warns(wcs.FITSFixedWarning):
        wcs.WCS(fits_name)


def test_dict_init():
    """
    Test that WCS can be initialized with a dict-like object
    """

    # Dictionary with no actual WCS, returns identity transform
    with ctx_for_v71_dateref_warnings():
        w = wcs.WCS({})

    xp, yp = w.wcs_world2pix(41.0, 2.0, 1)

    assert_array_almost_equal_nulp(xp, 41.0, 10)
    assert_array_almost_equal_nulp(yp, 2.0, 10)

    # Valid WCS
    hdr = {
        "CTYPE1": "GLON-CAR",
        "CTYPE2": "GLAT-CAR",
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CRPIX1": 1,
        "CRPIX2": 1,
        "CRVAL1": 40.0,
        "CRVAL2": 0.0,
        "CDELT1": -0.1,
        "CDELT2": 0.1,
    }
    if _WCSLIB_VER >= Version("7.1"):
        hdr["DATEREF"] = "1858-11-17"

    if _WCSLIB_VER >= Version("7.4"):
        ctx = pytest.warns(
            wcs.wcs.FITSFixedWarning,
            match=r"'datfix' made the change 'Set MJDREF to 0\.000000 from DATEREF'\.",
        )
    else:
        ctx = nullcontext()

    with ctx:
        w = wcs.WCS(hdr)

    xp, yp = w.wcs_world2pix(41.0, 2.0, 0)

    assert_array_almost_equal_nulp(xp, -10.0, 10)
    assert_array_almost_equal_nulp(yp, 20.0, 10)


def test_extra_kwarg():
    """
    Issue #444
    """
    w = wcs.WCS()
    with NumpyRNGContext(123456789):
        data = np.random.rand(100, 2)
        with pytest.raises(TypeError):
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
        result = w.wcs_pix2world(data[..., 0], data[..., 1], data[..., 2], 1)
        assert len(result) == 3


def test_preserve_shape():
    w = wcs.WCS(naxis=2)

    x = np.random.random((2, 3, 4))
    y = np.random.random((2, 3, 4))

    xw, yw = w.wcs_pix2world(x, y, 1)

    assert xw.shape == (2, 3, 4)
    assert yw.shape == (2, 3, 4)

    xp, yp = w.wcs_world2pix(x, y, 1)

    assert xp.shape == (2, 3, 4)
    assert yp.shape == (2, 3, 4)


def test_broadcasting():
    w = wcs.WCS(naxis=2)

    x = np.random.random((2, 3, 4))
    y = 1

    xp, yp = w.wcs_world2pix(x, y, 1)

    assert xp.shape == (2, 3, 4)
    assert yp.shape == (2, 3, 4)


def test_shape_mismatch():
    w = wcs.WCS(naxis=2)

    x = np.random.random((2, 3, 4))
    y = np.random.random((3, 2, 4))

    MESSAGE = r"Coordinate arrays are not broadcastable to each other"
    with pytest.raises(ValueError, match=MESSAGE):
        xw, yw = w.wcs_pix2world(x, y, 1)

    with pytest.raises(ValueError, match=MESSAGE):
        xp, yp = w.wcs_world2pix(x, y, 1)

    # There are some ambiguities that need to be worked around when
    # naxis == 1
    w = wcs.WCS(naxis=1)

    x = np.random.random((42, 1))
    xw = w.wcs_pix2world(x, 1)
    assert xw.shape == (42, 1)

    x = np.random.random((42,))
    (xw,) = w.wcs_pix2world(x, 1)
    assert xw.shape == (42,)


def test_invalid_shape():
    """Issue #1395"""
    MESSAGE = r"When providing two arguments, the array must be of shape [(]N, 2[)]"

    w = wcs.WCS(naxis=2)

    xy = np.random.random((2, 3))
    with pytest.raises(ValueError, match=MESSAGE):
        w.wcs_pix2world(xy, 1)

    xy = np.random.random((2, 1))
    with pytest.raises(ValueError, match=MESSAGE):
        w.wcs_pix2world(xy, 1)


def test_warning_about_defunct_keywords():
    header = get_pkg_data_contents("data/defunct_keywords.hdr", encoding="binary")
    if Version("7.4") <= _WCSLIB_VER < Version("7.6"):
        n_warn = 5
    else:
        n_warn = 4

    # Make sure the warnings come out every time...
    for _ in range(2):
        with pytest.warns(wcs.FITSFixedWarning) as w:
            wcs.WCS(header)

        assert len(w) == n_warn
        # 7.4 adds a fifth warning "'datfix' made the change 'Success'."
        for item in w[:4]:
            assert "PCi_ja" in str(item.message)


def test_warning_about_defunct_keywords_exception():
    header = get_pkg_data_contents("data/defunct_keywords.hdr", encoding="binary")
    with pytest.warns(wcs.FITSFixedWarning):
        wcs.WCS(header)


def test_to_header_string():
    # fmt: off
    hdrstr = (
        "WCSAXES =                    2 / Number of coordinate axes                      ",
        "CRPIX1  =                  0.0 / Pixel coordinate of reference point            ",
        "CRPIX2  =                  0.0 / Pixel coordinate of reference point            ",
        "CDELT1  =                  1.0 / Coordinate increment at reference point        ",
        "CDELT2  =                  1.0 / Coordinate increment at reference point        ",
        "CRVAL1  =                  0.0 / Coordinate value at reference point            ",
        "CRVAL2  =                  0.0 / Coordinate value at reference point            ",
        "LATPOLE =                 90.0 / [deg] Native latitude of celestial pole        ",
    )
    # fmt: on

    if _WCSLIB_VER >= Version("7.3"):
        # fmt: off
        hdrstr += (
            "MJDREF  =                  0.0 / [d] MJD of fiducial time                       ",
        )
        # fmt: on

    elif _WCSLIB_VER >= Version("7.1"):
        # fmt: off
        hdrstr += (
            "DATEREF = '1858-11-17'         / ISO-8601 fiducial time                         ",
            "MJDREFI =                  0.0 / [d] MJD of fiducial time, integer part         ",
            "MJDREFF =                  0.0 / [d] MJD of fiducial time, fractional part      "
        )
        # fmt: on

    hdrstr += ("END",)

    header_string = "".join(hdrstr)

    w = wcs.WCS()
    h0 = fits.Header.fromstring(w.to_header_string().strip())
    if "COMMENT" in h0:
        del h0["COMMENT"]
    if "" in h0:
        del h0[""]
    h1 = fits.Header.fromstring(header_string.strip())
    assert dict(h0) == dict(h1)


def test_to_fits():
    nrec = 11 if _WCSLIB_VER >= Version("7.1") else 8
    if _WCSLIB_VER < Version("7.1"):
        nrec = 8
    elif _WCSLIB_VER < Version("7.3"):
        nrec = 11
    else:
        nrec = 9

    w = wcs.WCS()
    header_string = w.to_header()
    wfits = w.to_fits()
    assert isinstance(wfits, fits.HDUList)
    assert isinstance(wfits[0], fits.PrimaryHDU)
    assert header_string == wfits[0].header[-nrec:]


def test_to_header_warning():
    fits_name = get_pkg_data_filename("data/sip.fits")
    with pytest.warns(wcs.FITSFixedWarning):
        x = wcs.WCS(fits_name)
    with pytest.warns(AstropyWarning, match="A_ORDER") as w:
        x.to_header()
    assert len(w) == 1


def test_no_comments_in_header():
    w = wcs.WCS()
    header = w.to_header()
    assert w.wcs.alt not in header
    assert "COMMENT" + w.wcs.alt.strip() not in header
    assert "COMMENT" not in header
    wkey = "P"
    header = w.to_header(key=wkey)
    assert wkey not in header
    assert "COMMENT" not in header
    assert "COMMENT" + w.wcs.alt.strip() not in header


def test_find_all_wcs_crash():
    """
    Causes a double free without a recent fix in wcslib_wrap.C
    """
    with open(get_pkg_data_filename("data/too_many_pv.hdr")) as fd:
        header = fd.read()
    # We have to set fix=False here, because one of the fixing tasks is to
    # remove redundant SCAMP distortion parameters when SIP distortion
    # parameters are also present.
    with pytest.raises(wcs.InvalidTransformError), pytest.warns(wcs.FITSFixedWarning):
        wcs.find_all_wcs(header, fix=False)


# NOTE: Warning bubbles up from C layer during wcs.validate() and
# is hard to catch, so we just ignore it.
@pytest.mark.filterwarnings("ignore")
def test_validate():
    results = wcs.validate(get_pkg_data_filename("data/validate.fits"))
    results_txt = sorted({x.strip() for x in repr(results).splitlines()})
    if _WCSLIB_VER >= Version("7.6"):
        filename = "data/validate.7.6.txt"
    elif _WCSLIB_VER >= Version("7.4"):
        filename = "data/validate.7.4.txt"
    elif _WCSLIB_VER >= Version("6.0"):
        filename = "data/validate.6.txt"
    elif _WCSLIB_VER >= Version("5.13"):
        filename = "data/validate.5.13.txt"
    elif _WCSLIB_VER >= Version("5.0"):
        filename = "data/validate.5.0.txt"
    else:
        filename = "data/validate.txt"
    with open(get_pkg_data_filename(filename)) as fd:
        lines = fd.readlines()
    assert sorted({x.strip() for x in lines}) == results_txt


@pytest.mark.filterwarnings("ignore")
def test_validate_wcs_tab():
    results = wcs.validate(get_pkg_data_filename("data/tab-time-last-axis.fits"))
    results_txt = sorted({x.strip() for x in repr(results).splitlines()})
    assert results_txt == [
        "",
        "HDU 0 (PRIMARY):",
        "HDU 1 (WCS-TABLE):",
        "No issues.",
        "WCS key ' ':",
    ]


def test_validate_with_2_wcses():
    # From Issue #2053
    with pytest.warns(AstropyUserWarning):
        results = wcs.validate(get_pkg_data_filename("data/2wcses.hdr"))

    assert "WCS key 'A':" in str(results)


def test_crpix_maps_to_crval():
    twcs = wcs.WCS(naxis=2)
    twcs.wcs.crval = [251.29, 57.58]
    twcs.wcs.cdelt = [1, 1]
    twcs.wcs.crpix = [507, 507]
    twcs.wcs.pc = np.array([[7.7e-6, 3.3e-5], [3.7e-5, -6.8e-6]])
    twcs._naxis = [1014, 1014]
    twcs.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
    a = np.array(
        [
            [0, 0, 5.33092692e-08, 3.73753773e-11, -2.02111473e-13],
            [0, 2.44084308e-05, 2.81394789e-11, 5.17856895e-13, 0.0],
            [-2.41334657e-07, 1.29289255e-10, 2.35753629e-14, 0.0, 0.0],
            [-2.37162007e-10, 5.43714947e-13, 0.0, 0.0, 0.0],
            [-2.81029767e-13, 0.0, 0.0, 0.0, 0.0],
        ]
    )
    b = np.array(
        [
            [0, 0, 2.99270374e-05, -2.38136074e-10, 7.23205168e-13],
            [0, -1.71073858e-07, 6.31243431e-11, -5.16744347e-14, 0.0],
            [6.95458963e-06, -3.08278961e-10, -1.75800917e-13, 0.0, 0.0],
            [3.51974159e-11, 5.60993016e-14, 0.0, 0.0, 0.0],
            [-5.92438525e-13, 0.0, 0.0, 0.0, 0.0],
        ]
    )
    twcs.sip = wcs.Sip(a, b, None, None, twcs.wcs.crpix)
    twcs.wcs.set()
    pscale = np.sqrt(wcs.utils.proj_plane_pixel_area(twcs))

    # test that CRPIX maps to CRVAL:
    assert_allclose(
        twcs.wcs_pix2world(*twcs.wcs.crpix, 1),
        twcs.wcs.crval,
        rtol=0.0,
        atol=1e-6 * pscale,
    )

    # test that CRPIX maps to CRVAL:
    assert_allclose(
        twcs.all_pix2world(*twcs.wcs.crpix, 1),
        twcs.wcs.crval,
        rtol=0.0,
        atol=1e-6 * pscale,
    )


def test_all_world2pix(
    fname=None,
    ext=0,
    tolerance=1.0e-4,
    origin=0,
    random_npts=25000,
    adaptive=False,
    maxiter=20,
    detect_divergence=True,
):
    """Test all_world2pix, iterative inverse of all_pix2world"""

    # Open test FITS file:
    if fname is None:
        fname = get_pkg_data_filename("data/j94f05bgq_flt.fits")
        ext = ("SCI", 1)
    if not os.path.isfile(fname):
        raise OSError(f"Input file '{fname:s}' to 'test_all_world2pix' not found.")
    h = fits.open(fname)
    w = wcs.WCS(h[ext].header, h)
    h.close()
    del h

    crpix = w.wcs.crpix
    ncoord = crpix.shape[0]

    # Assume that CRPIX is at the center of the image and that the image has
    # a power-of-2 number of pixels along each axis. Only use the central
    # 1/64 for this testing purpose:
    naxesi_l = list((7.0 / 16 * crpix).astype(int))
    naxesi_u = list((9.0 / 16 * crpix).astype(int))

    # Generate integer indices of pixels (image grid):
    img_pix = np.dstack(
        [i.flatten() for i in np.meshgrid(*map(range, naxesi_l, naxesi_u))]
    )[0]

    # Generate random data (in image coordinates):
    with NumpyRNGContext(123456789):
        rnd_pix = np.random.rand(random_npts, ncoord)

    # Scale random data to cover the central part of the image
    mwidth = 2 * (crpix * 1.0 / 8)
    rnd_pix = crpix - 0.5 * mwidth + (mwidth - 1) * rnd_pix

    # Reference pixel coordinates in image coordinate system (CS):
    test_pix = np.append(img_pix, rnd_pix, axis=0)
    # Reference pixel coordinates in sky CS using forward transformation:
    all_world = w.all_pix2world(test_pix, origin)

    try:
        runtime_begin = datetime.now()
        # Apply the inverse iterative process to pixels in world coordinates
        # to recover the pixel coordinates in image space.
        all_pix = w.all_world2pix(
            all_world,
            origin,
            tolerance=tolerance,
            adaptive=adaptive,
            maxiter=maxiter,
            detect_divergence=detect_divergence,
        )
        runtime_end = datetime.now()
    except wcs.wcs.NoConvergence as e:
        runtime_end = datetime.now()
        ndiv = 0
        if e.divergent is not None:
            ndiv = e.divergent.shape[0]
            print(f"There are {ndiv} diverging solutions.")
            print(f"Indices of diverging solutions:\n{e.divergent}")
            print(f"Diverging solutions:\n{e.best_solution[e.divergent]}\n")
            print(
                "Mean radius of the diverging solutions:"
                f" {np.mean(np.linalg.norm(e.best_solution[e.divergent], axis=1))}"
            )
            print(
                "Mean accuracy of the diverging solutions:"
                f" {np.mean(np.linalg.norm(e.accuracy[e.divergent], axis=1))}\n"
            )
        else:
            print("There are no diverging solutions.")

        nslow = 0
        if e.slow_conv is not None:
            nslow = e.slow_conv.shape[0]
            print(f"There are {nslow} slowly converging solutions.")
            print(f"Indices of slowly converging solutions:\n{e.slow_conv}")
            print(f"Slowly converging solutions:\n{e.best_solution[e.slow_conv]}\n")
        else:
            print("There are no slowly converging solutions.\n")

        print(
            f"There are {e.best_solution.shape[0] - ndiv - nslow} converged solutions."
        )
        print(f"Best solutions (all points):\n{e.best_solution}")
        print(f"Accuracy:\n{e.accuracy}\n")
        print(
            "\nFinished running 'test_all_world2pix' with errors.\n"
            f"ERROR: {e.args[0]}\nRun time: {runtime_end - runtime_begin}\n"
        )
        raise e

    # Compute differences between reference pixel coordinates and
    # pixel coordinates (in image space) recovered from reference
    # pixels in world coordinates:
    errors = np.sqrt(np.sum(np.power(all_pix - test_pix, 2), axis=1))
    meanerr = np.mean(errors)
    maxerr = np.amax(errors)
    print(
        "\nFinished running 'test_all_world2pix'.\n"
        f"Mean error = {meanerr:e}  (Max error = {maxerr:e})\n"
        f"Run time: {runtime_end - runtime_begin}\n"
    )

    assert maxerr < 2.0 * tolerance


def test_scamp_sip_distortion_parameters():
    """
    Test parsing of WCS parameters with redundant SIP and SCAMP distortion
    parameters.
    """
    header = get_pkg_data_contents("data/validate.fits", encoding="binary")
    with pytest.warns(wcs.FITSFixedWarning):
        w = wcs.WCS(header)
    # Just check that this doesn't raise an exception.
    w.all_pix2world(0, 0, 0)


def test_fixes2():
    """
    From github issue #1854
    """
    header = get_pkg_data_contents("data/nonstandard_units.hdr", encoding="binary")
    with pytest.raises(wcs.InvalidTransformError):
        wcs.WCS(header, fix=False)


def test_unit_normalization():
    """
    From github issue #1918
    """
    header = get_pkg_data_contents("data/unit.hdr", encoding="binary")
    w = wcs.WCS(header)
    assert w.wcs.cunit[2] == "m/s"


def test_footprint_to_file(tmp_path):
    """
    From github issue #1912
    """
    # Arbitrary keywords from real data
    hdr = {
        "CTYPE1": "RA---ZPN",
        "CRUNIT1": "deg",
        "CRPIX1": -3.3495999e02,
        "CRVAL1": 3.185790700000e02,
        "CTYPE2": "DEC--ZPN",
        "CRUNIT2": "deg",
        "CRPIX2": 3.0453999e03,
        "CRVAL2": 4.388538000000e01,
        "PV2_1": 1.0,
        "PV2_3": 220.0,
        "NAXIS1": 2048,
        "NAXIS2": 1024,
    }
    w = wcs.WCS(hdr)

    testfile = tmp_path / "test.txt"
    w.footprint_to_file(testfile)

    with open(testfile) as f:
        lines = f.readlines()

    assert len(lines) == 4
    assert lines[2] == "ICRS\n"
    assert "color=green" in lines[3]

    w.footprint_to_file(testfile, coordsys="FK5", color="red")

    with open(testfile) as f:
        lines = f.readlines()

    assert len(lines) == 4
    assert lines[2] == "FK5\n"
    assert "color=red" in lines[3]

    with pytest.raises(ValueError):
        w.footprint_to_file(testfile, coordsys="FOO")

    del hdr["NAXIS1"]
    del hdr["NAXIS2"]
    w = wcs.WCS(hdr)
    with pytest.warns(AstropyUserWarning):
        w.footprint_to_file(testfile)


# Ignore FITSFixedWarning about keyrecords following the END keyrecord were
# ignored, which comes from src/astropy_wcs.c . Only a blind catch like this
# seems to work when pytest warnings are turned into exceptions.
@pytest.mark.filterwarnings("ignore")
def test_validate_faulty_wcs():
    """
    From github issue #2053
    """
    h = fits.Header()
    # Illegal WCS:
    h["RADESYSA"] = "ICRS"
    h["PV2_1"] = 1.0
    hdu = fits.PrimaryHDU([[0]], header=h)
    hdulist = fits.HDUList([hdu])
    # Check that this doesn't raise a NameError exception
    wcs.validate(hdulist)


def test_error_message():
    header = get_pkg_data_contents("data/invalid_header.hdr", encoding="binary")

    # make WCS transformation invalid
    hdr = fits.Header.fromstring(header)
    del hdr["PV?_*"]
    hdr["PV1_1"] = 110
    hdr["PV1_2"] = 110
    hdr["PV2_1"] = -110
    hdr["PV2_2"] = -110
    with pytest.raises(wcs.InvalidTransformError):
        with pytest.warns(wcs.FITSFixedWarning):
            w = wcs.WCS(hdr, _do_set=False)
            w.all_pix2world([[536.0, 894.0]], 0)


def test_out_of_bounds():
    # See #2107
    header = get_pkg_data_contents("data/zpn-hole.hdr", encoding="binary")
    w = wcs.WCS(header)

    ra, dec = w.wcs_pix2world(110, 110, 0)

    assert np.isnan(ra)
    assert np.isnan(dec)

    ra, dec = w.wcs_pix2world(0, 0, 0)

    assert not np.isnan(ra)
    assert not np.isnan(dec)


def test_calc_footprint_1():
    fits = get_pkg_data_filename("data/sip.fits")
    with pytest.warns(wcs.FITSFixedWarning):
        w = wcs.WCS(fits)

        axes = (1000, 1051)
        ref = np.array(
            [
                [202.39314493, 47.17753352],
                [202.71885939, 46.94630488],
                [202.94631893, 47.15855022],
                [202.72053428, 47.37893142],
            ]
        )
        footprint = w.calc_footprint(axes=axes)
        assert_allclose(footprint, ref)


def test_calc_footprint_2():
    """Test calc_footprint without distortion."""
    fits = get_pkg_data_filename("data/sip.fits")
    with pytest.warns(wcs.FITSFixedWarning):
        w = wcs.WCS(fits)

        axes = (1000, 1051)
        ref = np.array(
            [
                [202.39265216, 47.17756518],
                [202.7469062, 46.91483312],
                [203.11487481, 47.14359319],
                [202.76092671, 47.40745948],
            ]
        )
        footprint = w.calc_footprint(axes=axes, undistort=False)
        assert_allclose(footprint, ref)


def test_calc_footprint_3():
    """Test calc_footprint with corner of the pixel."""
    w = wcs.WCS()
    w.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]
    w.wcs.crpix = [1.5, 5.5]
    w.wcs.cdelt = [-0.1, 0.1]
    axes = (2, 10)
    ref = np.array([[0.1, -0.5], [0.1, 0.5], [359.9, 0.5], [359.9, -0.5]])

    footprint = w.calc_footprint(axes=axes, undistort=False, center=False)
    assert_allclose(footprint, ref)


def test_sip():
    # See #2107
    header = get_pkg_data_contents("data/irac_sip.hdr", encoding="binary")
    w = wcs.WCS(header)

    x0, y0 = w.sip_pix2foc(200, 200, 0)

    assert_allclose(72, x0, 1e-3)
    assert_allclose(72, y0, 1e-3)

    x1, y1 = w.sip_foc2pix(x0, y0, 0)

    assert_allclose(200, x1, 1e-3)
    assert_allclose(200, y1, 1e-3)


def test_sub_3d_with_sip():
    # See #10527
    header = get_pkg_data_contents("data/irac_sip.hdr", encoding="binary")
    header = fits.Header.fromstring(header)
    header["NAXIS"] = 3
    header.set("NAXIS3", 64, after=header.index("NAXIS2"))
    w = wcs.WCS(header, naxis=2)
    assert w.naxis == 2


def test_printwcs(capsys):
    """
    Just make sure that it runs
    """
    h = get_pkg_data_contents("data/spectra/orion-freq-1.hdr", encoding="binary")
    with pytest.warns(wcs.FITSFixedWarning):
        w = wcs.WCS(h)
        w.printwcs()
        captured = capsys.readouterr()
        assert "WCS Keywords" in captured.out
    h = get_pkg_data_contents("data/3d_cd.hdr", encoding="binary")
    w = wcs.WCS(h)
    w.printwcs()
    captured = capsys.readouterr()
    assert "WCS Keywords" in captured.out


def test_invalid_spherical():
    header = """
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
    """

    f = io.StringIO(header)
    header = fits.Header.fromtextfile(f)

    w = wcs.WCS(header)
    x, y = w.wcs_world2pix(211, -26, 0)
    assert np.isnan(x) and np.isnan(y)


def test_no_iteration():
    """Regression test for #3066"""
    MESSAGE = "'{}' object is not iterable"

    w = wcs.WCS(naxis=2)

    with pytest.raises(TypeError, match=MESSAGE.format("WCS")):
        iter(w)

    class NewWCS(wcs.WCS):
        pass

    w = NewWCS(naxis=2)

    with pytest.raises(TypeError, match=MESSAGE.format("NewWCS")):
        iter(w)


@pytest.mark.skipif(
    _wcs.__version__[0] < "5", reason="TPV only works with wcslib 5.x or later"
)
def test_sip_tpv_agreement():
    sip_header = get_pkg_data_contents(
        os.path.join("data", "siponly.hdr"), encoding="binary"
    )
    tpv_header = get_pkg_data_contents(
        os.path.join("data", "tpvonly.hdr"), encoding="binary"
    )

    with (
        pytest.warns(wcs.FITSFixedWarning),
        pytest.warns(
            AstropyWarning, match="Some non-standard WCS keywords were excluded"
        ),
    ):
        w_sip = wcs.WCS(sip_header)
        w_tpv = wcs.WCS(tpv_header)

        assert_array_almost_equal(
            w_sip.all_pix2world([w_sip.wcs.crpix], 1),
            w_tpv.all_pix2world([w_tpv.wcs.crpix], 1),
        )

        w_sip2 = wcs.WCS(w_sip.to_header())
        w_tpv2 = wcs.WCS(w_tpv.to_header())

        assert_array_almost_equal(
            w_sip.all_pix2world([w_sip.wcs.crpix], 1),
            w_sip2.all_pix2world([w_sip.wcs.crpix], 1),
        )
        assert_array_almost_equal(
            w_tpv.all_pix2world([w_sip.wcs.crpix], 1),
            w_tpv2.all_pix2world([w_sip.wcs.crpix], 1),
        )
        assert_array_almost_equal(
            w_sip2.all_pix2world([w_sip.wcs.crpix], 1),
            w_tpv2.all_pix2world([w_tpv.wcs.crpix], 1),
        )


def test_tpv_ctype_sip():
    sip_header = fits.Header.fromstring(
        get_pkg_data_contents(os.path.join("data", "siponly.hdr"), encoding="binary")
    )
    tpv_header = fits.Header.fromstring(
        get_pkg_data_contents(os.path.join("data", "tpvonly.hdr"), encoding="binary")
    )
    sip_header.update(tpv_header)
    sip_header["CTYPE1"] = "RA---TAN-SIP"
    sip_header["CTYPE2"] = "DEC--TAN-SIP"

    with (
        pytest.warns(
            wcs.FITSFixedWarning,
            match="Removed redundant SCAMP distortion parameters "
            "because SIP parameters are also present",
        ),
        pytest.warns(
            wcs.FITSFixedWarning,
            match=".*RADECSYS keyword is deprecated, use RADESYSa",
        ),
        pytest.warns(wcs.FITSFixedWarning, match=".*Set MJD-OBS to .* from DATE-OBS"),
    ):
        w_sip = wcs.WCS(sip_header)

    assert w_sip.sip is not None


def test_tpv_ctype_tpv():
    sip_header = fits.Header.fromstring(
        get_pkg_data_contents(os.path.join("data", "siponly.hdr"), encoding="binary")
    )
    tpv_header = fits.Header.fromstring(
        get_pkg_data_contents(os.path.join("data", "tpvonly.hdr"), encoding="binary")
    )
    sip_header.update(tpv_header)
    sip_header["CTYPE1"] = "RA---TPV"
    sip_header["CTYPE2"] = "DEC--TPV"

    with (
        pytest.warns(
            wcs.FITSFixedWarning,
            match="Removed redundant SIP distortion parameters "
            "because CTYPE explicitly specifies TPV distortions",
        ),
        pytest.warns(
            wcs.FITSFixedWarning,
            match=".*RADECSYS keyword is deprecated, use RADESYSa",
        ),
        pytest.warns(wcs.FITSFixedWarning, match=".*Set MJD-OBS to .* from DATE-OBS"),
    ):
        w_sip = wcs.WCS(sip_header)

    assert w_sip.sip is None


def test_tpv_ctype_tan():
    sip_header = fits.Header.fromstring(
        get_pkg_data_contents(os.path.join("data", "siponly.hdr"), encoding="binary")
    )
    tpv_header = fits.Header.fromstring(
        get_pkg_data_contents(os.path.join("data", "tpvonly.hdr"), encoding="binary")
    )
    sip_header.update(tpv_header)
    sip_header["CTYPE1"] = "RA---TAN"
    sip_header["CTYPE2"] = "DEC--TAN"

    with (
        pytest.warns(
            wcs.FITSFixedWarning,
            match="Removed redundant SIP distortion parameters "
            "because SCAMP' PV distortions are also present",
        ),
        pytest.warns(
            wcs.FITSFixedWarning,
            match=".*RADECSYS keyword is deprecated, use RADESYSa",
        ),
        pytest.warns(wcs.FITSFixedWarning, match=".*Set MJD-OBS to .* from DATE-OBS"),
    ):
        w_sip = wcs.WCS(sip_header)

    assert w_sip.sip is None


def test_car_sip_with_pv():
    # https://github.com/astropy/astropy/issues/14255
    header_dict = {
        "SIMPLE": True,
        "BITPIX": -32,
        "NAXIS": 2,
        "NAXIS1": 1024,
        "NAXIS2": 1024,
        "CRPIX1": 512.0,
        "CRPIX2": 512.0,
        "CDELT1": 0.01,
        "CDELT2": 0.01,
        "CRVAL1": 120.0,
        "CRVAL2": 29.0,
        "CTYPE1": "RA---CAR-SIP",
        "CTYPE2": "DEC--CAR-SIP",
        "PV1_1": 120.0,
        "PV1_2": 29.0,
        "PV1_0": 1.0,
        "A_ORDER": 2,
        "A_2_0": 5.0e-4,
        "B_ORDER": 2,
        "B_2_0": 5.0e-4,
    }

    w = wcs.WCS(header_dict)

    assert w.sip is not None

    assert w.wcs.get_pv() == [(1, 1, 120.0), (1, 2, 29.0), (1, 0, 1.0)]

    assert np.allclose(
        w.all_pix2world(header_dict["CRPIX1"], header_dict["CRPIX2"], 1),
        [header_dict["CRVAL1"], header_dict["CRVAL2"]],
    )


@pytest.mark.skipif(
    _wcs.__version__[0] < "5", reason="TPV only works with wcslib 5.x or later"
)
def test_tpv_copy():
    # See #3904

    tpv_header = get_pkg_data_contents(
        os.path.join("data", "tpvonly.hdr"), encoding="binary"
    )

    with pytest.warns(wcs.FITSFixedWarning):
        w_tpv = wcs.WCS(tpv_header)

        ra, dec = w_tpv.wcs_pix2world([0, 100, 200], [0, -100, 200], 0)
        assert ra[0] != ra[1] and ra[1] != ra[2]
        assert dec[0] != dec[1] and dec[1] != dec[2]


def test_hst_wcs():
    path = get_pkg_data_filename("data/dist_lookup.fits.gz")

    with fits.open(path) as hdulist:
        # wcslib will complain about the distortion parameters if they
        # weren't correctly deleted from the header
        w = wcs.WCS(hdulist[1].header, hdulist)

        # Check pixel scale and area
        assert_quantity_allclose(
            w.proj_plane_pixel_scales(), [1.38484378e-05, 1.39758488e-05] * u.deg
        )
        assert_quantity_allclose(
            w.proj_plane_pixel_area(), 1.93085492e-10 * (u.deg * u.deg)
        )

        # Exercise the main transformation functions, mainly just for
        # coverage
        w.p4_pix2foc([0, 100, 200], [0, -100, 200], 0)
        w.det2im([0, 100, 200], [0, -100, 200], 0)

        w.cpdis1 = w.cpdis1
        w.cpdis2 = w.cpdis2

        w.det2im1 = w.det2im1
        w.det2im2 = w.det2im2

        w.sip = w.sip

        w.cpdis1.cdelt = w.cpdis1.cdelt
        w.cpdis1.crpix = w.cpdis1.crpix
        w.cpdis1.crval = w.cpdis1.crval
        w.cpdis1.data = w.cpdis1.data

        assert w.sip.a_order == 4
        assert w.sip.b_order == 4
        assert w.sip.ap_order == 0
        assert w.sip.bp_order == 0
        assert_array_equal(w.sip.crpix, [2048.0, 1024.0])
        wcs.WCS(hdulist[1].header, hdulist)


def test_cpdis_comments():
    path = get_pkg_data_filename("data/dist_lookup.fits.gz")

    f = fits.open(path)
    w = wcs.WCS(f[1].header, f)
    hdr = w.to_fits()[0].header
    f.close()

    wcscards = list(hdr["CPDIS*"].cards) + list(hdr["DP*"].cards)
    wcsdict = {k: (v, c) for k, v, c in wcscards}

    refcards = [
        ("CPDIS1", "LOOKUP", "Prior distortion function type"),
        ("DP1.EXTVER", 1.0, "Version number of WCSDVARR extension"),
        ("DP1.NAXES", 2.0, "Number of independent variables in CPDIS function"),
        ("DP1.AXIS.1", 1.0, "Axis number of the 1st variable in a CPDIS function"),
        ("DP1.AXIS.2", 2.0, "Axis number of the 2nd variable in a CPDIS function"),
        ("CPDIS2", "LOOKUP", "Prior distortion function type"),
        ("DP2.EXTVER", 2.0, "Version number of WCSDVARR extension"),
        ("DP2.NAXES", 2.0, "Number of independent variables in CPDIS function"),
        ("DP2.AXIS.1", 1.0, "Axis number of the 1st variable in a CPDIS function"),
        ("DP2.AXIS.2", 2.0, "Axis number of the 2nd variable in a CPDIS function"),
    ]

    assert len(wcsdict) == len(refcards)

    for k, v, c in refcards:
        assert wcsdict[k] == (v, c)


def test_d2im_comments():
    path = get_pkg_data_filename("data/ie6d07ujq_wcs.fits")

    f = fits.open(path)
    with pytest.warns(wcs.FITSFixedWarning):
        w = wcs.WCS(f[0].header, f)
    f.close()
    wcscards = list(w.to_fits()[0].header["D2IM*"].cards)
    wcsdict = {k: (v, c) for k, v, c in wcscards}

    refcards = [
        ("D2IMDIS1", "LOOKUP", "Detector to image correction type"),
        ("D2IM1.EXTVER", 1.0, "Version number of WCSDVARR extension"),
        ("D2IM1.NAXES", 2.0, "Number of independent variables in D2IM function"),
        ("D2IM1.AXIS.1", 1.0, "Axis number of the 1st variable in a D2IM function"),
        ("D2IM1.AXIS.2", 2.0, "Axis number of the 2nd variable in a D2IM function"),
        ("D2IMDIS2", "LOOKUP", "Detector to image correction type"),
        ("D2IM2.EXTVER", 2.0, "Version number of WCSDVARR extension"),
        ("D2IM2.NAXES", 2.0, "Number of independent variables in D2IM function"),
        ("D2IM2.AXIS.1", 1.0, "Axis number of the 1st variable in a D2IM function"),
        ("D2IM2.AXIS.2", 2.0, "Axis number of the 2nd variable in a D2IM function"),
        # ('D2IMERR1', 0.049, 'Maximum error of D2IM correction for axis 1'),
        # ('D2IMERR2', 0.035, 'Maximum error of D2IM correction for axis 2'),
        # ('D2IMEXT', 'iref$y7b1516hi_d2i.fits', ''),
    ]

    assert len(wcsdict) == len(refcards)

    for k, v, c in refcards:
        assert wcsdict[k] == (v, c)


def test_sip_broken():
    # This header caused wcslib to segfault because it has a SIP
    # specification in a non-default keyword
    hdr = get_pkg_data_contents("data/sip-broken.hdr")

    wcs.WCS(hdr)


def test_no_truncate_crval():
    """
    Regression test for https://github.com/astropy/astropy/issues/4612
    """
    w = wcs.WCS(naxis=3)
    w.wcs.crval = [50, 50, 2.12345678e11]
    w.wcs.cdelt = [1e-3, 1e-3, 1e8]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "FREQ"]
    w.wcs.set()

    header = w.to_header()
    for ii in range(3):
        assert header[f"CRVAL{ii + 1}"] == w.wcs.crval[ii]
        assert header[f"CDELT{ii + 1}"] == w.wcs.cdelt[ii]


def test_no_truncate_crval_try2():
    """
    Regression test for https://github.com/astropy/astropy/issues/4612
    """
    w = wcs.WCS(naxis=3)
    w.wcs.crval = [50, 50, 2.12345678e11]
    w.wcs.cdelt = [1e-5, 1e-5, 1e5]
    w.wcs.ctype = ["RA---SIN", "DEC--SIN", "FREQ"]
    w.wcs.cunit = ["deg", "deg", "Hz"]
    w.wcs.crpix = [1, 1, 1]
    w.wcs.restfrq = 2.34e11
    w.wcs.set()

    header = w.to_header()
    for ii in range(3):
        assert header[f"CRVAL{ii + 1}"] == w.wcs.crval[ii]
        assert header[f"CDELT{ii + 1}"] == w.wcs.cdelt[ii]


def test_no_truncate_crval_p17():
    """
    Regression test for https://github.com/astropy/astropy/issues/5162
    """
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [50.1234567890123456, 50.1234567890123456]
    w.wcs.cdelt = [1e-3, 1e-3]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.set()

    header = w.to_header()
    assert header["CRVAL1"] != w.wcs.crval[0]
    assert header["CRVAL2"] != w.wcs.crval[1]
    header = w.to_header(relax=wcs.WCSHDO_P17)
    assert header["CRVAL1"] == w.wcs.crval[0]
    assert header["CRVAL2"] == w.wcs.crval[1]


def test_no_truncate_using_compare():
    """
    Regression test for https://github.com/astropy/astropy/issues/4612

    This one uses WCS.wcs.compare and some slightly different values
    """
    w = wcs.WCS(naxis=3)
    w.wcs.crval = [2.409303333333e02, 50, 2.12345678e11]
    w.wcs.cdelt = [1e-3, 1e-3, 1e8]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "FREQ"]
    w.wcs.set()
    w2 = wcs.WCS(w.to_header())
    w.wcs.compare(w2.wcs)


def test_passing_ImageHDU():
    """
    Passing ImageHDU or PrimaryHDU and comparing it with
    wcs initialized from header. For #4493.
    """
    path = get_pkg_data_filename("data/validate.fits")
    with fits.open(path) as hdulist:
        with pytest.warns(wcs.FITSFixedWarning):
            wcs_hdu = wcs.WCS(hdulist[0])
            wcs_header = wcs.WCS(hdulist[0].header)
            assert wcs_hdu.wcs.compare(wcs_header.wcs)
            wcs_hdu = wcs.WCS(hdulist[1])
            wcs_header = wcs.WCS(hdulist[1].header)
            assert wcs_hdu.wcs.compare(wcs_header.wcs)


def test_inconsistent_sip():
    """
    Test for #4814
    """
    hdr = get_pkg_data_contents("data/sip-broken.hdr")
    ctx = ctx_for_v71_dateref_warnings()
    with ctx:
        w = wcs.WCS(hdr)
    with pytest.warns(AstropyWarning):
        newhdr = w.to_header(relax=None)
    # CTYPE should not include "-SIP" if relax is None
    with ctx:
        wnew = wcs.WCS(newhdr)
    assert all(not ctyp.endswith("-SIP") for ctyp in wnew.wcs.ctype)
    newhdr = w.to_header(relax=False)
    assert "A_0_2" not in newhdr
    # CTYPE should not include "-SIP" if relax is False
    with ctx:
        wnew = wcs.WCS(newhdr)
    assert all(not ctyp.endswith("-SIP") for ctyp in wnew.wcs.ctype)
    with pytest.warns(AstropyWarning):
        newhdr = w.to_header(key="C")
    assert "A_0_2" not in newhdr
    # Test writing header with a different key
    with ctx:
        wnew = wcs.WCS(newhdr, key="C")
    assert all(not ctyp.endswith("-SIP") for ctyp in wnew.wcs.ctype)
    with pytest.warns(AstropyWarning):
        newhdr = w.to_header(key=" ")
    # Test writing a primary WCS to header
    with ctx:
        wnew = wcs.WCS(newhdr)
    assert all(not ctyp.endswith("-SIP") for ctyp in wnew.wcs.ctype)
    # Test that "-SIP" is kept into CTYPE if relax=True and
    # "-SIP" was in the original header
    newhdr = w.to_header(relax=True)
    with ctx:
        wnew = wcs.WCS(newhdr)
    assert all(ctyp.endswith("-SIP") for ctyp in wnew.wcs.ctype)
    assert "A_0_2" in newhdr
    # Test that SIP coefficients are also written out.
    assert wnew.sip is not None
    # ######### broken header ###########
    # Test that "-SIP" is added to CTYPE if relax=True and
    # "-SIP" was not in the original header but SIP coefficients
    # are present.
    with ctx:
        w = wcs.WCS(hdr)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    newhdr = w.to_header(relax=True)
    with ctx:
        wnew = wcs.WCS(newhdr)
    assert all(ctyp.endswith("-SIP") for ctyp in wnew.wcs.ctype)


def test_bounds_check():
    """Test for #4957"""
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---CAR", "DEC--CAR"]
    w.wcs.cdelt = [10, 10]
    w.wcs.crval = [-90, 90]
    w.wcs.crpix = [1, 1]
    w.wcs.bounds_check(False, False)
    ra, dec = w.wcs_pix2world(300, 0, 0)
    assert_allclose(ra, -180)
    assert_allclose(dec, -30)


def test_naxis():
    w = wcs.WCS(naxis=2)
    w.wcs.crval = [1, 1]
    w.wcs.cdelt = [0.1, 0.1]
    w.wcs.crpix = [1, 1]
    w._naxis = [1000, 500]
    assert w.pixel_shape == (1000, 500)
    assert w.array_shape == (500, 1000)

    w.pixel_shape = (99, 59)
    assert w._naxis == [99, 59]

    w.array_shape = (45, 23)
    assert w._naxis == [23, 45]
    assert w.pixel_shape == (23, 45)

    w.pixel_shape = None
    assert w.pixel_bounds is None


def test_sip_with_altkey():
    """
    Test that when creating a WCS object using a key, CTYPE with
    that key is looked at and not the primary CTYPE.
    fix for #5443.
    """
    with fits.open(get_pkg_data_filename("data/sip.fits")) as f:
        with pytest.warns(wcs.FITSFixedWarning):
            w = wcs.WCS(f[0].header)
    # create a header with two WCSs.
    h1 = w.to_header(relax=True, key="A")
    h2 = w.to_header(relax=False)
    h1["CTYPE1A"] = "RA---SIN-SIP"
    h1["CTYPE2A"] = "DEC--SIN-SIP"
    h1.update(h2)
    with ctx_for_v71_dateref_warnings():
        w = wcs.WCS(h1, key="A")
    assert (w.wcs.ctype == np.array(["RA---SIN-SIP", "DEC--SIN-SIP"])).all()


def test_to_fits_1():
    """
    Test to_fits() with LookupTable distortion.
    """
    fits_name = get_pkg_data_filename("data/dist.fits")

    with (
        pytest.warns(AstropyDeprecationWarning),
        pytest.warns(
            wcs.FITSFixedWarning,
            match="The WCS transformation has more axes",
        ),
    ):
        w = wcs.WCS(fits_name)
    wfits = w.to_fits()
    assert isinstance(wfits, fits.HDUList)
    assert isinstance(wfits[0], fits.PrimaryHDU)
    assert isinstance(wfits[1], fits.ImageHDU)


def test_keyedsip():
    """
    Test sip reading with extra key.
    """
    hdr_name = get_pkg_data_filename("data/sip-broken.hdr")
    header = fits.Header.fromfile(hdr_name)
    del header["CRPIX1"]
    del header["CRPIX2"]

    w = wcs.WCS(header=header, key="A")
    assert isinstance(w.sip, wcs.Sip)
    assert w.sip.crpix[0] == 2048
    assert w.sip.crpix[1] == 1026


def test_zero_size_input():
    with fits.open(get_pkg_data_filename("data/sip.fits")) as f:
        with pytest.warns(wcs.FITSFixedWarning):
            w = wcs.WCS(f[0].header)

    inp = np.zeros((0, 2))
    assert_array_equal(inp, w.all_pix2world(inp, 0))
    assert_array_equal(inp, w.all_world2pix(inp, 0))

    inp = [], [1]
    result = w.all_pix2world([], [1], 0)
    assert_array_equal(inp[0], result[0])
    assert_array_equal(inp[1], result[1])

    result = w.all_world2pix([], [1], 0)
    assert_array_equal(inp[0], result[0])
    assert_array_equal(inp[1], result[1])


def test_scalar_inputs():
    """
    Issue #7845
    """
    wcsobj = wcs.WCS(naxis=1)
    result = wcsobj.all_pix2world(2, 1)
    assert_array_equal(result, [np.array(2.0)])
    assert result[0].shape == ()

    result = wcsobj.all_pix2world([2], 1)
    assert_array_equal(result, [np.array([2.0])])
    assert result[0].shape == (1,)


# Ignore RuntimeWarning raised on s390.
@pytest.mark.filterwarnings("ignore:.*invalid value encountered in.*")
def test_footprint_contains():
    """
    Test WCS.footprint_contains(skycoord)
    """

    header = """
WCSAXES =                    2 / Number of coordinate axes
CRPIX1  =               1045.0 / Pixel coordinate of reference point
CRPIX2  =               1001.0 / Pixel coordinate of reference point
PC1_1   =    -0.00556448550786 / Coordinate transformation matrix element
PC1_2   =   -0.001042120133257 / Coordinate transformation matrix element
PC2_1   =    0.001181477028705 / Coordinate transformation matrix element
PC2_2   =   -0.005590809742987 / Coordinate transformation matrix element
CDELT1  =                  1.0 / [deg] Coordinate increment at reference point
CDELT2  =                  1.0 / [deg] Coordinate increment at reference point
CUNIT1  = 'deg'                / Units of coordinate increment and value
CUNIT2  = 'deg'                / Units of coordinate increment and value
CTYPE1  = 'RA---TAN'           / TAN (gnomonic) projection + SIP distortions
CTYPE2  = 'DEC--TAN'           / TAN (gnomonic) projection + SIP distortions
CRVAL1  =      250.34971683647 / [deg] Coordinate value at reference point
CRVAL2  =      2.2808772582495 / [deg] Coordinate value at reference point
LONPOLE =                180.0 / [deg] Native longitude of celestial pole
LATPOLE =      2.2808772582495 / [deg] Native latitude of celestial pole
RADESYS = 'ICRS'               / Equatorial coordinate system
MJD-OBS =      58612.339199259 / [d] MJD of observation matching DATE-OBS
DATE-OBS= '2019-05-09T08:08:26.816Z' / ISO-8601 observation date matching MJD-OB
NAXIS   =                    2 / NAXIS
NAXIS1  =                 2136 / length of first array dimension
NAXIS2  =                 2078 / length of second array dimension
    """

    header = fits.Header.fromstring(header.strip(), "\n")
    test_wcs = wcs.WCS(header)

    hasCoord = test_wcs.footprint_contains(SkyCoord(254, 2, unit="deg"))
    assert hasCoord

    hasCoord = test_wcs.footprint_contains(SkyCoord(240, 2, unit="deg"))
    assert not hasCoord

    hasCoord = test_wcs.footprint_contains(SkyCoord(24, 2, unit="deg"))
    assert not hasCoord


def test_cunit():
    # Initializing WCS
    w1 = wcs.WCS(naxis=2)
    w2 = wcs.WCS(naxis=2)
    w3 = wcs.WCS(naxis=2)
    w4 = wcs.WCS(naxis=2)
    # Initializing the values of cunit
    w1.wcs.cunit = ["deg", "m/s"]
    w2.wcs.cunit = ["km/h", "km/h"]
    w3.wcs.cunit = ["deg", "m/s"]
    w4.wcs.cunit = ["deg", "deg"]

    # Equality checking a cunit with itself
    assert w1.wcs.cunit == w1.wcs.cunit
    assert not w1.wcs.cunit != w1.wcs.cunit
    # Equality checking of two different cunit object having same values
    assert w1.wcs.cunit == w3.wcs.cunit
    assert not w1.wcs.cunit != w3.wcs.cunit
    # Equality checking of two different cunit object having the same first unit
    # but different second unit (see #9154)
    assert not w1.wcs.cunit == w4.wcs.cunit
    assert w1.wcs.cunit != w4.wcs.cunit
    # Inequality checking of two different cunit object having different values
    assert not w1.wcs.cunit == w2.wcs.cunit
    assert w1.wcs.cunit != w2.wcs.cunit
    # Inequality checking of cunit with a list of literals
    assert not w1.wcs.cunit == [1, 2, 3]
    assert w1.wcs.cunit != [1, 2, 3]
    # Inequality checking with some characters
    assert not w1.wcs.cunit == ["a", "b", "c"]
    assert w1.wcs.cunit != ["a", "b", "c"]
    # Comparison is not implemented TypeError will raise
    with pytest.raises(TypeError):
        w1.wcs.cunit < w2.wcs.cunit  # noqa: B015


class TestWcsWithTime:
    def setup_method(self):
        if _WCSLIB_VER >= Version("7.1"):
            fname = get_pkg_data_filename("data/header_with_time_wcslib71.fits")
        else:
            fname = get_pkg_data_filename("data/header_with_time.fits")
        self.header = fits.Header.fromfile(fname)
        with pytest.warns(wcs.FITSFixedWarning):
            self.w = wcs.WCS(self.header, key="A")

    def test_keywods2wcsprm(self):
        """Make sure Wcsprm is populated correctly from the header."""

        ctype = [self.header[val] for val in self.header["CTYPE*"]]
        crval = [self.header[val] for val in self.header["CRVAL*"]]
        crpix = [self.header[val] for val in self.header["CRPIX*"]]
        cdelt = [self.header[val] for val in self.header["CDELT*"]]
        cunit = [self.header[val] for val in self.header["CUNIT*"]]
        assert list(self.w.wcs.ctype) == ctype
        time_axis_code = 4000 if _WCSLIB_VER >= Version("7.9") else 0
        assert list(self.w.wcs.axis_types) == [2200, 2201, 3300, time_axis_code]
        assert_allclose(self.w.wcs.crval, crval)
        assert_allclose(self.w.wcs.crpix, crpix)
        assert_allclose(self.w.wcs.cdelt, cdelt)
        assert list(self.w.wcs.cunit) == cunit

        naxis = self.w.naxis
        assert naxis == 4
        pc = np.zeros((naxis, naxis), dtype=np.float64)
        for i in range(1, 5):
            for j in range(1, 5):
                if i == j:
                    pc[i - 1, j - 1] = self.header.get(f"PC{i}_{j}A", 1)
                else:
                    pc[i - 1, j - 1] = self.header.get(f"PC{i}_{j}A", 0)
        assert_allclose(self.w.wcs.pc, pc)

        char_keys = [
            "timesys",
            "trefpos",
            "trefdir",
            "plephem",
            "timeunit",
            "dateref",
            "dateobs",
            "datebeg",
            "dateavg",
            "dateend",
        ]
        for key in char_keys:
            assert getattr(self.w.wcs, key) == self.header.get(key, "")

        num_keys = [
            "mjdref",
            "mjdobs",
            "mjdbeg",
            "mjdend",
            "jepoch",
            "bepoch",
            "tstart",
            "tstop",
            "xposure",
            "timsyer",
            "timrder",
            "timedel",
            "timepixr",
            "timeoffs",
            "telapse",
            "czphs",
            "cperi",
        ]

        for key in num_keys:
            if key.upper() == "MJDREF":
                hdrv = [
                    self.header.get("MJDREFIA", np.nan),
                    self.header.get("MJDREFFA", np.nan),
                ]
            else:
                hdrv = self.header.get(key, np.nan)
            assert_allclose(getattr(self.w.wcs, key), hdrv)

    def test_transforms(self):
        assert_allclose(self.w.all_pix2world(*self.w.wcs.crpix, 1), self.w.wcs.crval)


def test_invalid_coordinate_masking():
    # Regression test for an issue which caused all coordinates to be set to NaN
    # after a transformation rather than just the invalid ones as reported by
    # WCSLIB. A specific example of this is that when considering an all-sky
    # spectral cube with a spectral axis that is not correlated with the sky
    # axes, if transforming pixel coordinates that did not fall 'in' the sky,
    # the spectral world value was also masked even though that coordinate
    # was valid.

    w = wcs.WCS(naxis=3)
    w.wcs.ctype = "VELO_LSR", "GLON-CAR", "GLAT-CAR"
    w.wcs.crval = -20, 0, 0
    w.wcs.crpix = 1, 1441, 241
    w.wcs.cdelt = 1.3, -0.125, 0.125

    px = [-10, -10, 20]
    py = [-10, 10, 20]
    pz = [-10, 10, 20]

    wx, wy, wz = w.wcs_pix2world(px, py, pz, 0)

    # Before fixing this, wx used to return np.nan for the first element

    assert_allclose(wx, [-33, -33, 6])
    assert_allclose(wy, [np.nan, 178.75, 177.5])
    assert_allclose(wz, [np.nan, -28.75, -27.5])


def test_no_pixel_area():
    w = wcs.WCS(naxis=3)

    # Pixel area cannot be computed
    with pytest.raises(ValueError, match="Pixel area is defined only for 2D pixels"):
        w.proj_plane_pixel_area()

    # Pixel scales still possible
    assert_quantity_allclose(w.proj_plane_pixel_scales(), 1)


def test_distortion_header(tmp_path):
    """
    Test that plate distortion model is correctly described by `wcs.to_header()`
    and preserved when creating a Cutout2D from the image, writing it to FITS,
    and reading it back from the file.
    """
    path = get_pkg_data_filename("data/dss.14.29.56-62.41.05.fits.gz")
    cen = np.array((50, 50))
    size = np.array((20, 20))

    with fits.open(path) as hdulist:
        with pytest.warns(VerifyWarning), pytest.warns(wcs.FITSFixedWarning):
            w = wcs.WCS(hdulist[0].header)
        cut = Cutout2D(hdulist[0].data, position=cen, size=size, wcs=w)

    # This converts the DSS plate solution model with AMD[XY]n coefficients into a
    # Template Polynomial Distortion model (TPD.FWD.n coefficients);
    # not testing explicitly for the header keywords here.

    if _WCSLIB_VER < Version("7.4"):
        with pytest.warns(
            AstropyWarning, match="WCS contains a TPD distortion model in CQDIS"
        ):
            w0 = wcs.WCS(w.to_header_string())
        with pytest.warns(
            AstropyWarning, match="WCS contains a TPD distortion model in CQDIS"
        ):
            w1 = wcs.WCS(cut.wcs.to_header_string())
        if _WCSLIB_VER >= Version("7.1"):
            pytest.xfail("TPD coefficients incomplete with WCSLIB >= 7.1 < 7.4")
    else:
        w0 = wcs.WCS(w.to_header_string())
        w1 = wcs.WCS(cut.wcs.to_header_string())

    assert w.pixel_to_world(0, 0).separation(w0.pixel_to_world(0, 0)) < 1.0e-3 * u.mas
    assert w.pixel_to_world(*cen).separation(w0.pixel_to_world(*cen)) < 1.0e-3 * u.mas

    assert (
        w.pixel_to_world(*cen).separation(w1.pixel_to_world(*(size / 2)))
        < 1.0e-3 * u.mas
    )

    cutfile = tmp_path / "cutout.fits"
    fits.writeto(cutfile, cut.data, cut.wcs.to_header())

    with fits.open(cutfile) as hdulist:
        w2 = wcs.WCS(hdulist[0].header)

    assert (
        w.pixel_to_world(*cen).separation(w2.pixel_to_world(*(size / 2)))
        < 1.0e-3 * u.mas
    )


def test_pixlist_wcs_colsel():
    """
    Test selection of a specific pixel list WCS using ``colsel``. See #11412.
    """
    hdr_file = get_pkg_data_filename("data/chandra-pixlist-wcs.hdr")
    hdr = fits.Header.fromtextfile(hdr_file)
    with pytest.warns(wcs.FITSFixedWarning):
        w = wcs.WCS(hdr, keysel=["image", "pixel"], colsel=[11, 12])
    assert w.naxis == 2
    assert list(w.wcs.ctype) == ["RA---TAN", "DEC--TAN"]
    assert np.allclose(w.wcs.crval, [229.38051931869, -58.81108068885])
    assert np.allclose(w.wcs.pc, [[1, 0], [0, 1]])
    assert np.allclose(w.wcs.cdelt, [-0.00013666666666666, 0.00013666666666666])
    assert np.allclose(w.wcs.lonpole, 180.0)


@pytest.mark.skipif(
    _WCSLIB_VER < Version("7.8"),
    reason="TIME axis extraction only works with wcslib 7.8 or later",
)
def test_time_axis_selection():
    w = wcs.WCS(naxis=3)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "TIME"]
    w.wcs.set()
    assert list(w.sub([wcs.WCSSUB_TIME]).wcs.ctype) == ["TIME"]
    assert (
        w.wcs_pix2world([[1, 2, 3]], 0)[0, 2]
        == w.sub([wcs.WCSSUB_TIME]).wcs_pix2world([[3]], 0)[0, 0]
    )


@pytest.mark.skipif(
    _WCSLIB_VER < Version("7.8"),
    reason="TIME axis extraction only works with wcslib 7.8 or later",
)
def test_temporal():
    w = wcs.WCS(naxis=3)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "TIME"]
    w.wcs.set()
    assert w.has_temporal
    assert w.sub([wcs.WCSSUB_TIME]).is_temporal
    assert (
        w.wcs_pix2world([[1, 2, 3]], 0)[0, 2]
        == w.temporal.wcs_pix2world([[3]], 0)[0, 0]
    )


def test_swapaxes_same_val_roundtrip():
    w = wcs.WCS(naxis=3)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "FREQ"]
    w.wcs.crpix = [32.5, 16.5, 1.0]
    w.wcs.crval = [5.63, -72.05, 1.0]
    w.wcs.pc = [[5.9e-06, 1.3e-05, 0.0], [-1.2e-05, 5.0e-06, 0.0], [0.0, 0.0, 1.0]]
    w.wcs.cdelt = [1.0, 1.0, 1.0]
    w.wcs.set()
    axes_order = [3, 2, 1]
    axes_order0 = [i - 1 for i in axes_order]
    ws = w.sub(axes_order)
    imcoord = np.array([3, 5, 7])
    imcoords = imcoord[axes_order0]
    val_ref = w.wcs_pix2world([imcoord], 0)[0]
    val_swapped = ws.wcs_pix2world([imcoords], 0)[0]

    # check original axis and swapped give same results
    assert np.allclose(val_ref[axes_order0], val_swapped, rtol=0, atol=1e-8)

    # check round-tripping:
    assert np.allclose(w.wcs_world2pix([val_ref], 0)[0], imcoord, rtol=0, atol=1e-8)


def test_DistortionLookupTable():
    img_world_wcs = wcs.WCS(naxis=2)
    # A simple "pixel coordinates are world coordinates" WCS, to which we'll
    # add distortion lookup tables
    img_world_wcs.wcs.crpix = 1, 1
    img_world_wcs.wcs.crval = 0, 0
    img_world_wcs.wcs.cdelt = 1, 1

    # Create maps with zero distortion except at one particular pixel
    x_dist_array = np.zeros((25, 25))
    x_dist_array[10, 20] = 0.5
    map_x = wcs.DistortionLookupTable(
        x_dist_array.astype(np.float32), (5, 10), (10, 20), (2, 2)
    )

    y_dist_array = np.zeros((25, 25))
    y_dist_array[10, 5] = 0.7
    map_y = wcs.DistortionLookupTable(
        y_dist_array.astype(np.float32), (5, 10), (10, 20), (3, 3)
    )

    img_world_wcs.cpdis1 = map_x
    img_world_wcs.cpdis2 = map_y

    # The x distortion of 0.5 pixels should appear at a specific spot in the
    # image, and we need to work out what that is so we can check it. The
    # distortion is at array index (10, 20), which is a 1-based pixel
    # coordinate of (21, 11) and therefore at an offset of (16, 1) from the
    # lookup table's CRPIX of (5, 10). CDELT is 2 image pixels / distortion
    # pixel, so the distortion applies to the image pixel which is at an offset
    # of (32, 2) from the CRVAL of (10, 20), meaning we should see the
    # distortion at the 1-based image pixel coordinate of (42, 22).
    assert_allclose(map_x.get_offset(42, 22), 0.5)
    # And we should see it applied at the 0-based coordinate (41, 21) with our
    # simple img<->world wcs.
    assert_allclose(img_world_wcs.pixel_to_world_values(41, 21), [41.5, 21])
    # Similarly for the y distortion, the distortion is at array index (10, 5),
    # which is a 1-based pixel coordinate of (6, 11) and therefore at an offset
    # of (1, 1) from the lookup table's CRPIX of (5, 10). CDELT is 3 image
    # pixels / distortion pixel, so the distortion applies to the image pixel
    # which is at an offset of (3, 3) from the CRVAL of (10, 20), meaning we
    # should see the distortion at the 1-based image pixel coordinate of (13,
    # 23).
    assert_allclose(map_y.get_offset(13, 23), 0.7)
    # And we should see it applied at the 0-based coordinate (12, 22) with our
    # simple img<->world wcs.
    assert_allclose(img_world_wcs.pixel_to_world_values(12, 22), [12, 22.7])

    # Now check that when we move the image location by the equivalent of 1/2
    # distortion-array pixel, we see only half the distortion.
    for dx, dy in [(0.5, 0), (-0.5, 0), (0, 0.5), (0, -0.5)]:
        # Scale dx, dy by 2 for the CDELT in the x distortion table (since
        # we're looking for the x distortion).
        assert_allclose(
            img_world_wcs.pixel_to_world_values(41 + dx * 2, 21 + dy * 2),
            [41 + dx * 2 + 0.25, 21 + dy * 2],
        )
        # Scale dx, dy by 3 for the CDELT in the y distortion table (since
        # we're looking for the y distortion).
        assert_allclose(
            img_world_wcs.pixel_to_world_values(12 + dx * 3, 22 + dy * 3),
            [12 + dx * 3, 22 + dy * 3 + 0.35],
        )

    # Now check that when we move the image location by the equivalent of 1
    # distortion-array pixel, we see no distortion.
    for dx, dy in [(2, 0), (-2, 0), (0, 2), (0, -2)]:
        # Scale dx, dy by 2 for the CDELT in the x distortion table (since
        # we're looking for the x distortion).
        assert_allclose(
            img_world_wcs.pixel_to_world_values(41 + dx * 2, 21 + dy * 2),
            [41 + dx * 2, 21 + dy * 2],
        )
        # Scale dx, dy by 3 for the CDELT in the y distortion table (since
        # we're looking for the y distortion).
        assert_allclose(
            img_world_wcs.pixel_to_world_values(12 + dx * 3, 22 + dy * 3),
            [12 + dx * 3, 22 + dy * 3],
        )
