# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_allclose

from ...utils.data import get_pkg_data_contents, get_pkg_data_filename
from ...time import Time
from ... import units as u

from ..wcs import WCS, Sip, WCSSUB_LONGITUDE, WCSSUB_LATITUDE
from ..utils import (proj_plane_pixel_scales, proj_plane_pixel_area,
                     is_proj_plane_distorted,
                     non_celestial_pixel_scales, wcs_to_celestial_frame,
                     celestial_frame_to_wcs, skycoord_to_pixel,
                     pixel_to_skycoord, custom_wcs_to_frame_mappings,
                     custom_frame_to_wcs_mappings, add_stokes_axis_to_wcs)


def test_wcs_dropping():
    wcs = WCS(naxis=4)
    wcs.wcs.pc = np.zeros([4, 4])
    np.fill_diagonal(wcs.wcs.pc, np.arange(1, 5))
    pc = wcs.wcs.pc  # for later use below

    dropped = wcs.dropaxis(0)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([2, 3, 4]))
    dropped = wcs.dropaxis(1)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1, 3, 4]))
    dropped = wcs.dropaxis(2)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1, 2, 4]))
    dropped = wcs.dropaxis(3)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1, 2, 3]))

    wcs = WCS(naxis=4)
    wcs.wcs.cd = pc

    dropped = wcs.dropaxis(0)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([2, 3, 4]))
    dropped = wcs.dropaxis(1)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1, 3, 4]))
    dropped = wcs.dropaxis(2)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1, 2, 4]))
    dropped = wcs.dropaxis(3)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1, 2, 3]))


def test_wcs_swapping():
    wcs = WCS(naxis=4)
    wcs.wcs.pc = np.zeros([4, 4])
    np.fill_diagonal(wcs.wcs.pc, np.arange(1, 5))
    pc = wcs.wcs.pc  # for later use below

    swapped = wcs.swapaxes(0, 1)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([2, 1, 3, 4]))
    swapped = wcs.swapaxes(0, 3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([4, 2, 3, 1]))
    swapped = wcs.swapaxes(2, 3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([1, 2, 4, 3]))

    wcs = WCS(naxis=4)
    wcs.wcs.cd = pc

    swapped = wcs.swapaxes(0, 1)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([2, 1, 3, 4]))
    swapped = wcs.swapaxes(0, 3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([4, 2, 3, 1]))
    swapped = wcs.swapaxes(2, 3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([1, 2, 4, 3]))


@pytest.mark.parametrize('ndim', (2, 3))
def test_add_stokes(ndim):
    wcs = WCS(naxis=ndim)

    for ii in range(ndim + 1):
        outwcs = add_stokes_axis_to_wcs(wcs, ii)
        assert outwcs.wcs.naxis == ndim + 1
        assert outwcs.wcs.ctype[ii] == 'STOKES'
        assert outwcs.wcs.cname[ii] == 'STOKES'


def test_slice():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1, 1]
    mywcs.wcs.cdelt = [0.1, 0.1]
    mywcs.wcs.crpix = [1, 1]
    mywcs._naxis = [1000, 500]
    pscale = 0.1 # from cdelt

    slice_wcs = mywcs.slice([slice(1, None), slice(0, None)])
    assert np.all(slice_wcs.wcs.crpix == np.array([1, 0]))
    assert slice_wcs._naxis == [1000, 499]

    # test that CRPIX maps to CRVAL:
    assert_allclose(
        slice_wcs.wcs_pix2world(*slice_wcs.wcs.crpix, 1),
        slice_wcs.wcs.crval, rtol=0.0, atol=1e-6 * pscale
    )

    slice_wcs = mywcs.slice([slice(1, None, 2), slice(0, None, 4)])
    assert np.all(slice_wcs.wcs.crpix == np.array([0.625, 0.25]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4, 0.2]))
    assert slice_wcs._naxis == [250, 250]

    slice_wcs = mywcs.slice([slice(None, None, 2), slice(0, None, 2)])
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2, 0.2]))
    assert slice_wcs._naxis == [500, 250]

    # Non-integral values do not alter the naxis attribute
    slice_wcs = mywcs.slice([slice(50.), slice(20.)])
    assert slice_wcs._naxis == [1000, 500]
    slice_wcs = mywcs.slice([slice(50.), slice(20)])
    assert slice_wcs._naxis == [20, 500]
    slice_wcs = mywcs.slice([slice(50), slice(20.5)])
    assert slice_wcs._naxis == [1000, 50]


def test_slice_with_sip():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1, 1]
    mywcs.wcs.cdelt = [0.1, 0.1]
    mywcs.wcs.crpix = [1, 1]
    mywcs._naxis = [1000, 500]
    mywcs.wcs.ctype = ['RA---TAN-SIP', 'DEC--TAN-SIP']
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
    mywcs.sip = Sip(a, b, None, None, mywcs.wcs.crpix)
    mywcs.wcs.set()
    pscale = 0.1 # from cdelt

    slice_wcs = mywcs.slice([slice(1, None), slice(0, None)])
    # test that CRPIX maps to CRVAL:
    assert_allclose(
        slice_wcs.all_pix2world(*slice_wcs.wcs.crpix, 1),
        slice_wcs.wcs.crval, rtol=0.0, atol=1e-6 * pscale
    )

    slice_wcs = mywcs.slice([slice(1, None, 2), slice(0, None, 4)])
    # test that CRPIX maps to CRVAL:
    assert_allclose(
        slice_wcs.all_pix2world(*slice_wcs.wcs.crpix, 1),
        slice_wcs.wcs.crval, rtol=0.0, atol=1e-6 * pscale
    )


def test_slice_getitem():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1, 1]
    mywcs.wcs.cdelt = [0.1, 0.1]
    mywcs.wcs.crpix = [1, 1]

    slice_wcs = mywcs[1::2, 0::4]
    assert np.all(slice_wcs.wcs.crpix == np.array([0.625, 0.25]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4, 0.2]))

    mywcs.wcs.crpix = [2, 2]
    slice_wcs = mywcs[1::2, 0::4]
    assert np.all(slice_wcs.wcs.crpix == np.array([0.875, 0.75]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4, 0.2]))

    # Default: numpy order
    slice_wcs = mywcs[1::2]
    assert np.all(slice_wcs.wcs.crpix == np.array([2, 0.75]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.1, 0.2]))


def test_slice_fitsorder():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1, 1]
    mywcs.wcs.cdelt = [0.1, 0.1]
    mywcs.wcs.crpix = [1, 1]

    slice_wcs = mywcs.slice([slice(1, None), slice(0, None)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0, 1]))

    slice_wcs = mywcs.slice([slice(1, None, 2), slice(0, None, 4)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0.25, 0.625]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2, 0.4]))

    slice_wcs = mywcs.slice([slice(1, None, 2)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0.25, 1]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2, 0.1]))


def test_invalid_slice():
    mywcs = WCS(naxis=2)

    with pytest.raises(ValueError) as exc:
        mywcs[0]
    assert exc.value.args[0] == ("Cannot downsample a WCS with indexing.  Use "
                                 "wcs.sub or wcs.dropaxis if you want to remove "
                                 "axes.")

    with pytest.raises(ValueError) as exc:
        mywcs[0, ::2]
    assert exc.value.args[0] == ("Cannot downsample a WCS with indexing.  Use "
                                 "wcs.sub or wcs.dropaxis if you want to remove "
                                 "axes.")


def test_axis_names():
    mywcs = WCS(naxis=4)
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'VOPT-LSR', 'STOKES']

    assert mywcs.axis_type_names == ['RA', 'DEC', 'VOPT', 'STOKES']

    mywcs.wcs.cname = ['RA', 'DEC', 'VOPT', 'STOKES']

    assert mywcs.axis_type_names == ['RA', 'DEC', 'VOPT', 'STOKES']


def test_celestial():
    mywcs = WCS(naxis=4)
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'VOPT', 'STOKES']
    cel = mywcs.celestial
    assert tuple(cel.wcs.ctype) == ('RA---TAN', 'DEC--TAN')
    assert cel.axis_type_names == ['RA', 'DEC']


def test_wcs_to_celestial_frame():

    # Import astropy.coordinates here to avoid circular imports
    from ...coordinates.builtin_frames import ICRS, ITRS, FK5, FK4, Galactic

    mywcs = WCS(naxis=2)
    mywcs.wcs.set()
    with pytest.raises(ValueError) as exc:
        assert wcs_to_celestial_frame(mywcs) is None
    assert exc.value.args[0] == "Could not determine celestial frame corresponding to the specified WCS object"

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['XOFFSET', 'YOFFSET']
    mywcs.wcs.set()
    with pytest.raises(ValueError):
        assert wcs_to_celestial_frame(mywcs) is None

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    mywcs.wcs.set()
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, ICRS)

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    mywcs.wcs.equinox = 1987.
    mywcs.wcs.set()
    print(mywcs.to_header())
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, FK5)
    assert frame.equinox == Time(1987., format='jyear')

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    mywcs.wcs.equinox = 1982
    mywcs.wcs.set()
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, FK4)
    assert frame.equinox == Time(1982., format='byear')

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['GLON-SIN', 'GLAT-SIN']
    mywcs.wcs.set()
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, Galactic)

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['TLON-CAR', 'TLAT-CAR']
    mywcs.wcs.dateobs = '2017-08-17T12:41:04.430'
    mywcs.wcs.set()
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, ITRS)
    assert frame.obstime == Time('2017-08-17T12:41:04.430')

    for equinox in [np.nan, 1987, 1982]:
        mywcs = WCS(naxis=2)
        mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        mywcs.wcs.radesys = 'ICRS'
        mywcs.wcs.equinox = equinox
        mywcs.wcs.set()
        frame = wcs_to_celestial_frame(mywcs)
        assert isinstance(frame, ICRS)

    # Flipped order
    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['DEC--TAN', 'RA---TAN']
    mywcs.wcs.set()
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, ICRS)

    # More than two dimensions
    mywcs = WCS(naxis=3)
    mywcs.wcs.ctype = ['DEC--TAN', 'VELOCITY', 'RA---TAN']
    mywcs.wcs.set()
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, ICRS)

    mywcs = WCS(naxis=3)
    mywcs.wcs.ctype = ['GLAT-CAR', 'VELOCITY', 'GLON-CAR']
    mywcs.wcs.set()
    frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, Galactic)


def test_wcs_to_celestial_frame_extend():

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['XOFFSET', 'YOFFSET']
    mywcs.wcs.set()
    with pytest.raises(ValueError):
        wcs_to_celestial_frame(mywcs)

    class OffsetFrame:
        pass

    def identify_offset(wcs):
        if wcs.wcs.ctype[0].endswith('OFFSET') and wcs.wcs.ctype[1].endswith('OFFSET'):
            return OffsetFrame()

    with custom_wcs_to_frame_mappings(identify_offset):
        frame = wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, OffsetFrame)

    # Check that things are back to normal after the context manager
    with pytest.raises(ValueError):
        wcs_to_celestial_frame(mywcs)


def test_celestial_frame_to_wcs():

    # Import astropy.coordinates here to avoid circular imports
    from ...coordinates import ICRS, ITRS, FK5, FK4, FK4NoETerms, Galactic, BaseCoordinateFrame

    class FakeFrame(BaseCoordinateFrame):
        pass

    frame = FakeFrame()
    with pytest.raises(ValueError) as exc:
        celestial_frame_to_wcs(frame)
    assert exc.value.args[0] == ("Could not determine WCS corresponding to "
                                 "the specified coordinate frame.")

    frame = ICRS()
    mywcs = celestial_frame_to_wcs(frame)
    mywcs.wcs.set()
    assert tuple(mywcs.wcs.ctype) == ('RA---TAN', 'DEC--TAN')
    assert mywcs.wcs.radesys == 'ICRS'
    assert np.isnan(mywcs.wcs.equinox)
    assert mywcs.wcs.lonpole == 180
    assert mywcs.wcs.latpole == 0

    frame = FK5(equinox='J1987')
    mywcs = celestial_frame_to_wcs(frame)
    assert tuple(mywcs.wcs.ctype) == ('RA---TAN', 'DEC--TAN')
    assert mywcs.wcs.radesys == 'FK5'
    assert mywcs.wcs.equinox == 1987.

    frame = FK4(equinox='B1982')
    mywcs = celestial_frame_to_wcs(frame)
    assert tuple(mywcs.wcs.ctype) == ('RA---TAN', 'DEC--TAN')
    assert mywcs.wcs.radesys == 'FK4'
    assert mywcs.wcs.equinox == 1982.

    frame = FK4NoETerms(equinox='B1982')
    mywcs = celestial_frame_to_wcs(frame)
    assert tuple(mywcs.wcs.ctype) == ('RA---TAN', 'DEC--TAN')
    assert mywcs.wcs.radesys == 'FK4-NO-E'
    assert mywcs.wcs.equinox == 1982.

    frame = Galactic()
    mywcs = celestial_frame_to_wcs(frame)
    assert tuple(mywcs.wcs.ctype) == ('GLON-TAN', 'GLAT-TAN')
    assert mywcs.wcs.radesys == ''
    assert np.isnan(mywcs.wcs.equinox)

    frame = Galactic()
    mywcs = celestial_frame_to_wcs(frame, projection='CAR')
    assert tuple(mywcs.wcs.ctype) == ('GLON-CAR', 'GLAT-CAR')
    assert mywcs.wcs.radesys == ''
    assert np.isnan(mywcs.wcs.equinox)

    frame = Galactic()
    mywcs = celestial_frame_to_wcs(frame, projection='CAR')
    mywcs.wcs.crval = [100, -30]
    mywcs.wcs.set()
    assert_allclose((mywcs.wcs.lonpole, mywcs.wcs.latpole), (180, 60))

    frame = ITRS(obstime=Time('2017-08-17T12:41:04.43'))
    mywcs = celestial_frame_to_wcs(frame, projection='CAR')
    assert tuple(mywcs.wcs.ctype) == ('TLON-CAR', 'TLAT-CAR')
    assert mywcs.wcs.radesys == 'ITRS'
    assert mywcs.wcs.dateobs == '2017-08-17T12:41:04.430'

    frame = ITRS()
    mywcs = celestial_frame_to_wcs(frame, projection='CAR')
    assert tuple(mywcs.wcs.ctype) == ('TLON-CAR', 'TLAT-CAR')
    assert mywcs.wcs.radesys == 'ITRS'
    assert mywcs.wcs.dateobs == Time('J2000').utc.isot


def test_celestial_frame_to_wcs_extend():

    class OffsetFrame:
        pass

    frame = OffsetFrame()

    with pytest.raises(ValueError):
        celestial_frame_to_wcs(frame)

    def identify_offset(frame, projection=None):
        if isinstance(frame, OffsetFrame):
            wcs = WCS(naxis=2)
            wcs.wcs.ctype = ['XOFFSET', 'YOFFSET']
            return wcs

    with custom_frame_to_wcs_mappings(identify_offset):
        mywcs = celestial_frame_to_wcs(frame)
    assert tuple(mywcs.wcs.ctype) == ('XOFFSET', 'YOFFSET')

    # Check that things are back to normal after the context manager
    with pytest.raises(ValueError):
        celestial_frame_to_wcs(frame)


def test_pixscale_nodrop():
    mywcs = WCS(naxis=2)
    mywcs.wcs.cdelt = [0.1, 0.2]
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    assert_almost_equal(proj_plane_pixel_scales(mywcs), (0.1, 0.2))

    mywcs.wcs.cdelt = [-0.1, 0.2]
    assert_almost_equal(proj_plane_pixel_scales(mywcs), (0.1, 0.2))


def test_pixscale_withdrop():
    mywcs = WCS(naxis=3)
    mywcs.wcs.cdelt = [0.1, 0.2, 1]
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'VOPT']
    assert_almost_equal(proj_plane_pixel_scales(mywcs.celestial), (0.1, 0.2))

    mywcs.wcs.cdelt = [-0.1, 0.2, 1]
    assert_almost_equal(proj_plane_pixel_scales(mywcs.celestial), (0.1, 0.2))


def test_pixscale_cd():
    mywcs = WCS(naxis=2)
    mywcs.wcs.cd = [[-0.1, 0], [0, 0.2]]
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    assert_almost_equal(proj_plane_pixel_scales(mywcs), (0.1, 0.2))


@pytest.mark.parametrize('angle',
                         (30, 45, 60, 75))
def test_pixscale_cd_rotated(angle):
    mywcs = WCS(naxis=2)
    rho = np.radians(angle)
    scale = 0.1
    mywcs.wcs.cd = [[scale * np.cos(rho), -scale * np.sin(rho)],
                    [scale * np.sin(rho), scale * np.cos(rho)]]
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    assert_almost_equal(proj_plane_pixel_scales(mywcs), (0.1, 0.1))


@pytest.mark.parametrize('angle',
                         (30, 45, 60, 75))
def test_pixscale_pc_rotated(angle):
    mywcs = WCS(naxis=2)
    rho = np.radians(angle)
    scale = 0.1
    mywcs.wcs.cdelt = [-scale, scale]
    mywcs.wcs.pc = [[np.cos(rho), -np.sin(rho)],
                    [np.sin(rho), np.cos(rho)]]
    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    assert_almost_equal(proj_plane_pixel_scales(mywcs), (0.1, 0.1))


@pytest.mark.parametrize(('cdelt', 'pc', 'pccd'),
                         (([0.1, 0.2], np.eye(2), np.diag([0.1, 0.2])),
                          ([0.1, 0.2, 0.3], np.eye(3), np.diag([0.1, 0.2, 0.3])),
                          ([1, 1, 1], np.diag([0.1, 0.2, 0.3]), np.diag([0.1, 0.2, 0.3]))))
def test_pixel_scale_matrix(cdelt, pc, pccd):

    mywcs = WCS(naxis=(len(cdelt)))
    mywcs.wcs.cdelt = cdelt
    mywcs.wcs.pc = pc

    assert_almost_equal(mywcs.pixel_scale_matrix, pccd)


@pytest.mark.parametrize(('ctype', 'cel'),
                         ((['RA---TAN', 'DEC--TAN'], True),
                          (['RA---TAN', 'DEC--TAN', 'FREQ'], False),
                          (['RA---TAN', 'FREQ'], False),))
def test_is_celestial(ctype, cel):
    mywcs = WCS(naxis=len(ctype))
    mywcs.wcs.ctype = ctype

    assert mywcs.is_celestial == cel


@pytest.mark.parametrize(('ctype', 'cel'),
                         ((['RA---TAN', 'DEC--TAN'], True),
                          (['RA---TAN', 'DEC--TAN', 'FREQ'], True),
                          (['RA---TAN', 'FREQ'], False),))
def test_has_celestial(ctype, cel):
    mywcs = WCS(naxis=len(ctype))
    mywcs.wcs.ctype = ctype

    assert mywcs.has_celestial == cel


@pytest.mark.parametrize(('cdelt', 'pc', 'cd'),
                         ((np.array([0.1, 0.2]), np.eye(2), np.eye(2)),
                          (np.array([1, 1]), np.diag([0.1, 0.2]), np.eye(2)),
                          (np.array([0.1, 0.2]), np.eye(2), None),
                          (np.array([0.1, 0.2]), None, np.eye(2)),
                          ))
def test_noncelestial_scale(cdelt, pc, cd):

    mywcs = WCS(naxis=2)
    if cd is not None:
        mywcs.wcs.cd = cd
    if pc is not None:
        mywcs.wcs.pc = pc
    mywcs.wcs.cdelt = cdelt

    mywcs.wcs.ctype = ['RA---TAN', 'FREQ']

    ps = non_celestial_pixel_scales(mywcs)

    assert_almost_equal(ps.to_value(u.deg), np.array([0.1, 0.2]))


@pytest.mark.parametrize('mode', ['all', 'wcs'])
def test_skycoord_to_pixel(mode):

    # Import astropy.coordinates here to avoid circular imports
    from ...coordinates import SkyCoord

    header = get_pkg_data_contents('maps/1904-66_TAN.hdr', encoding='binary')
    wcs = WCS(header)

    ref = SkyCoord(0.1 * u.deg, -89. * u.deg, frame='icrs')

    xp, yp = skycoord_to_pixel(ref, wcs, mode=mode)

    # WCS is in FK5 so we need to transform back to ICRS
    new = pixel_to_skycoord(xp, yp, wcs, mode=mode).transform_to('icrs')

    assert_allclose(new.ra.degree, ref.ra.degree)
    assert_allclose(new.dec.degree, ref.dec.degree)

    # Make sure you can specify a different class using ``cls`` keyword
    class SkyCoord2(SkyCoord):
        pass

    new2 = pixel_to_skycoord(xp, yp, wcs, mode=mode,
                             cls=SkyCoord2).transform_to('icrs')

    assert new2.__class__ is SkyCoord2
    assert_allclose(new2.ra.degree, ref.ra.degree)
    assert_allclose(new2.dec.degree, ref.dec.degree)


def test_skycoord_to_pixel_swapped():

    # Regression test for a bug that caused skycoord_to_pixel and
    # pixel_to_skycoord to not work correctly if the axes were swapped in the
    # WCS.

    # Import astropy.coordinates here to avoid circular imports
    from ...coordinates import SkyCoord

    header = get_pkg_data_contents('maps/1904-66_TAN.hdr', encoding='binary')
    wcs = WCS(header)

    wcs_swapped = wcs.sub([WCSSUB_LATITUDE, WCSSUB_LONGITUDE])

    ref = SkyCoord(0.1 * u.deg, -89. * u.deg, frame='icrs')

    xp1, yp1 = skycoord_to_pixel(ref, wcs)
    xp2, yp2 = skycoord_to_pixel(ref, wcs_swapped)

    assert_allclose(xp1, xp2)
    assert_allclose(yp1, yp2)

    # WCS is in FK5 so we need to transform back to ICRS
    new1 = pixel_to_skycoord(xp1, yp1, wcs).transform_to('icrs')
    new2 = pixel_to_skycoord(xp1, yp1, wcs_swapped).transform_to('icrs')

    assert_allclose(new1.ra.degree, new2.ra.degree)
    assert_allclose(new1.dec.degree, new2.dec.degree)


def test_is_proj_plane_distorted():
    # non-orthogonal CD:
    wcs = WCS(naxis=2)
    wcs.wcs.cd = [[-0.1, 0], [0, 0.2]]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    assert(is_proj_plane_distorted(wcs))

    # almost orthogonal CD:
    wcs.wcs.cd = [[0.1 + 2.0e-7, 1.7e-7], [1.2e-7, 0.1 - 1.3e-7]]
    assert(not is_proj_plane_distorted(wcs))

    # real case:
    header = get_pkg_data_filename('data/sip.fits')
    wcs = WCS(header)
    assert(is_proj_plane_distorted(wcs))


@pytest.mark.parametrize('mode', ['all', 'wcs'])
def test_skycoord_to_pixel_distortions(mode):

    # Import astropy.coordinates here to avoid circular imports
    from ...coordinates import SkyCoord

    header = get_pkg_data_filename('data/sip.fits')
    wcs = WCS(header)

    ref = SkyCoord(202.50 * u.deg, 47.19 * u.deg, frame='icrs')

    xp, yp = skycoord_to_pixel(ref, wcs, mode=mode)

    # WCS is in FK5 so we need to transform back to ICRS
    new = pixel_to_skycoord(xp, yp, wcs, mode=mode).transform_to('icrs')

    assert_allclose(new.ra.degree, ref.ra.degree)
    assert_allclose(new.dec.degree, ref.dec.degree)
