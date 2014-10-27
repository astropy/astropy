# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

from ...utils.data import get_pkg_data_contents, get_pkg_data_filename
from ...wcs import WCS
from .. import utils
from ..utils import celestial_pixel_scale, non_celestial_pixel_scales
from ...tests.helper import pytest, catch_warnings
from ...utils.exceptions import AstropyUserWarning
from ... import units as u

import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_allclose

def test_wcs_dropping():
    wcs = WCS(naxis=4)
    wcs.wcs.pc = np.zeros([4,4])
    np.fill_diagonal(wcs.wcs.pc, np.arange(1,5))
    pc = wcs.wcs.pc # for later use below

    dropped = wcs.dropaxis(0)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([2,3,4]))
    dropped = wcs.dropaxis(1)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,3,4]))
    dropped = wcs.dropaxis(2)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,4]))
    dropped = wcs.dropaxis(3)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,3]))

    wcs = WCS(naxis=4)
    wcs.wcs.cd = pc

    dropped = wcs.dropaxis(0)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([2,3,4]))
    dropped = wcs.dropaxis(1)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,3,4]))
    dropped = wcs.dropaxis(2)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,4]))
    dropped = wcs.dropaxis(3)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,3]))

def test_wcs_swapping():
    wcs = WCS(naxis=4)
    wcs.wcs.pc = np.zeros([4,4])
    np.fill_diagonal(wcs.wcs.pc, np.arange(1,5))
    pc = wcs.wcs.pc # for later use below

    swapped = wcs.swapaxes(0,1)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([2,1,3,4]))
    swapped = wcs.swapaxes(0,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([4,2,3,1]))
    swapped = wcs.swapaxes(2,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([1,2,4,3]))

    wcs = WCS(naxis=4)
    wcs.wcs.cd = pc

    swapped = wcs.swapaxes(0,1)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([2,1,3,4]))
    swapped = wcs.swapaxes(0,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([4,2,3,1]))
    swapped = wcs.swapaxes(2,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([1,2,4,3]))

@pytest.mark.parametrize('ndim',(2,3))
def test_add_stokes(ndim):
    wcs = WCS(naxis=ndim)

    for ii in range(ndim+1):
        outwcs = utils.add_stokes_axis_to_wcs(wcs,ii)
        assert outwcs.wcs.naxis == ndim+1
        assert outwcs.wcs.ctype[ii] == 'STOKES'
        assert outwcs.wcs.cname[ii] == 'STOKES'

def test_slice():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1,1]
    mywcs.wcs.cdelt = [0.1,0.1]
    mywcs.wcs.crpix = [1,1]

    slice_wcs = mywcs.slice([slice(1,None),slice(0,None)])
    assert np.all(slice_wcs.wcs.crpix == np.array([1,0]))

    slice_wcs = mywcs.slice([slice(1,None,2),slice(0,None,4)])
    assert np.all(slice_wcs.wcs.crpix == np.array([0.625, 0.25]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4,0.2]))

    slice_wcs = mywcs.slice([slice(None,None,2),slice(0,None,2)])
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2,0.2]))

def test_slice_getitem():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1,1]
    mywcs.wcs.cdelt = [0.1,0.1]
    mywcs.wcs.crpix = [1,1]

    slice_wcs = mywcs[1::2, 0::4]
    assert np.all(slice_wcs.wcs.crpix == np.array([0.625,0.25]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4,0.2]))

    mywcs.wcs.crpix = [2,2]
    slice_wcs = mywcs[1::2, 0::4]
    assert np.all(slice_wcs.wcs.crpix == np.array([0.875,0.75]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4,0.2]))

    # Default: numpy order
    slice_wcs = mywcs[1::2]
    assert np.all(slice_wcs.wcs.crpix == np.array([2,0.75]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.1,0.2]))

def test_slice_fitsorder():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1,1]
    mywcs.wcs.cdelt = [0.1,0.1]
    mywcs.wcs.crpix = [1,1]

    slice_wcs = mywcs.slice([slice(1,None),slice(0,None)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0,1]))

    slice_wcs = mywcs.slice([slice(1,None,2),slice(0,None,4)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0.25,0.625]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2,0.4]))

    slice_wcs = mywcs.slice([slice(1,None,2)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0.25,1]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2,0.1]))

def test_invalid_slice():
    mywcs = WCS(naxis=2)

    with pytest.raises(ValueError) as exc:
        mywcs[0]
    assert exc.value.args[0] == ("Cannot downsample a WCS with indexing.  Use "
                                 "wcs.sub or wcs.dropaxis if you want to remove "
                                 "axes.")

    with pytest.raises(ValueError) as exc:
        mywcs[0,::2]
    assert exc.value.args[0] == ("Cannot downsample a WCS with indexing.  Use "
                                 "wcs.sub or wcs.dropaxis if you want to remove "
                                 "axes.")

def test_axis_names():
    mywcs = WCS(naxis=4)
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN','VOPT-LSR','STOKES']

    assert mywcs.axis_type_names == ['RA','DEC','VOPT','STOKES']

    mywcs.wcs.cname = ['RA','DEC','VOPT','STOKES']

    assert mywcs.axis_type_names == ['RA','DEC','VOPT','STOKES']

def test_celestial():
    mywcs = WCS(naxis=4)
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN','VOPT','STOKES']
    cel = mywcs.celestial
    assert list(cel.wcs.ctype) == ['RA---TAN','DEC--TAN']
    assert cel.axis_type_names == ['RA','DEC']

def test_wcs_to_celestial_frame():

    from ...coordinates.builtin_frames import ICRS, FK5, FK4, Galactic
    from ...time import Time

    mywcs = WCS(naxis=2)
    with pytest.raises(ValueError) as exc:
        assert utils.wcs_to_celestial_frame(mywcs) is None
    assert exc.value.args[0] == "Could not determine celestial frame corresponding to the specified WCS object"

    mywcs.wcs.ctype = ['XOFFSET', 'YOFFSET']
    with pytest.raises(ValueError):
        assert utils.wcs_to_celestial_frame(mywcs) is None

    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    frame = utils.wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, ICRS)

    mywcs.wcs.equinox = 1987.
    frame = utils.wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, FK5)
    assert frame.equinox == Time(1987., format='jyear')

    mywcs.wcs.equinox = 1982
    frame = utils.wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, FK4)
    assert frame.equinox == Time(1982., format='byear')

    mywcs.wcs.equinox = np.nan
    mywcs.wcs.ctype = ['GLON-SIN', 'GLAT-SIN']
    frame = utils.wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, Galactic)

    mywcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    mywcs.wcs.radesys = 'ICRS'

    for equinox in [np.nan, 1987, 1982]:
        mywcs.wcs.equinox = equinox
        frame = utils.wcs_to_celestial_frame(mywcs)
        assert isinstance(frame, ICRS)

    # Flipped order
    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['DEC--TAN', 'RA---TAN']
    assert isinstance(frame, ICRS)

    # More than two dimensions
    mywcs = WCS(naxis=3)
    mywcs.wcs.ctype = ['DEC--TAN', 'VELOCITY', 'RA---TAN']
    assert isinstance(frame, ICRS)


def test_wcs_to_celestial_frame_extend():

    from ...coordinates.builtin_frames import ICRS, FK5, FK4, Galactic
    from ...time import Time

    mywcs = WCS(naxis=2)
    mywcs.wcs.ctype = ['XOFFSET', 'YOFFSET']
    with pytest.raises(ValueError):
        utils.wcs_to_celestial_frame(mywcs)

    class OffsetFrame(object):
        pass

    def identify_offset(wcs):
        if wcs.wcs.ctype[0].endswith('OFFSET') and wcs.wcs.ctype[1].endswith('OFFSET'):
            return OffsetFrame()

    from ..utils import custom_frame_mappings

    with custom_frame_mappings(identify_offset):
        frame = utils.wcs_to_celestial_frame(mywcs)
    assert isinstance(frame, OffsetFrame)

    # Check that things are back to normal after the context manager
    with pytest.raises(ValueError):
        utils.wcs_to_celestial_frame(mywcs)

def test_pixscale_nodrop():
    mywcs = WCS(naxis=2)
    mywcs.wcs.cdelt = [0.1,0.1]
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN']
    assert_almost_equal(celestial_pixel_scale(mywcs).to(u.deg).value, 0.1)

    mywcs.wcs.cdelt = [-0.1,0.1]
    assert_almost_equal(celestial_pixel_scale(mywcs).to(u.deg).value, 0.1)

def test_pixscale_withdrop():
    mywcs = WCS(naxis=3)
    mywcs.wcs.cdelt = [0.1,0.1,1]
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN','VOPT']
    assert_almost_equal(celestial_pixel_scale(mywcs).to(u.deg).value, 0.1)

    mywcs.wcs.cdelt = [-0.1,0.1,1]
    assert_almost_equal(celestial_pixel_scale(mywcs).to(u.deg).value, 0.1)


def test_pixscale_cd():
    mywcs = WCS(naxis=2)
    mywcs.wcs.cd = [[-0.1,0],[0,0.1]]
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN']
    assert_almost_equal(celestial_pixel_scale(mywcs).to(u.deg).value, 0.1)

def test_pixscale_warning(recwarn):
    mywcs = WCS(naxis=2)
    mywcs.wcs.cd = [[-0.1,0],[0,0.1]]
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN']

    with catch_warnings(AstropyUserWarning) as warning_lines:

        celestial_pixel_scale(mywcs)
        assert ("Pixel sizes may very over the image for "
                "projection class TAN"
                in str(warning_lines[0].message))

def test_pixscale_asymmetric():
    mywcs = WCS(naxis=2)
    mywcs.wcs.cd = [[-0.2,0],[0,0.1]]
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN']

    with pytest.raises(ValueError) as exc:
        celestial_pixel_scale(mywcs)
    assert exc.value.args[0] == "Pixels are not square: 'pixel scale' is ambiguous"

@pytest.mark.parametrize('angle',
                         (30,45,60,75))
def test_pixscale_cd_rotated(angle):
    mywcs = WCS(naxis=2)
    rho = angle/180.*np.pi
    scale = 0.1
    mywcs.wcs.cd = [[scale*np.cos(rho), -scale*np.sin(rho)],
                    [scale*np.sin(rho), scale*np.cos(rho)]]
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN']
    assert_almost_equal(celestial_pixel_scale(mywcs).to(u.deg).value, 0.1)

@pytest.mark.parametrize('angle',
                         (30,45,60,75))
def test_pixscale_pc_rotated(angle):
    mywcs = WCS(naxis=2)
    rho = angle/180.*np.pi
    scale = 0.1
    mywcs.wcs.cdelt = [-scale, scale]
    mywcs.wcs.pc = [[np.cos(rho), -np.sin(rho)],
                    [np.sin(rho), np.cos(rho)]]
    mywcs.wcs.ctype = ['RA---TAN','DEC--TAN']
    assert_almost_equal(celestial_pixel_scale(mywcs).to(u.deg).value, 0.1)

@pytest.mark.parametrize(('cdelt','pc','pccd'),
                         (([0.1,0.2], np.eye(2), np.diag([0.1,0.2])),
                          ([0.1,0.2,0.3], np.eye(3), np.diag([0.1,0.2,0.3])),
                          ([1,1,1], np.diag([0.1,0.2,0.3]), np.diag([0.1,0.2,0.3])),
                         ))
def test_pixel_scale_matrix(cdelt, pc, pccd):

    mywcs = WCS(naxis=(len(cdelt)))
    mywcs.wcs.cdelt = cdelt
    mywcs.wcs.pc = pc

    assert_almost_equal(mywcs.pixel_scale_matrix, pccd)

@pytest.mark.parametrize(('ctype', 'cel'),
                         ((['RA---TAN','DEC--TAN'], True),
                          (['RA---TAN','DEC--TAN','FREQ'], False),
                          (['RA---TAN','FREQ'], False),))
def test_is_celestial(ctype, cel):
    mywcs = WCS(naxis=len(ctype))
    mywcs.wcs.ctype = ctype

    assert mywcs.is_celestial == cel

@pytest.mark.parametrize(('ctype', 'cel'),
                         ((['RA---TAN','DEC--TAN'], True),
                          (['RA---TAN','DEC--TAN','FREQ'], True),
                          (['RA---TAN','FREQ'], False),))
def test_has_celestial(ctype, cel):
    mywcs = WCS(naxis=len(ctype))
    mywcs.wcs.ctype = ctype

    assert mywcs.has_celestial == cel

@pytest.mark.parametrize(('cdelt','pc','cd'),
                         ((np.array([0.1,0.2]), np.eye(2), np.eye(2)),
                          (np.array([1,1]), np.diag([0.1,0.2]), np.eye(2)),
                          (np.array([0.1,0.2]), np.eye(2), None),
                          (np.array([0.1,0.2]), None, np.eye(2)),
                          )
                        )
def test_noncelestial_scale(cdelt, pc, cd):
    mywcs = WCS(naxis=2)
    if cd is not None:
        mywcs.wcs.cd = cd
    if pc is not None:
        mywcs.wcs.pc = pc
    mywcs.wcs.cdelt = cdelt

    mywcs.wcs.ctype = ['RA---TAN','FREQ']

    ps = non_celestial_pixel_scales(mywcs)

    assert_almost_equal(ps.to(u.deg).value, np.array([0.1,0.2]))

@pytest.mark.parametrize('mode', ['all', 'wcs'])
def test_skycoord_to_pixel(mode):

    from ... import units as u
    from ...coordinates import SkyCoord
    from ..utils import skycoord_to_pixel, pixel_to_skycoord

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


@pytest.mark.parametrize('mode', ['all', 'wcs'])
def test_skycoord_to_pixel_distortions(mode):

    from ... import units as u
    from ...coordinates import SkyCoord
    from ..utils import skycoord_to_pixel, pixel_to_skycoord

    header = get_pkg_data_filename('data/sip.fits')
    wcs = WCS(header)

    ref = SkyCoord(202.50 * u.deg, 47.19 * u.deg, frame='icrs')

    xp, yp = skycoord_to_pixel(ref, wcs, mode=mode)

    # WCS is in FK5 so we need to transform back to ICRS
    new = pixel_to_skycoord(xp, yp, wcs, mode=mode).transform_to('icrs')

    assert_allclose(new.ra.degree, ref.ra.degree)
    assert_allclose(new.dec.degree, ref.dec.degree)
