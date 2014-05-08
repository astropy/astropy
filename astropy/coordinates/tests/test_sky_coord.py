# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for the SkyCoord class.  Note that there are also SkyCoord tests in
test_api_ape5.py
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import functools

import numpy as np

from ... import units as u
from ...tests.helper import pytest
from ...coordinates import (ICRS, FK4, FK5, FK4NoETerms, Galactic, SkyCoord, Angle,
                            SphericalRepresentation, CartesianRepresentation)
from ...time import Time

RA = 1.0 * u.deg
DEC = 2.0 * u.deg
C_ICRS = ICRS(RA, DEC)
C_FK5 = C_ICRS.transform_to(FK5)
J2001 = Time('J2001', scale='utc')

allclose = functools.partial(np.allclose, rtol=0.0, atol=1e-8)


def test_transform_to():
    for frame in (FK5, FK5(equinox=Time('J1975.0')),
                  FK4, FK4(equinox=Time('J1975.0')),
                  SkyCoord(RA, DEC, 'fk4', equinox='J1980')):
        c_frame = C_ICRS.transform_to(frame)
        s_icrs = SkyCoord(RA, DEC, frame='icrs')
        s_frame = s_icrs.transform_to(frame)
        assert allclose(c_frame.ra, s_frame.ra)
        assert allclose(c_frame.dec, s_frame.dec)
        assert allclose(c_frame.distance, s_frame.distance)


# set up for parametrized test
rt_sets = []
rt_frames = [ICRS, FK4, FK5, FK4NoETerms, Galactic]
for rt_frame0 in rt_frames:
    for rt_frame1 in rt_frames:
            for equinox0 in (None, 'J1975.0'):
                for obstime0 in (None, 'J1980.0'):
                    for equinox1 in (None, 'J1975.0'):
                        for obstime1 in (None, 'J1980.0'):
                            rt_sets.append([rt_frame0, rt_frame1,
                                            equinox0, equinox1,
                                            obstime0, obstime1])
rt_args = 'frame0,frame1,equinox0,equinox1,obstime0,obstime1'
@pytest.mark.parametrize(rt_args, rt_sets)
def test_round_tripping(frame0, frame1, equinox0, equinox1, obstime0, obstime1):
    """
    Test round tripping out and back using transform_to in every combination.
    """
    attrs0 = {'equinox': equinox0, 'obstime': obstime0}
    attrs1 = {'equinox': equinox1, 'obstime': obstime1}

    # Remove None values
    attrs0 = dict((k, v) for k, v in attrs0.items() if v is not None)
    attrs1 = dict((k, v) for k, v in attrs1.items() if v is not None)

    # Go out and back
    sc = SkyCoord(frame0, RA, DEC, **attrs0)

    # Keep only frame attributes for frame1
    attrs1 = dict((attr, val) for attr, val in attrs1.items()
                  if attr in frame1.frame_attr_names)
    sc2 = sc.transform_to(frame1(**attrs1))

    # When coming back only keep frame0 attributes for transform_to
    attrs0 = dict((attr, val) for attr, val in attrs0.items()
                  if attr in frame0.frame_attr_names)
    # also, if any are None, fill in with defaults
    for attrnm in frame0.frame_attr_names:
        if attrs0.get(attrnm, None) is None:
            if attrnm=='obstime' and frame0.frame_attr_names[attrnm] is None:
                attrs0[attrnm] = attrs0['equinox']
            else:
                attrs0[attrnm] = frame0.frame_attr_names[attrnm]
    sc_rt = sc2.transform_to(frame0(**attrs0))

    if frame0 is Galactic:
        assert allclose(sc.l, sc_rt.l)
        assert allclose(sc.b, sc_rt.b)
    else:
        assert allclose(sc.ra, sc_rt.ra)
        assert allclose(sc.dec, sc_rt.dec)
    if equinox0:
        assert Time(sc.equinox) == Time(sc_rt.equinox)
    if obstime0:
        assert Time(sc.obstime) == Time(sc_rt.obstime)




def test_coord_init_string():
    """
    Spherical or Cartesian represenation input coordinates.
    """
    sc = SkyCoord('1d 2d')
    assert allclose(sc.ra, 1 * u.deg)
    assert allclose(sc.dec, 2 * u.deg)

    sc = SkyCoord('1d', '2d')
    assert allclose(sc.ra, 1 * u.deg)
    assert allclose(sc.dec, 2 * u.deg)

    sc = SkyCoord(u'1°2′3″', u'2°3′4″')
    assert allclose(sc.ra, Angle(u'1°2′3″'))
    assert allclose(sc.dec, Angle(u'2°3′4″'))

    sc = SkyCoord(u'1°2′3″ 2°3′4″')
    assert allclose(sc.ra, Angle(u'1°2′3″'))
    assert allclose(sc.dec, Angle(u'2°3′4″'))

    with pytest.raises(ValueError) as err:
        SkyCoord('1d 2d 3d')
    assert "Cannot parse longitude and latitude" in str(err)


def test_coord_init_list():
    """
    Spherical or Cartesian representation input coordinates.
    """
    sc = SkyCoord([('1d', '2d'),
                   (1 * u.deg, 2 * u.deg),
                   '1d 2d',
                   (u'1°', u'2°'),
                   u'1° 2°'], unit='deg')
    assert allclose(sc.ra, Angle('1d'))
    assert allclose(sc.dec, Angle('2d'))

    with pytest.raises(ValueError) as err:
        SkyCoord(['1d 2d 3d'])
    assert "Cannot parse longitude and latitude" in str(err)

    with pytest.raises(ValueError) as err:
        SkyCoord([('1d', '2d', '3d')])
    assert "Cannot parse longitude and latitude" in str(err)

    sc = SkyCoord([1 * u.deg, 1 * u.deg], [2 * u.deg, 2 * u.deg])
    assert allclose(sc.ra, Angle('1d'))
    assert allclose(sc.dec, Angle('2d'))

    with pytest.raises(ValueError) as err:
        SkyCoord([1 * u.deg, 2 * u.deg])  # this list is taken as RA w/ missing dec
    assert "Cannot parse longitude and latitude" in str(err)


def test_coord_init_representation():
    """
    Spherical or Cartesian represenation input coordinates.
    """
    coord = SphericalRepresentation(lon=8 * u.deg, lat=5 * u.deg, distance=1 * u.kpc)
    sc = SkyCoord(coord, 'icrs')
    assert allclose(sc.ra, coord.lon)
    assert allclose(sc.dec, coord.lat)
    assert allclose(sc.distance, coord.distance)

    with pytest.raises(ValueError) as err:
        SkyCoord(coord, 'icrs', ra='1d')
    assert "conflicts with keyword argument 'ra'" in str(err)

    coord = CartesianRepresentation(1, 2, 3)
    sc = SkyCoord(coord, 'icrs')
    sc_cart = sc.represent_as(CartesianRepresentation)
    assert allclose(sc_cart.x, 1.0)
    assert allclose(sc_cart.y, 2.0)
    assert allclose(sc_cart.z, 3.0)


def test_frame_init():
    """
    Different ways of providing the frame.
    """
    sc = SkyCoord(RA, DEC, frame='icrs')
    assert sc.frame == 'icrs'

    sc = SkyCoord(RA, DEC, frame=ICRS)
    assert sc.frame == 'icrs'

    sc = SkyCoord(RA, DEC, 'icrs')
    assert sc.frame == 'icrs'

    sc = SkyCoord(RA, DEC, ICRS)
    assert sc.frame == 'icrs'

    sc = SkyCoord('icrs', RA, DEC)
    assert sc.frame == 'icrs'

    sc = SkyCoord(ICRS, RA, DEC)
    assert sc.frame == 'icrs'

    sc = SkyCoord(sc)
    assert sc.frame == 'icrs'

    sc = SkyCoord(C_ICRS)
    assert sc.frame == 'icrs'

    SkyCoord(C_ICRS, frame='icrs')
    assert sc.frame == 'icrs'

    with pytest.raises(ValueError) as err:
        SkyCoord(C_ICRS, frame='galactic')
    assert 'Cannot override frame=' in str(err)


def test_attr_inheritance():
    """
    When initializing from an existing coord the preferred attrs like
    equinox should be inherited to the SkyCoord.  If there is a conflict
    then raise an exception.
    """
    sc = SkyCoord('icrs', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    sc2 = SkyCoord(sc)
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc2 = SkyCoord(sc._coord)  # Doesn't have equinox there so we get FK4 defaults
    assert sc2.equinox != sc.equinox
    assert sc2.obstime != sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc = SkyCoord('fk4', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    sc2 = SkyCoord(sc)
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc2 = SkyCoord(sc._coord)  # sc._coord has equinox, obstime
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)


def test_attr_conflicts():
    """
    Check conflicts resolution between coordinate attributes and init kwargs.
    """
    sc = SkyCoord('icrs', 1, 2, unit='deg', equinox='J1999', obstime='J2001')

    # OK if attrs both specified but with identical values
    SkyCoord(sc, equinox='J1999', obstime='J2001')

    # OK because sc._coord doesn't have obstime
    SkyCoord(sc._coord, equinox='J1999', obstime='J2100')

    # Not OK if attrs don't match
    with pytest.raises(ValueError) as err:
        SkyCoord(sc, equinox='J1999', obstime='J2002')
    assert "Coordinate attribute 'obstime'=" in str(err)

    # Same game but with fk4 which has equinox and obstime frame attrs
    sc = SkyCoord('fk4', 1, 2, unit='deg', equinox='J1999', obstime='J2001')

    # OK if attrs both specified but with identical values
    SkyCoord(sc, equinox='J1999', obstime='J2001')

    # Not OK if SkyCoord attrs don't match
    with pytest.raises(ValueError) as err:
        SkyCoord(sc, equinox='J1999', obstime='J2002')
    assert "Coordinate attribute 'obstime'=" in str(err)

    # Not OK because sc._coord has different attrs
    with pytest.raises(ValueError) as err:
        SkyCoord(sc._coord, equinox='J1999', obstime='J2002')
    assert "Coordinate attribute 'obstime'=" in str(err)


def test_frame_attr_getattr():
    """
    When accessing frame attributes like equinox, the value should come
    from self._coord when that object has the relevant attribute, otherwise
    from self.
    """
    sc = SkyCoord('icrs', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    assert sc.equinox == 'J1999'  # Just the raw value (not validated)
    assert sc.obstime == 'J2001'

    sc = SkyCoord('fk4', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    assert sc.equinox == Time('J1999')  # Coming from the self._coord object
    assert sc.obstime == Time('J2001')

    sc = SkyCoord('fk4', 1, 2, unit='deg', equinox='J1999')
    assert sc.equinox == Time('J1999')
    assert sc.obstime == Time('J1999')


def test_to_string():
    """
    Basic testing of converting SkyCoord to strings.  This just tests
    for a single input coordinate and and 1-element list.  It does not
    test the underlying `Angle.to_string` method itself.
    """
    coord = '1h2m3s 1d2m3s'
    for wrap in (lambda x: x, lambda x: [x]):
        sc = SkyCoord(wrap(coord))
        assert sc.to_string() == wrap('15.5125 1.03417')
        assert sc.to_string('dms') == wrap('15d30m45s 1d02m03s')
        assert sc.to_string('hmsdms') == wrap('01h02m03s +01d02m03s')
        with_kwargs = sc.to_string('hmsdms', precision=3, pad=True, alwayssign=True)
        assert with_kwargs == wrap('+01h02m03.000s +01d02m03.000s')


def test_seps():
    sc1 = SkyCoord('icrs', 0*u.deg, 1*u.deg)
    sc2 = SkyCoord('icrs', 0*u.deg, 2*u.deg)

    sep = sc1.separation(sc2)

    assert (sep - 1*u.deg)/u.deg < 1e-10

    with pytest.raises(ValueError):
        sc1.separation_3d(sc2)

    sc3 = SkyCoord('icrs', 1*u.deg, 1*u.deg, distance=1*u.kpc)
    sc4 = SkyCoord('icrs', 1*u.deg, 1*u.deg, distance=2*u.kpc)
    sep3d = sc3.separation_3d(sc4)

    assert sep3d == 1*u.kpc


def test_repr():
    sc1 = SkyCoord('icrs', 0*u.deg, 1*u.deg)
    sc2 = SkyCoord('icrs', 1*u.deg, 1*u.deg, distance=1*u.kpc)

    assert repr(sc1) == '<SkyCoord (ICRS): ra=0.0 deg, dec=1.0 deg>'
    assert repr(sc2) == '<SkyCoord (ICRS): ra=1.0 deg, dec=1.0 deg, distance=1.0 kpc>'

    sc3 = SkyCoord('icrs', 0.1*u.deg, [1, 2.5]*u.deg)
    assert repr(sc3) == ('<SkyCoord (ICRS): (ra, dec) in deg\n'
                         '    [(0.10000000000002274, 1.0), (0.10000000000002274, 2.5)]>')


def test_ops():
    """
    Tests miscellaneous operations like `len`
    """
    sc = SkyCoord('icrs', 0*u.deg, 1*u.deg)
    sc_arr = SkyCoord('icrs', 0*u.deg, [1, 2]*u.deg)
    sc_empty = SkyCoord('icrs', []*u.deg, []*u.deg)

    assert sc.isscalar
    assert not sc_arr.isscalar
    assert not sc_empty.isscalar

    with pytest.raises(TypeError):
        len(sc)
    assert len(sc_arr) == 2
    assert len(sc_empty) == 0

    assert bool(sc)
    assert bool(sc_arr)
    assert not bool(sc_empty)


