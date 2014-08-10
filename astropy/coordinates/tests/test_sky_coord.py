# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for the SkyCoord class.  Note that there are also SkyCoord tests in
test_api_ape5.py
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy
import functools

import numpy as np
from numpy import testing as npt

from ... import units as u
from ...tests.helper import pytest, catch_warnings
from ..representation import REPRESENTATION_CLASSES
from ...coordinates import (ICRS, FK4, FK5, Galactic, SkyCoord, Angle,
                            SphericalRepresentation, CartesianRepresentation)
from ...coordinates import Latitude, Longitude
from ...time import Time
from ...utils.exceptions import AstropyDeprecationWarning

RA = 1.0 * u.deg
DEC = 2.0 * u.deg
C_ICRS = ICRS(RA, DEC)
C_FK5 = C_ICRS.transform_to(FK5)
J2001 = Time('J2001', scale='utc')

allclose = functools.partial(np.allclose, rtol=0.0, atol=1e-8)

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


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
rt_frames = [ICRS, FK4, FK5, Galactic]
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
                  if attr in frame1.get_frame_attr_names())
    sc2 = sc.transform_to(frame1(**attrs1))

    # When coming back only keep frame0 attributes for transform_to
    attrs0 = dict((attr, val) for attr, val in attrs0.items()
                  if attr in frame0.get_frame_attr_names())
    # also, if any are None, fill in with defaults
    for attrnm in frame0.get_frame_attr_names():
        if attrs0.get(attrnm, None) is None:
            if attrnm == 'obstime' and frame0.get_frame_attr_names()[attrnm] is None:
                if 'equinox' in attrs0:
                    attrs0[attrnm] = attrs0['equinox']
            else:
                attrs0[attrnm] = frame0.get_frame_attr_names()[attrnm]
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

    sc = SkyCoord('1°2′3″', '2°3′4″')
    assert allclose(sc.ra, Angle('1°2′3″'))
    assert allclose(sc.dec, Angle('2°3′4″'))

    sc = SkyCoord('1°2′3″ 2°3′4″')
    assert allclose(sc.ra, Angle('1°2′3″'))
    assert allclose(sc.dec, Angle('2°3′4″'))

    with pytest.raises(ValueError) as err:
        SkyCoord('1d 2d 3d')
    assert "Cannot parse longitude and latitude" in str(err)

    sc1 = SkyCoord('8 00 00 +5 00 00.0', unit=(u.hour, u.deg), frame='icrs')
    assert isinstance(sc1, SkyCoord)
    assert allclose(sc1.ra, Angle(120 * u.deg))
    assert allclose(sc1.dec, Angle(5 * u.deg))

    sc2 = SkyCoord('8 00 -5 00 00.0', unit=(u.hour, u.deg), frame='icrs')
    assert isinstance(sc2, SkyCoord)
    assert allclose(sc2.ra, Angle(120 * u.deg))
    assert allclose(sc2.dec, Angle(-5 * u.deg))

    sc3 = SkyCoord('8 00 -5 00.6', unit=(u.hour, u.deg), frame='icrs')
    assert isinstance(sc3, SkyCoord)
    assert allclose(sc3.ra, Angle(120 * u.deg))
    assert allclose(sc3.dec, Angle(-5.01 * u.deg))

    sc5 = SkyCoord('8h00.6m -5d00.6m', unit=(u.hour, u.deg), frame='icrs')
    assert isinstance(sc5, SkyCoord)
    assert allclose(sc5.ra, Angle(120.15 * u.deg))
    assert allclose(sc5.dec, Angle(-5.01 * u.deg))


def test_coord_init_unit():
    """
    Test variations of the unit keyword.
    """
    for unit in ('deg', 'deg,deg', ' deg , deg ', u.deg, (u.deg, u.deg),
                 np.array(['deg', 'deg'])):
        sc = SkyCoord(1, 2, unit=unit)
        assert allclose(sc.ra, Angle(1 * u.deg))
        assert allclose(sc.dec, Angle(2 * u.deg))

    for unit in ('hourangle', 'hourangle,hourangle', ' hourangle , hourangle ',
                 u.hourangle, [u.hourangle, u.hourangle]):
        sc = SkyCoord(1, 2, unit=unit)
        assert allclose(sc.ra, Angle(15 * u.deg))
        assert allclose(sc.dec, Angle(30 * u.deg))

    for unit in ('hourangle,deg', (u.hourangle, u.deg)):
        sc = SkyCoord(1, 2, unit=unit)
        assert allclose(sc.ra, Angle(15 * u.deg))
        assert allclose(sc.dec, Angle(2 * u.deg))

    for unit in ('deg,deg,deg,deg', [u.deg, u.deg, u.deg, u.deg], None):
        with pytest.raises(ValueError) as err:
            SkyCoord(1, 2, unit=unit)
        assert 'Unit keyword must have one to three unit values' in str(err)

    for unit in ('m', (u.m, u.deg), ''):
        with pytest.raises(u.UnitsError) as err:
            SkyCoord(1, 2, unit=unit)


def test_coord_init_list():
    """
    Spherical or Cartesian representation input coordinates.
    """
    sc = SkyCoord([('1d', '2d'),
                   (1 * u.deg, 2 * u.deg),
                   '1d 2d',
                   ('1°', '2°'),
                   '1° 2°'], unit='deg')
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
    assert "One or more elements of input sequence does not have a length" in str(err)


def test_coord_init_array():
    """
    Input in the form of a list array or numpy array
    """
    for a in (['1 2', '3 4'],
              [['1', '2'], ['3', '4']],
              [[1, 2], [3, 4]]):
        sc = SkyCoord(a, unit='deg')
        assert allclose(sc.ra - [1, 3] * u.deg, 0)
        assert allclose(sc.dec - [2, 4] * u.deg, 0)

        sc = SkyCoord(np.array(a), unit='deg')
        assert allclose(sc.ra - [1, 3] * u.deg, 0)
        assert allclose(sc.dec - [2, 4] * u.deg, 0)


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

    coord = CartesianRepresentation(1 * u.one, 2 * u.one, 3 * u.one)
    sc = SkyCoord(coord, 'icrs')
    sc_cart = sc.represent_as(CartesianRepresentation)
    assert allclose(sc_cart.x, 1.0)
    assert allclose(sc_cart.y, 2.0)
    assert allclose(sc_cart.z, 3.0)


FRAME_DEPRECATION_WARNING = ("Passing a frame as a positional argument is now "
                             "deprecated, use the frame= keyword argument "
                             "instead.")

def test_frame_init():
    """
    Different ways of providing the frame.
    """

    sc = SkyCoord(RA, DEC, frame='icrs')
    assert sc.frame.name == 'icrs'

    sc = SkyCoord(RA, DEC, frame=ICRS)
    assert sc.frame.name == 'icrs'

    with catch_warnings(AstropyDeprecationWarning) as w:
        sc = SkyCoord(RA, DEC, 'icrs')
    assert sc.frame.name == 'icrs'
    assert len(w) == 1
    assert str(w[0].message) == FRAME_DEPRECATION_WARNING

    with catch_warnings(AstropyDeprecationWarning) as w:
        sc = SkyCoord(RA, DEC, ICRS)
    assert sc.frame.name == 'icrs'
    assert len(w) == 1
    assert str(w[0].message) == FRAME_DEPRECATION_WARNING

    with catch_warnings(AstropyDeprecationWarning) as w:
        sc = SkyCoord('icrs', RA, DEC)
    assert sc.frame.name == 'icrs'
    assert len(w) == 1
    assert str(w[0].message) == FRAME_DEPRECATION_WARNING

    with catch_warnings(AstropyDeprecationWarning) as w:
        sc = SkyCoord(ICRS, RA, DEC)
    assert sc.frame.name == 'icrs'
    assert len(w) == 1
    assert str(w[0].message) == FRAME_DEPRECATION_WARNING

    sc = SkyCoord(sc)
    assert sc.frame.name == 'icrs'

    sc = SkyCoord(C_ICRS)
    assert sc.frame.name == 'icrs'

    SkyCoord(C_ICRS, frame='icrs')
    assert sc.frame.name == 'icrs'

    with pytest.raises(ValueError) as err:
        SkyCoord(C_ICRS, frame='galactic')
    assert 'Cannot override frame=' in str(err)


def test_attr_inheritance():
    """
    When initializing from an existing coord the representation attrs like
    equinox should be inherited to the SkyCoord.  If there is a conflict
    then raise an exception.
    """
    sc = SkyCoord(1, 2, frame='icrs', unit='deg', equinox='J1999', obstime='J2001')
    sc2 = SkyCoord(sc)
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc2 = SkyCoord(sc.frame)  # Doesn't have equinox there so we get FK4 defaults
    assert sc2.equinox != sc.equinox
    assert sc2.obstime != sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc = SkyCoord(1, 2, frame='fk4', unit='deg', equinox='J1999', obstime='J2001')
    sc2 = SkyCoord(sc)
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)

    sc2 = SkyCoord(sc.frame)  # sc.frame has equinox, obstime
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime
    assert allclose(sc2.ra, sc.ra)
    assert allclose(sc2.dec, sc.dec)
    assert allclose(sc2.distance, sc.distance)


def test_attr_conflicts():
    """
    Check conflicts resolution between coordinate attributes and init kwargs.
    """
    sc = SkyCoord(1, 2, frame='icrs', unit='deg', equinox='J1999', obstime='J2001')

    # OK if attrs both specified but with identical values
    SkyCoord(sc, equinox='J1999', obstime='J2001')

    # OK because sc.frame doesn't have obstime
    SkyCoord(sc.frame, equinox='J1999', obstime='J2100')

    # Not OK if attrs don't match
    with pytest.raises(ValueError) as err:
        SkyCoord(sc, equinox='J1999', obstime='J2002')
    assert "Coordinate attribute 'obstime'=" in str(err)

    # Same game but with fk4 which has equinox and obstime frame attrs
    sc = SkyCoord(1, 2, frame='fk4', unit='deg', equinox='J1999', obstime='J2001')

    # OK if attrs both specified but with identical values
    SkyCoord(sc, equinox='J1999', obstime='J2001')

    # Not OK if SkyCoord attrs don't match
    with pytest.raises(ValueError) as err:
        SkyCoord(sc, equinox='J1999', obstime='J2002')
    assert "Coordinate attribute 'obstime'=" in str(err)

    # Not OK because sc.frame has different attrs
    with pytest.raises(ValueError) as err:
        SkyCoord(sc.frame, equinox='J1999', obstime='J2002')
    assert "Coordinate attribute 'obstime'=" in str(err)


def test_frame_attr_getattr():
    """
    When accessing frame attributes like equinox, the value should come
    from self.frame when that object has the relevant attribute, otherwise
    from self.
    """
    sc = SkyCoord(1, 2, frame='icrs', unit='deg', equinox='J1999', obstime='J2001')
    assert sc.equinox == 'J1999'  # Just the raw value (not validated)
    assert sc.obstime == 'J2001'

    sc = SkyCoord(1, 2, frame='fk4', unit='deg', equinox='J1999', obstime='J2001')
    assert sc.equinox == Time('J1999')  # Coming from the self.frame object
    assert sc.obstime == Time('J2001')

    sc = SkyCoord(1, 2, frame='fk4', unit='deg', equinox='J1999')
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
    sc1 = SkyCoord(0 * u.deg, 1 * u.deg, frame='icrs')
    sc2 = SkyCoord(0 * u.deg, 2 * u.deg, frame='icrs')

    sep = sc1.separation(sc2)

    assert (sep - 1 * u.deg)/u.deg < 1e-10

    with pytest.raises(ValueError):
        sc1.separation_3d(sc2)

    sc3 = SkyCoord(1 * u.deg, 1 * u.deg, distance=1 * u.kpc, frame='icrs')
    sc4 = SkyCoord(1 * u.deg, 1 * u.deg, distance=2 * u.kpc, frame='icrs')
    sep3d = sc3.separation_3d(sc4)

    assert sep3d == 1 * u.kpc


def test_repr():
    # Repr tests must use exact floating point vals because Python 2.6
    # outputs values like 0.1 as 0.1000000000001.  No workaround found.
    sc1 = SkyCoord(0 * u.deg, 1 * u.deg, frame='icrs')
    sc2 = SkyCoord(1 * u.deg, 1 * u.deg, frame='icrs', distance=1 * u.kpc)

    assert repr(sc1) == '<SkyCoord (ICRS): ra=0.0 deg, dec=1.0 deg>'
    assert repr(sc2) == '<SkyCoord (ICRS): ra=1.0 deg, dec=1.0 deg, distance=1.0 kpc>'

    sc3 = SkyCoord(0.25 * u.deg, [1, 2.5] * u.deg, frame='icrs')
    assert repr(sc3) == ('<SkyCoord (ICRS): (ra, dec) in deg\n'
                         '    [(0.25, 1.0), (0.25, 2.5)]>')

    sc_default = SkyCoord(0 * u.deg, 1 * u.deg)
    assert repr(sc_default) == '<SkyCoord (ICRS): ra=0.0 deg, dec=1.0 deg>'


def test_ops():
    """
    Tests miscellaneous operations like `len`
    """
    sc = SkyCoord(0 * u.deg, 1 * u.deg, frame='icrs')
    sc_arr = SkyCoord(0 * u.deg, [1, 2] * u.deg, frame='icrs')
    sc_empty = SkyCoord([] * u.deg, [] * u.deg, frame='icrs')

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

    assert sc_arr[0].isscalar
    assert len(sc_arr[:1]) == 1
    with pytest.raises(TypeError):
        assert sc[0:]  # scalar, so it shouldn't be indexable


def test_none_transform():
    """
    Ensure that transforming from a SkyCoord with no frame provided works like
    ICRS
    """
    sc = SkyCoord(0 * u.deg, 1 * u.deg)
    sc_arr = SkyCoord(0 * u.deg, [1, 2] * u.deg)

    sc2 = sc.transform_to(ICRS)
    assert sc.ra == sc2.ra and sc.dec == sc2.dec

    sc5 = sc.transform_to('fk5')
    assert sc5.ra == sc2.transform_to('fk5').ra

    sc_arr2 = sc_arr.transform_to(ICRS)
    sc_arr5 = sc_arr.transform_to('fk5')
    npt.assert_array_equal(sc_arr5.ra, sc_arr2.transform_to('fk5').ra)


def test_position_angle():
    c1 = SkyCoord(0*u.deg, 0*u.deg)

    c2 = SkyCoord(1*u.deg, 0*u.deg)
    npt.assert_allclose(c1.position_angle(c2) - 90.0 * u.deg, 0)

    c3 = SkyCoord(1*u.deg, 0.1*u.deg)
    assert c1.position_angle(c3) < 90*u.deg

    c4 = SkyCoord(0*u.deg, 1*u.deg)
    npt.assert_allclose(c1.position_angle(c4), 0)

    carr1 = SkyCoord(0*u.deg, [0, 1, 2]*u.deg)
    carr2 = SkyCoord([-1, -2, -3]*u.deg, [0.1, 1.1, 2.1]*u.deg)

    res = carr1.position_angle(carr2)
    assert res.shape == (3,)
    assert np.all(res < 360*u.degree)
    assert np.all(res > 270*u.degree)

    cicrs = SkyCoord(0*u.deg, 0*u.deg, frame='icrs')
    cfk5 = SkyCoord(1*u.deg, 0*u.deg, frame='fk5')
    # because of the frame transform, it's just a *bit* more than 90 degrees
    assert cicrs.position_angle(cfk5) > 90.0 * u.deg
    assert cicrs.position_angle(cfk5) < 91.0 * u.deg


def test_table_to_coord():
    """
    Checks "end-to-end" use of `Table` with `SkyCoord` - the `Quantity`
    initializer is the intermediary that translate the table columns into
    something coordinates understands.

    (Regression test for #1762 )
    """
    from ...table import Table, Column

    t = Table()
    t.add_column(Column(data=[1, 2, 3], name='ra', unit=u.deg))
    t.add_column(Column(data=[4, 5, 6], name='dec', unit=u.deg))

    c = SkyCoord(t['ra'], t['dec'])

    assert allclose(c.ra.to(u.deg), [1, 2, 3])
    assert allclose(c.dec.to(u.deg), [4, 5, 6])


def assert_quantities_allclose(coord, q1s, attrs):
    """
    Compare two tuples of quantities.  This assumes that the values in q1 are of
    order(1) and uses atol=1e-13, rtol=0.  It also asserts that the units of the
    two quantities are the *same*, in order to check that the representation
    output has the expected units.
    """
    q2s = [getattr(coord, attr) for attr in attrs]
    assert len(q1s) == len(q2s)
    for q1, q2 in zip(q1s, q2s):
        assert q1.shape == q2.shape
        dq = q1 - q2
        assert np.allclose(dq.value, 0.0, rtol=0, atol=1e-13)


# Sets of inputs corresponding to Galactic frame
base_unit_attr_sets = [
    ('spherical', u.karcsec, u.karcsec, u.kpc, Latitude, 'l', 'b', 'distance'),
    ('unitspherical', u.karcsec, u.karcsec, None, Latitude, 'l', 'b', None),
    ('physicsspherical', u.karcsec, u.karcsec, u.kpc, Angle, 'phi', 'theta', 'r'),
    ('cartesian', u.km, u.km, u.km, u.Quantity, 'w', 'u', 'v'),
    ('cylindrical', u.km, u.karcsec, u.km, Angle, 'rho', 'phi', 'z')
]

units_attr_sets = []
for base_unit_attr_set in base_unit_attr_sets:
    repr_name = base_unit_attr_set[0]
    for representation in (repr_name, REPRESENTATION_CLASSES[repr_name]):
        for c1, c2, c3 in ((1, 2, 3), ([1], [2], [3])):
            for arrayify in True, False:
                if arrayify:
                    c1 = np.array(c1)
                    c2 = np.array(c2)
                    c3 = np.array(c3)
                units_attr_sets.append(base_unit_attr_set + (representation, c1, c2, c3))
units_attr_args = 'repr_name,unit1,unit2,unit3,cls2,attr1,attr2,attr3,representation,c1,c2,c3'


@pytest.mark.parametrize(units_attr_args,
                         (x for x in units_attr_sets if x[0] != 'unitspherical'))
def test_skycoord_three_components(repr_name, unit1, unit2, unit3, cls2, attr1, attr2, attr3,
                                   representation, c1, c2, c3):
    """
    Tests positional inputs using components (COMP1, COMP2, COMP3)
    and various representations.  Use weird units and Galactic frame.
    """
    sc = SkyCoord(Galactic, c1, c2, c3, unit=(unit1, unit2, unit3),
                  representation=representation)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2, c3*unit3),
                               (attr1, attr2, attr3))

    sc = SkyCoord(1000*c1*u.Unit(unit1/1000), cls2(c2, unit=unit2),
                  1000*c3*u.Unit(unit3/1000), Galactic,
                  unit=(unit1, unit2, unit3), representation=representation)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2, c3*unit3),
                               (attr1, attr2, attr3))

    kwargs = {attr3: c3}
    sc = SkyCoord(Galactic, c1, c2, unit=(unit1, unit2, unit3),
                  representation=representation, **kwargs)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2, c3*unit3),
                               (attr1, attr2, attr3))

    kwargs = {attr1: c1, attr2: c2, attr3: c3}
    sc = SkyCoord(Galactic, unit=(unit1, unit2, unit3),
                  representation=representation, **kwargs)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2, c3*unit3),
                               (attr1, attr2, attr3))


@pytest.mark.parametrize(units_attr_args,
                         (x for x in units_attr_sets
                          if x[0] in ('spherical', 'unitspherical')))
def test_skycoord_spherical_two_components(repr_name, unit1, unit2, unit3, cls2,
                                           attr1, attr2, attr3, representation, c1, c2, c3):
    """
    Tests positional inputs using components (COMP1, COMP2) for spherical
    representations.  Use weird units and Galactic frame.
    """
    sc = SkyCoord(Galactic, c1, c2, unit=(unit1, unit2),
                  representation=representation)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2),
                               (attr1, attr2))

    sc = SkyCoord(1000*c1*u.Unit(unit1/1000), cls2(c2, unit=unit2),
                  Galactic,
                  unit=(unit1, unit2, unit3), representation=representation)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2),
                               (attr1, attr2))

    kwargs = {attr1: c1, attr2: c2}
    sc = SkyCoord(Galactic, unit=(unit1, unit2),
                  representation=representation, **kwargs)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2),
                               (attr1, attr2))


@pytest.mark.parametrize(units_attr_args,
                         (x for x in units_attr_sets if x[0] != 'unitspherical'))
def test_galactic_three_components(repr_name, unit1, unit2, unit3, cls2, attr1, attr2, attr3,
                                   representation, c1, c2, c3):
    """
    Tests positional inputs using components (COMP1, COMP2, COMP3)
    and various representations.  Use weird units and Galactic frame.
    """
    sc = Galactic(1000*c1*u.Unit(unit1/1000), cls2(c2, unit=unit2),
                  1000*c3*u.Unit(unit3/1000), representation=representation)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2, c3*unit3),
                               (attr1, attr2, attr3))

    kwargs = {attr3: c3*unit3}
    sc = Galactic(c1*unit1, c2*unit2,
                  representation=representation, **kwargs)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2, c3*unit3),
                               (attr1, attr2, attr3))

    kwargs = {attr1: c1*unit1, attr2: c2*unit2, attr3: c3*unit3}
    sc = Galactic(representation=representation, **kwargs)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2, c3*unit3),
                               (attr1, attr2, attr3))


@pytest.mark.parametrize(units_attr_args,
                         (x for x in units_attr_sets
                          if x[0] in ('spherical', 'unitspherical')))
def test_galactic_spherical_two_components(repr_name, unit1, unit2, unit3, cls2,
                                           attr1, attr2, attr3, representation, c1, c2, c3):
    """
    Tests positional inputs using components (COMP1, COMP2) for spherical
    representations.  Use weird units and Galactic frame.
    """

    sc = Galactic(1000*c1*u.Unit(unit1/1000), cls2(c2, unit=unit2), representation=representation)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2), (attr1, attr2))

    sc = Galactic(c1*unit1, c2*unit2, representation=representation)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2), (attr1, attr2))

    kwargs = {attr1: c1*unit1, attr2: c2*unit2}
    sc = Galactic(representation=representation, **kwargs)
    assert_quantities_allclose(sc, (c1*unit1, c2*unit2), (attr1, attr2))


@pytest.mark.parametrize('repr_name,unit1,unit2,unit3,cls2,attr1,attr2,attr3',
                         (x for x in base_unit_attr_sets if x[0] != 'unitspherical'))
def test_skycoord_coordinate_input(repr_name, unit1, unit2, unit3, cls2, attr1, attr2, attr3):
    c1, c2, c3 = 1, 2, 3
    sc = SkyCoord([(c1, c2, c3)], unit=(unit1, unit2, unit3), representation=repr_name,
                  frame='galactic')
    assert_quantities_allclose(sc, ([c1]*unit1, [c2]*unit2, [c3]*unit3), (attr1, attr2, attr3))

    c1, c2, c3 = 1*unit1, 2*unit2, 3*unit3
    sc = SkyCoord([(c1, c2, c3)], representation=repr_name, frame='galactic')
    assert_quantities_allclose(sc, ([1]*unit1, [2]*unit2, [3]*unit3), (attr1, attr2, attr3))


def test_skycoord_string_coordinate_input():
    sc = SkyCoord('01 02 03 +02 03 04', unit='deg', representation='unitspherical')
    assert_quantities_allclose(sc, (Angle('01:02:03', unit='deg'),
                                    Angle('02:03:04', unit='deg')),
                               ('ra', 'dec'))
    sc = SkyCoord(['01 02 03 +02 03 04'], unit='deg', representation='unitspherical')
    assert_quantities_allclose(sc, (Angle(['01:02:03'], unit='deg'),
                                    Angle(['02:03:04'], unit='deg')),
                               ('ra', 'dec'))


def test_units():
    sc = SkyCoord(1, 2, 3, unit='m', representation='cartesian')  # All get meters
    assert sc.x.unit is u.m
    assert sc.y.unit is u.m
    assert sc.z.unit is u.m

    sc = SkyCoord(1, 2*u.km, 3, unit='m', representation='cartesian')  # All get u.m
    assert sc.x.unit is u.m
    assert sc.y.unit is u.m
    assert sc.z.unit is u.m

    sc = SkyCoord(1, 2, 3, unit=u.m, representation='cartesian')  # All get u.m
    assert sc.x.unit is u.m
    assert sc.y.unit is u.m
    assert sc.z.unit is u.m

    sc = SkyCoord(1, 2, 3, unit='m, km, pc', representation='cartesian')
    assert_quantities_allclose(sc, (1*u.m, 2*u.km, 3*u.pc), ('x', 'y', 'z'))

    with pytest.raises(u.UnitsError) as err:
        SkyCoord(1, 2, 3, unit=(u.m, u.m), representation='cartesian')
    assert 'should have matching physical types' in str(err)

    SkyCoord(1, 2, 3, unit=(u.m, u.km, u.pc), representation='cartesian')
    assert_quantities_allclose(sc, (1*u.m, 2*u.km, 3*u.pc), ('x', 'y', 'z'))


@pytest.mark.xfail
def test_units_known_fail():
    # should fail but doesn't => corner case oddity
    with pytest.raises(u.UnitsError):
        SkyCoord(1, 2, 3, unit=u.deg, representation='spherical')


def test_nodata_failure():
    with pytest.raises(ValueError):
        SkyCoord()


@pytest.mark.parametrize('mode, origin', [('wcs', 0),
                                          ('all', 0),
                                          ('all', 1)])
def test_wcs_methods(mode, origin):
    from ...wcs import WCS
    from ...utils.data import get_pkg_data_contents
    from ...wcs.utils import pixel_to_skycoord

    header = get_pkg_data_contents('../../wcs/tests/maps/1904-66_TAN.hdr', encoding='binary')
    wcs = WCS(header)

    ref = SkyCoord(0.1 * u.deg, -89. * u.deg, frame='icrs')

    xp, yp = ref.to_pixel(wcs, mode=mode, origin=origin)

    # WCS is in FK5 so we need to transform back to ICRS
    new = pixel_to_skycoord(xp, yp, wcs, mode=mode, origin=origin).transform_to('icrs')

    npt.assert_allclose(new.ra.degree, ref.ra.degree)
    npt.assert_allclose(new.dec.degree, ref.dec.degree)

    #also try to round-trip with `from_pixel`
    scnew = SkyCoord.from_pixel(xp, yp, wcs, mode=mode, origin=origin).transform_to('icrs')
    npt.assert_allclose(scnew.ra.degree, ref.ra.degree)
    npt.assert_allclose(scnew.dec.degree, ref.dec.degree)

    #Also make sure the right type comes out
    class SkyCoord2(SkyCoord):
        pass
    scnew2 = SkyCoord2.from_pixel(xp, yp, wcs, mode=mode, origin=origin)
    assert scnew.__class__ is SkyCoord
    assert scnew2.__class__ is SkyCoord2


def test_frame_attr_transform_inherit():
    """
    Test that frame attributes get inherited as expected during transform.
    Driven by #3106.
    """
    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK5)
    c2 = c.transform_to(FK4)
    assert c2.equinox.value == 'B1950.000'
    assert c2.obstime.value == 'B1950.000'

    c2 = c.transform_to(FK4(equinox='J1975', obstime='J1980'))
    assert c2.equinox.value == 'J1975.000'
    assert c2.obstime.value == 'J1980.000'

    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK4)
    c2 = c.transform_to(FK5)
    assert c2.equinox.value == 'J2000.000'
    assert c2.obstime is None

    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK4, obstime='J1980')
    c2 = c.transform_to(FK5)
    assert c2.equinox.value == 'J2000.000'
    assert c2.obstime.value == 'J1980.000'

    c = SkyCoord(1 * u.deg, 2 * u.deg, frame=FK4, equinox='J1975', obstime='J1980')
    c2 = c.transform_to(FK5)
    assert c2.equinox.value == 'J1975.000'
    assert c2.obstime.value == 'J1980.000'

    c2 = c.transform_to(FK5(equinox='J1990'))
    assert c2.equinox.value == 'J1990.000'
    assert c2.obstime.value == 'J1980.000'


def test_deepcopy():
    c1 = SkyCoord(1 * u.deg, 2 * u.deg)
    c2 = copy.copy(c1)
    c3 = copy.deepcopy(c1)

    c4 = SkyCoord([1, 2] * u.m, [2, 3] *u.m, [3, 4] * u.m, representation='cartesian', frame='fk5',
                  obstime='J1999.9', equinox='J1988.8')
    c5 = copy.deepcopy(c4)
    assert np.all(c5.x == c4.x)  # and y and z
    assert c5.frame.name == c4.frame.name
    assert c5.obstime == c4.obstime
    assert c5.equinox == c4.equinox
    assert c5.representation == c4.representation


def test_immutable():
    c1 = SkyCoord(1 * u.deg, 2 * u.deg)
    with pytest.raises(AttributeError):
        c1.ra = 3.0

    c1.foo = 42
    assert c1.foo == 42


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_search_around():
    """
    Test the search_around_* methods

    Here we don't actually test the values are right, just that the methods of
    SkyCoord work.  The accuracy tests are in ``test_matching.py``
    """
    from ...utils import NumpyRNGContext

    with NumpyRNGContext(987654321):
        sc1 = SkyCoord(np.random.rand(20) * 360.*u.degree,
                      (np.random.rand(20) * 180. - 90.)*u.degree)
        sc2 = SkyCoord(np.random.rand(100) * 360. * u.degree,
                      (np.random.rand(100) * 180. - 90.)*u.degree)

        sc1ds = SkyCoord(ra=sc1.ra, dec=sc1.dec, distance=np.random.rand(20)*u.kpc)
        sc2ds = SkyCoord(ra=sc2.ra, dec=sc2.dec, distance=np.random.rand(100)*u.kpc)

    idx1_sky, idx2_sky, d2d_sky, d3d_sky = sc1.search_around_sky(sc2, 10*u.deg)
    idx1_3d, idx2_3d, d2d_3d, d3d_3d = sc1ds.search_around_3d(sc2ds, 250*u.pc)


def test_init_with_frame_instance_keyword():

    # Frame instance
    c1 = SkyCoord(3 * u.deg, 4 * u.deg,
                  frame=FK5(equinox='J2010'))
    assert c1.equinox == Time('J2010')

    # Frame instance with data (data gets ignored)
    c2 = SkyCoord(3 * u.deg, 4 * u.deg,
                 frame=FK5(1. * u.deg, 2 * u.deg,
                 equinox='J2010'))
    assert c2.equinox == Time('J2010')
    assert allclose(c2.ra.degree, 3)
    assert allclose(c2.dec.degree, 4)

    # SkyCoord instance
    c3 = SkyCoord(3 * u.deg, 4 * u.deg, frame=c1)
    assert c3.equinox == Time('J2010')

    # Check duplicate arguments
    with pytest.raises(ValueError) as exc:
        c = SkyCoord(3 * u.deg, 4 * u.deg, frame=FK5(equinox='J2010'), equinox='J2001')
    assert exc.value.args[0] == ("cannot specify frame attribute "
                                 "'equinox' directly in SkyCoord "
                                 "since a frame instance was passed in")


def test_init_with_frame_instance_positional():

    # Frame instance
    with pytest.raises(ValueError) as exc:
        c1 = SkyCoord(3 * u.deg, 4 * u.deg, FK5(equinox='J2010'))
    assert exc.value.args[0] == ("FK5 instance cannot be passed as a "
                                 "positional argument for the frame, "
                                 "pass it using the frame= keyword "
                                 "instead.")

    # Positional frame instance with data raises exception
    with pytest.raises(ValueError) as exc:
        SkyCoord(3 * u.deg, 4 * u.deg, FK5(1. * u.deg, 2 * u.deg, equinox='J2010'))
    assert exc.value.args[0] == ("FK5 instance cannot be passed as a "
                                 "positional argument for the frame, "
                                 "pass it using the frame= keyword "
                                 "instead.")

    # Positional SkyCoord instance (for frame) raises exception
    with pytest.raises(ValueError) as exc:
        SkyCoord(3 * u.deg, 4 * u.deg, SkyCoord(1. * u.deg, 2 * u.deg, equinox='J2010'))
    assert exc.value.args[0] == ("SkyCoord instance cannot be passed as a "
                                 "positional argument for the frame, "
                                 "pass it using the frame= keyword "
                                 "instead.")


def test_guess_from_table():
    from ...table import Table, Column
    from ...utils import NumpyRNGContext

    tab = Table()
    with NumpyRNGContext(987654321):
        tab.add_column(Column(data=np.random.rand(1000),unit='deg',name='RA[J2000]'))
        tab.add_column(Column(data=np.random.rand(1000),unit='deg',name='DEC[J2000]'))

    sc = SkyCoord.guess_from_table(tab)
    npt.assert_array_equal(sc.ra.deg, tab['RA[J2000]'])
    npt.assert_array_equal(sc.dec.deg, tab['DEC[J2000]'])

    # try without units in the table
    tab['RA[J2000]'].unit = None
    tab['DEC[J2000]'].unit = None
    # should fail if not given explicitly
    with pytest.raises(u.UnitsError):
        sc2 = SkyCoord.guess_from_table(tab)

    # but should work if provided
    sc2 = SkyCoord.guess_from_table(tab, unit=u.deg)
    npt.assert_array_equal(sc.ra.deg, tab['RA[J2000]'])
    npt.assert_array_equal(sc.dec.deg, tab['DEC[J2000]'])

    # should fail if two options are available - ambiguity bad!
    tab.add_column(Column(data=np.random.rand(1000), name='RA_J1900'))
    with pytest.raises(ValueError) as excinfo:
        sc3 = SkyCoord.guess_from_table(tab, unit=u.deg)
    assert 'J1900' in excinfo.value.args[0] and 'J2000' in excinfo.value.args[0]

    # should also fail if user specifies something already in the table, but
    # should succeed even if the user has to give one of the components
    tab.remove_column('RA_J1900')
    with pytest.raises(ValueError):
        sc3 = SkyCoord.guess_from_table(tab, ra=tab['RA[J2000]'], unit=u.deg)

    oldra = tab['RA[J2000]']
    tab.remove_column('RA[J2000]')
    sc3 = SkyCoord.guess_from_table(tab, ra=oldra, unit=u.deg)
    npt.assert_array_equal(sc3.ra.deg, oldra)
    npt.assert_array_equal(sc3.dec.deg, tab['DEC[J2000]'])

    # check a few non-ICRS/spherical systems
    x, y, z = np.arange(3).reshape(3, 1) * u.pc
    l, b = np.arange(2).reshape(2, 1) * u.deg

    tabcart = Table([x, y, z], names=('x', 'y', 'z'))
    tabgal = Table([b, l], names=('b', 'l'))

    sc_cart = SkyCoord.guess_from_table(tabcart, representation='cartesian')
    npt.assert_array_equal(sc_cart.x, x)
    npt.assert_array_equal(sc_cart.y, y)
    npt.assert_array_equal(sc_cart.z, z)

    sc_gal = SkyCoord.guess_from_table(tabgal, frame='galactic')
    npt.assert_array_equal(sc_gal.l, l)
    npt.assert_array_equal(sc_gal.b, b)

    #also try some column names that *end* with the attribute name
    tabgal['b'].name = 'gal_b'
    tabgal['l'].name = 'gal_l'
    SkyCoord.guess_from_table(tabgal, frame='galactic')

    tabgal['gal_b'].name = 'blob'
    tabgal['gal_l'].name = 'central'
    with pytest.raises(ValueError):
        SkyCoord.guess_from_table(tabgal, frame='galactic')
