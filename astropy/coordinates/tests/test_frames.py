# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_allclose

from ... import units as u
from ...tests.helper import pytest
from .. import representation

def test_create_data_frames():
    from ..builtin_frames import ICRS

    #from repr
    i1 = ICRS(representation.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc))
    i2 = ICRS(representation.UnitSphericalRepresentation(lon=1*u.deg, lat=2*u.deg))

    #from preferred name
    i3 = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.kpc)
    i4 = ICRS(ra=1*u.deg, dec=2*u.deg)

    assert i1.data.lat == i3.data.lat
    assert i1.data.lon == i3.data.lon
    assert i1.data.distance == i3.data.distance

    assert i2.data.lat == i4.data.lat
    assert i2.data.lon == i4.data.lon

    #now make sure the preferred names work as properties
    assert_allclose(i1.ra, i3.ra)
    assert_allclose(i2.ra, i4.ra)
    assert_allclose(i1.distance, i3.distance)

    with pytest.raises(AttributeError):
        i1.ra = [11.]*u.deg


def test_create_orderered_data():
    from ..builtin_frames import ICRS, Galactic, AltAz

    TOL = 1e-10*u.deg

    i = ICRS(1*u.deg, 2*u.deg)
    assert (i.ra - 1*u.deg) < TOL
    assert (i.dec - 2*u.deg) < TOL

    g = Galactic(1*u.deg, 2*u.deg)
    assert (g.l - 1*u.deg) < TOL
    assert (g.b - 2*u.deg) < TOL

    a = AltAz(1*u.deg, 2*u.deg)
    assert (a.az - 1*u.deg) < TOL
    assert (a.alt - 2*u.deg) < TOL

    with pytest.raises(TypeError):
        ICRS(1*u.deg, 2*u.deg, 1*u.deg, 2*u.deg)

    with pytest.raises(TypeError):
        sph = representation.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc)
        ICRS(sph, 1*u.deg, 2*u.deg)


def test_create_nodata_frames():
    from ..builtin_frames import ICRS, FK4, FK5

    i = ICRS()
    assert len(i.frame_attr_names) == 0

    f5 = FK5()
    assert f5.equinox == FK5.frame_attr_names['equinox']

    f4 = FK4()
    assert f4.equinox == FK4.frame_attr_names['equinox']

    #obstime is special because it's a property that uses equinox if obstime is not set
    assert f4.obstime in (FK4.frame_attr_names['obstime'], FK4.frame_attr_names['equinox'])


def test_frame_repr():
    from ..builtin_frames import ICRS, FK5

    i = ICRS()
    assert repr(i) == '<ICRS Frame>'

    f5 = FK5()
    assert repr(f5).startswith('<FK5 Frame: equinox=')

    i2 = ICRS(ra=1*u.deg, dec=2*u.deg)
    i3 = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.kpc)

    assert repr(i2) == '<ICRS Coordinate: ra=1.0 deg, dec=2.0 deg>'
    assert repr(i3) == ('<ICRS Coordinate: ra=1.0 deg, dec=2.0 deg, '
                        'distance=3.0 kpc>')


def test_converting_units():
    import re
    from ..builtin_frames import ICRS, FK5
    from ..representation import SphericalRepresentation

    # this is a regular expression that with split (see below) removes what's
    # the decimal point  to fix rounding problems
    rexrepr = re.compile(r'(.*?=\d\.).*?( .*?=\d\.).*?( .*)')

    # Use values that aren't subject to rounding down to X.9999...
    i2 = ICRS(ra=1.1*u.deg, dec=2.1*u.deg)

    #converting from FK5 to ICRS and back changes the *internal* representation,
    # but it should still come out in the preferred form

    i4 = i2.transform_to(FK5).transform_to(ICRS)
    ri2 = ''.join(rexrepr.split(repr(i2)))
    ri4 = ''.join(rexrepr.split(repr(i4)))
    assert ri2 == ri4
    assert i2.data.lon.unit != i4.data.lon.unit  # Internal repr changed

    #but that *shouldn't* hold if we turn off units for the representation
    class FakeICRS(ICRS):
        _frame_specific_representation_info = {
            'spherical': {'names': ('ra', 'dec', 'distance'),
                          'units': (None, None, None)},
            'unitspherical': {'names': ('ra', 'dec'),
                              'units': (None, None)}}


    fi = FakeICRS(i4.data)
    ri2 = ''.join(rexrepr.split(repr(i2)))
    rfi = ''.join(rexrepr.split(repr(fi)))
    rfi = re.sub('FakeICRS', 'ICRS', rfi)  # Force frame name to match
    assert ri2 != rfi

    # the attributes should also get the right units
    assert i2.ra.unit == i4.ra.unit
    # unless no units
    assert i2.ra.unit != fi.ra.unit


def test_realizing():
    from ..builtin_frames import ICRS, FK5
    from ...time import Time

    rep = representation.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc)

    i = ICRS()
    i2 = i.realize_frame(rep)

    assert not i.has_data
    assert i2.has_data

    f = FK5(equinox=Time('J2001', scale='utc'))
    f2 = f.realize_frame(rep)

    assert not f.has_data
    assert f2.has_data

    assert f2.equinox == f.equinox
    assert f2.equinox != FK5.frame_attr_names['equinox']


def test_getitem():
    from ..builtin_frames import ICRS

    rep = representation.SphericalRepresentation(
        [1, 2, 3]*u.deg, [4, 5, 6]*u.deg, [7, 8, 9]*u.kpc)

    i = ICRS(rep)
    assert len(i.ra) == 3

    iidx = i[1:]
    assert len(iidx.ra) == 2

    iidx2 = i[0]
    assert iidx2.ra.isscalar

def test_transform():
    """
    This test just makes sure the transform architecture works, but does *not*
    actually test all the builtin transforms themselves are accurate
    """
    from ..builtin_frames import ICRS, FK4, FK5, Galactic
    from ...time import Time

    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg)
    f = i.transform_to(FK5)
    i2 = f.transform_to(ICRS)

    assert i2.data.__class__ == representation.UnitSphericalRepresentation

    assert_allclose(i.ra, i2.ra)
    assert_allclose(i.dec, i2.dec)


    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[5, 6]*u.kpc)
    f = i.transform_to(FK5)
    i2 = f.transform_to(ICRS)

    assert i2.data.__class__ != representation.UnitSphericalRepresentation


    f = FK5(ra=1*u.deg, dec=2*u.deg, equinox=Time('J2001', scale='utc'))
    f4 = f.transform_to(FK4)
    f4_2 = f.transform_to(FK4(equinox=f.equinox))

    #make sure attributes are copied over correctly
    assert f4.equinox == FK4.frame_attr_names['equinox']
    assert f4_2.equinox == f.equinox


    #make sure self-transforms also work
    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg)
    i2 = i.transform_to(ICRS)

    assert_allclose(i.ra, i2.ra)
    assert_allclose(i.dec, i2.dec)

    f = FK5(ra=1*u.deg, dec=2*u.deg, equinox=Time('J2001', scale='utc'))
    f2 = f.transform_to(FK5)  # default equinox, so should be *different*
    assert f2.equinox == FK5().equinox
    with pytest.raises(AssertionError):
        assert_allclose(f.ra, f2.ra)
    with pytest.raises(AssertionError):
        assert_allclose(f.dec, f2.dec)


    #finally, check Galactic round-tripping
    i1 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg)
    i2 = i1.transform_to(Galactic).transform_to(ICRS)

    assert_allclose(i1.ra, i2.ra)
    assert_allclose(i1.dec, i2.dec)

def test_sep():
    from ..builtin_frames import ICRS

    i1 = ICRS(ra=0*u.deg, dec=1*u.deg)
    i2 = ICRS(ra=0*u.deg, dec=2*u.deg)

    sep = i1.separation(i2)
    assert sep.deg == 1

    i3 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[5, 6]*u.kpc)
    i4 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[4, 5]*u.kpc)

    sep3d = i3.separation_3d(i4)
    assert_allclose(sep3d.to(u.kpc).value, np.array([1, 1]))


def test_time_inputs():
    """
    Test validation and conversion of inputs for equinox and obstime attributes.
    """
    from ...time import Time
    from ..builtin_frames import FK4

    c = FK4(1 * u.deg, 2 * u.deg, equinox='J2001.5', obstime='2000-01-01 12:00:00')
    assert c.equinox == Time('J2001.5')
    assert c.obstime == Time('2000-01-01 12:00:00')

    with pytest.raises(ValueError) as err:
        c = FK4(1 * u.deg, 2 * u.deg, equinox=1.5)
    assert 'Invalid time input' in str(err)

    with pytest.raises(ValueError) as err:
        c = FK4(1 * u.deg, 2 * u.deg, obstime='hello')
    assert 'Invalid time input' in str(err)

    with pytest.raises(ValueError) as err:
        c = FK4(1 * u.deg, 2 * u.deg, obstime=['J2000', 'J2001'])
    assert "must be a single (scalar) value" in str(err)


def test_is_frame_attr_default():
    """
    Check that the `is_frame_attr_default` machinery works as expected
    """
    from ...time import Time
    from ..builtin_frames import FK5

    c1 = FK5(ra=1*u.deg, dec=1*u.deg)
    c2 = FK5(ra=1*u.deg, dec=1*u.deg, equinox=FK5.frame_attr_names['equinox'])
    c3 = FK5(ra=1*u.deg, dec=1*u.deg, equinox=Time('J2001.5'))

    assert c1.equinox == c2.equinox
    assert c1.equinox != c3.equinox

    assert c1.is_frame_attr_default('equinox')
    assert not c2.is_frame_attr_default('equinox')
    assert not c3.is_frame_attr_default('equinox')

    c4 = c1.realize_frame(representation.UnitSphericalRepresentation(3*u.deg, 4*u.deg))
    c5 = c2.realize_frame(representation.UnitSphericalRepresentation(3*u.deg, 4*u.deg))

    assert c4.is_frame_attr_default('equinox')
    assert not c5.is_frame_attr_default('equinox')


def test_altaz_attributes():
    from ...time import Time
    from .. import EarthLocation, AltAz

    aa = AltAz(1*u.deg, 2*u.deg)
    assert aa.obstime is None
    assert aa.location is None

    aa2 = AltAz(1*u.deg, 2*u.deg, obstime='J2000')
    assert aa2.obstime == Time('J2000')

    aa3 = AltAz(1*u.deg, 2*u.deg, location=EarthLocation(0*u.deg, 0*u.deg, 0*u.m))
    assert isinstance(aa3.location, EarthLocation)


def test_representation():
    """
    Test the getter and setter properties for `representation`
    """
    from ..builtin_frames import ICRS

    # Create the frame object.
    icrs = ICRS(ra=1*u.deg, dec=1*u.deg)
    data = icrs.data

    # Create some representation objects.
    icrs_cart = icrs.cartesian
    icrs_spher = icrs.spherical

    # Testing when `_representation` set to `CartesianRepresentation`.
    icrs.representation = representation.CartesianRepresentation
    
    assert icrs.representation == representation.CartesianRepresentation
    assert icrs_cart.x == icrs.x
    assert icrs_cart.y == icrs.y
    assert icrs_cart.z == icrs.z
    assert icrs.data == data

    # Testing that an ICRS object in CartesianRepresentation must not have spherical attributes.
    for attr in ('ra', 'dec', 'distance'):
        with pytest.raises(AttributeError) as err:
            getattr(icrs, attr)
        assert 'object has no attribute' in str(err)

    # Testing when `_representation` set to `CylindricalRepresentation`.
    icrs.representation = representation.CylindricalRepresentation

    assert icrs.representation == representation.CylindricalRepresentation
    assert icrs.data == data

    # Testing setter input using text argument for spherical.
    icrs.representation = 'spherical'

    assert icrs.representation is representation.SphericalRepresentation
    assert icrs_spher.lat == icrs.dec
    assert icrs_spher.lon == icrs.ra
    assert icrs_spher.distance == icrs.distance
    assert icrs.data == data

    # Testing that an ICRS object in SphericalRepresentation must not have cartesian attributes.
    for attr in ('x', 'y', 'z'):
        with pytest.raises(AttributeError) as err:
            getattr(icrs, attr)
        assert 'object has no attribute' in str(err)

    # Testing setter input using text argument for cylindrical.
    icrs.representation = 'cylindrical'

    assert icrs.representation is representation.CylindricalRepresentation
    assert icrs.data == data

    with pytest.raises(ValueError) as err:
        icrs.representation = 'WRONG'
    assert 'but must be a BaseRepresentation class' in str(err)

    with pytest.raises(ValueError) as err:
        icrs.representation = ICRS
    assert 'but must be a BaseRepresentation class' in str(err)


def test_dynamic_attrs():
    from ..builtin_frames import ICRS
    c = ICRS(1*u.deg, 2*u.deg)
    assert 'ra' in dir(c)
    assert 'dec' in dir(c)

    with pytest.raises(AttributeError) as err:
        c.blahblah
    assert "object has no attribute 'blahblah'" in str(err)

    with pytest.raises(AttributeError) as err:
        c.ra = 1
    assert "Cannot set any frame attribute" in str(err)

    c.blahblah = 1
    assert c.blahblah == 1
