# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


from copy import deepcopy
import numpy as np

from astropy import units as u
from astropy.tests.helper import (catch_warnings, pytest,
                             assert_quantity_allclose as assert_allclose)
from astropy.utils import OrderedDescriptorContainer
from astropy.utils.compat import NUMPY_LT_1_14
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import representation as r
from astropy.coordinates.representation import REPRESENTATION_CLASSES
from astropy.units import allclose


from .test_representation import unitphysics  # this fixture is used below


def setup_function(func):
    func.REPRESENTATION_CLASSES_ORIG = deepcopy(REPRESENTATION_CLASSES)


def teardown_function(func):
    REPRESENTATION_CLASSES.clear()
    REPRESENTATION_CLASSES.update(func.REPRESENTATION_CLASSES_ORIG)


def test_frame_attribute_descriptor():
    """ Unit tests of the Attribute descriptor """
    from astropy.coordinates.attributes import Attribute

    class TestAttributes(metaclass=OrderedDescriptorContainer):
        attr_none = Attribute()
        attr_2 = Attribute(default=2)
        attr_3_attr2 = Attribute(default=3, secondary_attribute='attr_2')
        attr_none_attr2 = Attribute(default=None, secondary_attribute='attr_2')
        attr_none_nonexist = Attribute(default=None, secondary_attribute='nonexist')

    t = TestAttributes()

    # Defaults
    assert t.attr_none is None
    assert t.attr_2 == 2
    assert t.attr_3_attr2 == 3
    assert t.attr_none_attr2 == t.attr_2
    assert t.attr_none_nonexist is None  # No default and non-existent secondary attr

    # Setting values via '_'-prefixed internal vars (as would normally done in __init__)
    t._attr_none = 10
    assert t.attr_none == 10

    t._attr_2 = 20
    assert t.attr_2 == 20
    assert t.attr_3_attr2 == 3
    assert t.attr_none_attr2 == t.attr_2

    t._attr_none_attr2 = 40
    assert t.attr_none_attr2 == 40

    # Make sure setting values via public attribute fails
    with pytest.raises(AttributeError) as err:
        t.attr_none = 5
    assert 'Cannot set frame attribute' in str(err)


def test_frame_subclass_attribute_descriptor():
    from astropy.coordinates.builtin_frames import FK4
    from astropy.coordinates.attributes import Attribute, TimeAttribute
    from astropy.time import Time

    _EQUINOX_B1980 = Time('B1980', scale='tai')

    class MyFK4(FK4):
        # equinox inherited from FK4, obstime overridden, and newattr is new
        obstime = TimeAttribute(default=_EQUINOX_B1980)
        newattr = Attribute(default='newattr')

    mfk4 = MyFK4()
    assert mfk4.equinox.value == 'B1950.000'
    assert mfk4.obstime.value == 'B1980.000'
    assert mfk4.newattr == 'newattr'

    assert set(mfk4.get_frame_attr_names()) == set(['equinox', 'obstime', 'newattr'])

    mfk4 = MyFK4(equinox='J1980.0', obstime='J1990.0', newattr='world')
    assert mfk4.equinox.value == 'J1980.000'
    assert mfk4.obstime.value == 'J1990.000'
    assert mfk4.newattr == 'world'


def test_create_data_frames():
    from astropy.coordinates.builtin_frames import ICRS

    # from repr
    i1 = ICRS(r.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc))
    i2 = ICRS(r.UnitSphericalRepresentation(lon=1*u.deg, lat=2*u.deg))

    # from preferred name
    i3 = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.kpc)
    i4 = ICRS(ra=1*u.deg, dec=2*u.deg)

    assert i1.data.lat == i3.data.lat
    assert i1.data.lon == i3.data.lon
    assert i1.data.distance == i3.data.distance

    assert i2.data.lat == i4.data.lat
    assert i2.data.lon == i4.data.lon

    # now make sure the preferred names work as properties
    assert_allclose(i1.ra, i3.ra)
    assert_allclose(i2.ra, i4.ra)
    assert_allclose(i1.distance, i3.distance)

    with pytest.raises(AttributeError):
        i1.ra = [11.]*u.deg


def test_create_orderered_data():
    from astropy.coordinates.builtin_frames import ICRS, Galactic, AltAz

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
        sph = r.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc)
        ICRS(sph, 1*u.deg, 2*u.deg)


def test_create_nodata_frames():
    from astropy.coordinates.builtin_frames import ICRS, FK4, FK5

    i = ICRS()
    assert len(i.get_frame_attr_names()) == 0

    f5 = FK5()
    assert f5.equinox == FK5.get_frame_attr_names()['equinox']

    f4 = FK4()
    assert f4.equinox == FK4.get_frame_attr_names()['equinox']

    # obstime is special because it's a property that uses equinox if obstime is not set
    assert f4.obstime in (FK4.get_frame_attr_names()['obstime'],
                          FK4.get_frame_attr_names()['equinox'])


def test_no_data_nonscalar_frames():
    from astropy.coordinates.builtin_frames import AltAz
    from astropy.time import Time
    a1 = AltAz(obstime=Time('2012-01-01') + np.arange(10.) * u.day,
               temperature=np.ones((3, 1)) * u.deg_C)
    assert a1.obstime.shape == (3, 10)
    assert a1.temperature.shape == (3, 10)
    assert a1.shape == (3, 10)
    with pytest.raises(ValueError) as exc:
        AltAz(obstime=Time('2012-01-01') + np.arange(10.) * u.day,
              temperature=np.ones((3,)) * u.deg_C)
    assert 'inconsistent shapes' in str(exc)


def test_frame_repr():
    from astropy.coordinates.builtin_frames import ICRS, FK5

    i = ICRS()
    assert repr(i) == '<ICRS Frame>'

    f5 = FK5()
    assert repr(f5).startswith('<FK5 Frame (equinox=')

    i2 = ICRS(ra=1*u.deg, dec=2*u.deg)
    i3 = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.kpc)

    assert repr(i2) == ('<ICRS Coordinate: (ra, dec) in deg\n'
                        '    ({})>').format(' 1.,  2.' if NUMPY_LT_1_14
                                             else '1., 2.')
    assert repr(i3) == ('<ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)\n'
                        '    ({})>').format(' 1.,  2.,  3.' if NUMPY_LT_1_14
                                            else '1., 2., 3.')

    # try with arrays
    i2 = ICRS(ra=[1.1, 2.1]*u.deg, dec=[2.1, 3.1]*u.deg)
    i3 = ICRS(ra=[1.1, 2.1]*u.deg, dec=[-15.6, 17.1]*u.deg, distance=[11., 21.]*u.kpc)

    assert repr(i2) == ('<ICRS Coordinate: (ra, dec) in deg\n'
                        '    [{}]>').format('( 1.1,  2.1), ( 2.1,  3.1)'
                                            if NUMPY_LT_1_14 else
                                            '(1.1, 2.1), (2.1, 3.1)')

    if NUMPY_LT_1_14:
        assert repr(i3) == ('<ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)\n'
                            '    [( 1.1, -15.6,  11.), ( 2.1,  17.1,  21.)]>')
    else:
        assert repr(i3) == ('<ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)\n'
                            '    [(1.1, -15.6, 11.), (2.1,  17.1, 21.)]>')


def test_frame_repr_vels():
    from astropy.coordinates.builtin_frames import ICRS

    i = ICRS(ra=1*u.deg, dec=2*u.deg,
             pm_ra_cosdec=1*u.marcsec/u.yr, pm_dec=2*u.marcsec/u.yr)

    # unit comes out as mas/yr because of the preferred units defined in the
    # frame RepresentationMapping
    assert repr(i) == ('<ICRS Coordinate: (ra, dec) in deg\n'
                       '    ({0})\n'
                       ' (pm_ra_cosdec, pm_dec) in mas / yr\n'
                       '    ({0})>').format(' 1.,  2.' if NUMPY_LT_1_14 else
                                            '1., 2.')


def test_converting_units():
    import re
    from astropy.coordinates.baseframe import RepresentationMapping
    from astropy.coordinates.builtin_frames import ICRS, FK5

    # this is a regular expression that with split (see below) removes what's
    # the decimal point  to fix rounding problems
    rexrepr = re.compile(r'(.*?=\d\.).*?( .*?=\d\.).*?( .*)')

    # Use values that aren't subject to rounding down to X.9999...
    i2 = ICRS(ra=2.*u.deg, dec=2.*u.deg)
    i2_many = ICRS(ra=[2., 4.]*u.deg, dec=[2., -8.1]*u.deg)

    # converting from FK5 to ICRS and back changes the *internal* representation,
    # but it should still come out in the preferred form

    i4 = i2.transform_to(FK5).transform_to(ICRS)
    i4_many = i2_many.transform_to(FK5).transform_to(ICRS)

    ri2 = ''.join(rexrepr.split(repr(i2)))
    ri4 = ''.join(rexrepr.split(repr(i4)))
    assert ri2 == ri4
    assert i2.data.lon.unit != i4.data.lon.unit  # Internal repr changed

    ri2_many = ''.join(rexrepr.split(repr(i2_many)))
    ri4_many = ''.join(rexrepr.split(repr(i4_many)))

    assert ri2_many == ri4_many
    assert i2_many.data.lon.unit != i4_many.data.lon.unit  # Internal repr changed

    # but that *shouldn't* hold if we turn off units for the representation
    class FakeICRS(ICRS):
        frame_specific_representation_info = {
            'spherical': [RepresentationMapping('lon', 'ra', u.hourangle),
                          RepresentationMapping('lat', 'dec', None),
                          RepresentationMapping('distance', 'distance')]  # should fall back to default of None unit
        }

    fi = FakeICRS(i4.data)
    ri2 = ''.join(rexrepr.split(repr(i2)))
    rfi = ''.join(rexrepr.split(repr(fi)))
    rfi = re.sub('FakeICRS', 'ICRS', rfi)  # Force frame name to match
    assert ri2 != rfi

    # the attributes should also get the right units
    assert i2.dec.unit == i4.dec.unit
    # unless no/explicitly given units
    assert i2.dec.unit != fi.dec.unit
    assert i2.ra.unit != fi.ra.unit
    assert fi.ra.unit == u.hourangle


def test_representation_info():
    from astropy.coordinates.baseframe import RepresentationMapping
    from astropy.coordinates.builtin_frames import ICRS

    class NewICRS1(ICRS):
        frame_specific_representation_info = {
            r.SphericalRepresentation: [
                RepresentationMapping('lon', 'rara', u.hourangle),
                RepresentationMapping('lat', 'decdec', u.degree),
                RepresentationMapping('distance', 'distance', u.kpc)]
        }

    i1 = NewICRS1(rara=10*u.degree, decdec=-12*u.deg, distance=1000*u.pc,
                  pm_rara_cosdecdec=100*u.mas/u.yr,
                  pm_decdec=17*u.mas/u.yr,
                  radial_velocity=10*u.km/u.s)
    assert allclose(i1.rara, 10*u.deg)
    assert i1.rara.unit == u.hourangle
    assert allclose(i1.decdec, -12*u.deg)
    assert allclose(i1.distance, 1000*u.pc)
    assert i1.distance.unit == u.kpc
    assert allclose(i1.pm_rara_cosdecdec, 100*u.mas/u.yr)
    assert allclose(i1.pm_decdec, 17*u.mas/u.yr)

    # this should auto-set the names of UnitSpherical:
    i1.set_representation_cls(r.UnitSphericalRepresentation,
                              s=r.UnitSphericalCosLatDifferential)
    assert allclose(i1.rara, 10*u.deg)
    assert allclose(i1.decdec, -12*u.deg)
    assert allclose(i1.pm_rara_cosdecdec, 100*u.mas/u.yr)
    assert allclose(i1.pm_decdec, 17*u.mas/u.yr)

    # For backwards compatibility, we also support the string name in the
    # representation info dictionary:
    class NewICRS2(ICRS):
        frame_specific_representation_info = {
            'spherical': [
                RepresentationMapping('lon', 'ang1', u.hourangle),
                RepresentationMapping('lat', 'ang2', u.degree),
                RepresentationMapping('distance', 'howfar', u.kpc)]
        }

    i2 = NewICRS2(ang1=10*u.degree, ang2=-12*u.deg, howfar=1000*u.pc)
    assert allclose(i2.ang1, 10*u.deg)
    assert i2.ang1.unit == u.hourangle
    assert allclose(i2.ang2, -12*u.deg)
    assert allclose(i2.howfar, 1000*u.pc)
    assert i2.howfar.unit == u.kpc

    # Test that the differential kwargs get overridden
    class NewICRS3(ICRS):
        frame_specific_representation_info = {
            r.SphericalCosLatDifferential: [
                RepresentationMapping('d_lon_coslat', 'pm_ang1', u.hourangle/u.year),
                RepresentationMapping('d_lat', 'pm_ang2'),
                RepresentationMapping('d_distance', 'vlos', u.kpc/u.Myr)]
        }

    i3 = NewICRS3(lon=10*u.degree, lat=-12*u.deg, distance=1000*u.pc,
                  pm_ang1=1*u.mas/u.yr, pm_ang2=2*u.mas/u.yr,
                  vlos=100*u.km/u.s)
    assert allclose(i3.pm_ang1, 1*u.mas/u.yr)
    assert i3.pm_ang1.unit == u.hourangle/u.year
    assert allclose(i3.pm_ang2, 2*u.mas/u.yr)
    assert allclose(i3.vlos, 100*u.km/u.s)
    assert i3.vlos.unit == u.kpc/u.Myr


def test_realizing():
    from astropy.coordinates.builtin_frames import ICRS, FK5
    from astropy.time import Time

    rep = r.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc)

    i = ICRS()
    i2 = i.realize_frame(rep)

    assert not i.has_data
    assert i2.has_data

    f = FK5(equinox=Time('J2001', scale='utc'))
    f2 = f.realize_frame(rep)

    assert not f.has_data
    assert f2.has_data

    assert f2.equinox == f.equinox
    assert f2.equinox != FK5.get_frame_attr_names()['equinox']

    # Check that a nicer error message is returned:
    with pytest.raises(TypeError) as excinfo:
        f.realize_frame(f.representation_type)

    assert ('Class passed as data instead of a representation' in
            excinfo.value.args[0])

def test_replicating():
    from astropy.coordinates.builtin_frames import ICRS, AltAz
    from astropy.time import Time

    i = ICRS(ra=[1]*u.deg, dec=[2]*u.deg)

    icopy = i.replicate(copy=True)
    irepl = i.replicate(copy=False)
    i.data._lat[:] = 0*u.deg
    assert np.all(i.data.lat == irepl.data.lat)
    assert np.all(i.data.lat != icopy.data.lat)

    iclone = i.replicate_without_data()
    assert i.has_data
    assert not iclone.has_data

    aa = AltAz(alt=1*u.deg, az=2*u.deg, obstime=Time('J2000'))
    aaclone = aa.replicate_without_data(obstime=Time('J2001'))
    assert not aaclone.has_data
    assert aa.obstime != aaclone.obstime
    assert aa.pressure == aaclone.pressure
    assert aa.obswl == aaclone.obswl


def test_getitem():
    from astropy.coordinates.builtin_frames import ICRS

    rep = r.SphericalRepresentation(
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
    from astropy.coordinates.builtin_frames import ICRS, FK4, FK5, Galactic
    from astropy.time import Time

    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg)
    f = i.transform_to(FK5)
    i2 = f.transform_to(ICRS)

    assert i2.data.__class__ == r.UnitSphericalRepresentation

    assert_allclose(i.ra, i2.ra)
    assert_allclose(i.dec, i2.dec)

    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[5, 6]*u.kpc)
    f = i.transform_to(FK5)
    i2 = f.transform_to(ICRS)

    assert i2.data.__class__ != r.UnitSphericalRepresentation

    f = FK5(ra=1*u.deg, dec=2*u.deg, equinox=Time('J2001', scale='utc'))
    f4 = f.transform_to(FK4)
    f4_2 = f.transform_to(FK4(equinox=f.equinox))

    # make sure attributes are copied over correctly
    assert f4.equinox == FK4.get_frame_attr_names()['equinox']
    assert f4_2.equinox == f.equinox

    # make sure self-transforms also work
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

    # finally, check Galactic round-tripping
    i1 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg)
    i2 = i1.transform_to(Galactic).transform_to(ICRS)

    assert_allclose(i1.ra, i2.ra)
    assert_allclose(i1.dec, i2.dec)


def test_transform_to_nonscalar_nodata_frame():
    # https://github.com/astropy/astropy/pull/5254#issuecomment-241592353
    from astropy.coordinates.builtin_frames import ICRS, FK5
    from astropy.time import Time
    times = Time('2016-08-23') + np.linspace(0, 10, 12)*u.day
    coo1 = ICRS(ra=[[0.], [10.], [20.]]*u.deg,
                dec=[[-30.], [30.], [60.]]*u.deg)
    coo2 = coo1.transform_to(FK5(equinox=times))
    assert coo2.shape == (3, 12)


def test_sep():
    from astropy.coordinates.builtin_frames import ICRS

    i1 = ICRS(ra=0*u.deg, dec=1*u.deg)
    i2 = ICRS(ra=0*u.deg, dec=2*u.deg)

    sep = i1.separation(i2)
    assert sep.deg == 1

    i3 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[5, 6]*u.kpc)
    i4 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[4, 5]*u.kpc)

    sep3d = i3.separation_3d(i4)
    assert_allclose(sep3d.to(u.kpc), np.array([1, 1])*u.kpc)

    # check that it works even with velocities
    i5 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[5, 6]*u.kpc,
              pm_ra_cosdec=[1, 2]*u.mas/u.yr, pm_dec=[3, 4]*u.mas/u.yr,
              radial_velocity=[5, 6]*u.km/u.s)
    i6 = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[7, 8]*u.kpc,
              pm_ra_cosdec=[1, 2]*u.mas/u.yr, pm_dec=[3, 4]*u.mas/u.yr,
              radial_velocity=[5, 6]*u.km/u.s)

    sep3d = i5.separation_3d(i6)
    assert_allclose(sep3d.to(u.kpc), np.array([2, 2])*u.kpc)

def test_time_inputs():
    """
    Test validation and conversion of inputs for equinox and obstime attributes.
    """
    from astropy.time import Time
    from astropy.coordinates.builtin_frames import FK4

    c = FK4(1 * u.deg, 2 * u.deg, equinox='J2001.5', obstime='2000-01-01 12:00:00')
    assert c.equinox == Time('J2001.5')
    assert c.obstime == Time('2000-01-01 12:00:00')

    with pytest.raises(ValueError) as err:
        c = FK4(1 * u.deg, 2 * u.deg, equinox=1.5)
    assert 'Invalid time input' in str(err)

    with pytest.raises(ValueError) as err:
        c = FK4(1 * u.deg, 2 * u.deg, obstime='hello')
    assert 'Invalid time input' in str(err)

    # A vector time should work if the shapes match, but we don't automatically
    # broadcast the basic data (just like time).
    FK4([1, 2] * u.deg, [2, 3] * u.deg, obstime=['J2000', 'J2001'])
    with pytest.raises(ValueError) as err:
        FK4(1 * u.deg, 2 * u.deg, obstime=['J2000', 'J2001'])
    assert 'shape' in str(err)


def test_is_frame_attr_default():
    """
    Check that the `is_frame_attr_default` machinery works as expected
    """
    from astropy.time import Time
    from astropy.coordinates.builtin_frames import FK5

    c1 = FK5(ra=1*u.deg, dec=1*u.deg)
    c2 = FK5(ra=1*u.deg, dec=1*u.deg, equinox=FK5.get_frame_attr_names()['equinox'])
    c3 = FK5(ra=1*u.deg, dec=1*u.deg, equinox=Time('J2001.5'))

    assert c1.equinox == c2.equinox
    assert c1.equinox != c3.equinox

    assert c1.is_frame_attr_default('equinox')
    assert not c2.is_frame_attr_default('equinox')
    assert not c3.is_frame_attr_default('equinox')

    c4 = c1.realize_frame(r.UnitSphericalRepresentation(3*u.deg, 4*u.deg))
    c5 = c2.realize_frame(r.UnitSphericalRepresentation(3*u.deg, 4*u.deg))

    assert c4.is_frame_attr_default('equinox')
    assert not c5.is_frame_attr_default('equinox')


def test_altaz_attributes():
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, AltAz

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
    from astropy.coordinates.builtin_frames import ICRS

    # Create the frame object.
    icrs = ICRS(ra=1*u.deg, dec=1*u.deg)
    data = icrs.data

    # Create some representation objects.
    icrs_cart = icrs.cartesian
    icrs_spher = icrs.spherical

    # Testing when `_representation` set to `CartesianRepresentation`.
    icrs.representation_type = r.CartesianRepresentation

    assert icrs.representation_type == r.CartesianRepresentation
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
    icrs.representation_type = r.CylindricalRepresentation

    assert icrs.representation_type == r.CylindricalRepresentation
    assert icrs.data == data

    # Testing setter input using text argument for spherical.
    icrs.representation_type = 'spherical'

    assert icrs.representation_type is r.SphericalRepresentation
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
    icrs.representation_type = 'cylindrical'

    assert icrs.representation_type is r.CylindricalRepresentation
    assert icrs.data == data

    with pytest.raises(ValueError) as err:
        icrs.representation_type = 'WRONG'
    assert 'but must be a BaseRepresentation class' in str(err)

    with pytest.raises(ValueError) as err:
        icrs.representation_type = ICRS
    assert 'but must be a BaseRepresentation class' in str(err)


def test_represent_as():
    from astropy.coordinates.builtin_frames import ICRS

    icrs = ICRS(ra=1*u.deg, dec=1*u.deg)

    cart1 = icrs.represent_as('cartesian')
    cart2 = icrs.represent_as(r.CartesianRepresentation)

    cart1.x == cart2.x
    cart1.y == cart2.y
    cart1.z == cart2.z

    # now try with velocities
    icrs = ICRS(ra=0*u.deg, dec=0*u.deg, distance=10*u.kpc,
                pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
                radial_velocity=1*u.km/u.s)

    # single string
    rep2 = icrs.represent_as('cylindrical')
    assert isinstance(rep2, r.CylindricalRepresentation)
    assert isinstance(rep2.differentials['s'], r.CylindricalDifferential)

    # single class with positional in_frame_units, verify that warning raised
    with catch_warnings() as w:
        icrs.represent_as(r.CylindricalRepresentation, False)
        assert len(w) == 1
        assert w[0].category == AstropyWarning
        assert 'argument position' in str(w[0].message)

    # TODO: this should probably fail in the future once we figure out a better
    # workaround for dealing with UnitSphericalRepresentation's with
    # RadialDifferential's
    # two classes
    # rep2 = icrs.represent_as(r.CartesianRepresentation,
    #                          r.SphericalCosLatDifferential)
    # assert isinstance(rep2, r.CartesianRepresentation)
    # assert isinstance(rep2.differentials['s'], r.SphericalCosLatDifferential)

    with pytest.raises(ValueError):
        icrs.represent_as('odaigahara')


def test_shorthand_representations():
    from astropy.coordinates.builtin_frames import ICRS

    rep = r.CartesianRepresentation([1, 2, 3]*u.pc)
    dif = r.CartesianDifferential([1, 2, 3]*u.km/u.s)
    rep = rep.with_differentials(dif)

    icrs = ICRS(rep)

    sph = icrs.spherical
    assert isinstance(sph, r.SphericalRepresentation)
    assert isinstance(sph.differentials['s'], r.SphericalDifferential)

    sph = icrs.sphericalcoslat
    assert isinstance(sph, r.SphericalRepresentation)
    assert isinstance(sph.differentials['s'], r.SphericalCosLatDifferential)


def test_dynamic_attrs():
    from astropy.coordinates.builtin_frames import ICRS
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


def test_nodata_error():
    from astropy.coordinates.builtin_frames import ICRS

    i = ICRS()
    with pytest.raises(ValueError) as excinfo:
        i.data

    assert 'does not have associated data' in str(excinfo.value)


def test_len0_data():
    from astropy.coordinates.builtin_frames import ICRS

    i = ICRS([]*u.deg, []*u.deg)
    assert i.has_data
    repr(i)


def test_quantity_attributes():
    from astropy.coordinates.builtin_frames import GCRS

    # make sure we can create a GCRS frame with valid inputs
    GCRS(obstime='J2002', obsgeoloc=[1, 2, 3]*u.km, obsgeovel=[4, 5, 6]*u.km/u.s)

    # make sure it fails for invalid lovs or vels
    with pytest.raises(TypeError):
        GCRS(obsgeoloc=[1, 2, 3])  # no unit
    with pytest.raises(u.UnitsError):
        GCRS(obsgeoloc=[1, 2, 3]*u.km/u.s)  # incorrect unit
    with pytest.raises(ValueError):
        GCRS(obsgeoloc=[1, 3]*u.km)  # incorrect shape


@pytest.mark.remote_data
def test_eloc_attributes():
    from astropy.coordinates import AltAz, ITRS, GCRS, EarthLocation

    el = EarthLocation(lon=12.3*u.deg, lat=45.6*u.deg, height=1*u.km)
    it = ITRS(r.SphericalRepresentation(lon=12.3*u.deg, lat=45.6*u.deg, distance=1*u.km))
    gc = GCRS(ra=12.3*u.deg, dec=45.6*u.deg, distance=6375*u.km)

    el1 = AltAz(location=el).location
    assert isinstance(el1, EarthLocation)
    # these should match *exactly* because the EarthLocation
    assert el1.lat == el.lat
    assert el1.lon == el.lon
    assert el1.height == el.height

    el2 = AltAz(location=it).location
    assert isinstance(el2, EarthLocation)
    # these should *not* match because giving something in Spherical ITRS is
    # *not* the same as giving it as an EarthLocation: EarthLocation is on an
    # elliptical geoid. So the longitude should match (because flattening is
    # only along the z-axis), but latitude should not. Also, height is relative
    # to the *surface* in EarthLocation, but the ITRS distance is relative to
    # the center of the Earth
    assert not allclose(el2.lat, it.spherical.lat)
    assert allclose(el2.lon, it.spherical.lon)
    assert el2.height < -6000*u.km

    el3 = AltAz(location=gc).location
    # GCRS inputs implicitly get transformed to ITRS and then onto
    # EarthLocation's elliptical geoid. So both lat and lon shouldn't match
    assert isinstance(el3, EarthLocation)
    assert not allclose(el3.lat, gc.dec)
    assert not allclose(el3.lon, gc.ra)
    assert np.abs(el3.height) < 500*u.km


def test_equivalent_frames():
    from astropy.coordinates import SkyCoord
    from astropy.coordinates.builtin_frames import ICRS, FK4, FK5, AltAz

    i = ICRS()
    i2 = ICRS(1*u.deg, 2*u.deg)
    assert i.is_equivalent_frame(i)
    assert i.is_equivalent_frame(i2)
    with pytest.raises(TypeError):
        assert i.is_equivalent_frame(10)
    with pytest.raises(TypeError):
        assert i2.is_equivalent_frame(SkyCoord(i2))

    f1 = FK5()
    f2 = FK5(1*u.deg, 2*u.deg, equinox='J2000')
    f3 = FK5(equinox='J2010')
    f4 = FK4(equinox='J2010')

    assert f1.is_equivalent_frame(f1)
    assert not i.is_equivalent_frame(f1)
    assert f1.is_equivalent_frame(f2)
    assert not f1.is_equivalent_frame(f3)
    assert not f3.is_equivalent_frame(f4)

    aa1 = AltAz()
    aa2 = AltAz(obstime='J2010')

    assert aa2.is_equivalent_frame(aa2)
    assert not aa1.is_equivalent_frame(i)
    assert not aa1.is_equivalent_frame(aa2)


def test_representation_subclass():

    # Regression test for #3354

    from astropy.coordinates.builtin_frames import FK5

    # Normally when instantiating a frame without a distance the frame will try
    # and use UnitSphericalRepresentation internally instead of
    # SphericalRepresentation.
    frame = FK5(representation_type=r.SphericalRepresentation, ra=32 * u.deg, dec=20 * u.deg)
    assert type(frame._data) == r.UnitSphericalRepresentation
    assert frame.representation_type == r.SphericalRepresentation

    # If using a SphericalRepresentation class this used to not work, so we
    # test here that this is now fixed.
    class NewSphericalRepresentation(r.SphericalRepresentation):
        attr_classes = r.SphericalRepresentation.attr_classes

    frame = FK5(representation_type=NewSphericalRepresentation, lon=32 * u.deg, lat=20 * u.deg)
    assert type(frame._data) == r.UnitSphericalRepresentation
    assert frame.representation_type == NewSphericalRepresentation

    # A similar issue then happened in __repr__ with subclasses of
    # SphericalRepresentation.
    assert repr(frame) == ("<FK5 Coordinate (equinox=J2000.000): (lon, lat) in deg\n"
                           "    ({})>").format(' 32.,  20.' if NUMPY_LT_1_14
                                               else '32., 20.')

    # A more subtle issue is when specifying a custom
    # UnitSphericalRepresentation subclass for the data and
    # SphericalRepresentation or a subclass for the representation.

    class NewUnitSphericalRepresentation(r.UnitSphericalRepresentation):
        attr_classes = r.UnitSphericalRepresentation.attr_classes

        def __repr__(self):
            return "<NewUnitSphericalRepresentation: spam spam spam>"

    frame = FK5(NewUnitSphericalRepresentation(lon=32 * u.deg, lat=20 * u.deg),
                representation_type=NewSphericalRepresentation)

    assert repr(frame) == "<FK5 Coordinate (equinox=J2000.000):  spam spam spam>"


def test_getitem_representation():
    """
    Make sure current representation survives __getitem__ even if different
    from data representation.
    """
    from astropy.coordinates.builtin_frames import ICRS
    c = ICRS([1, 1] * u.deg, [2, 2] * u.deg)
    c.representation_type = 'cartesian'
    assert c[0].representation_type is r.CartesianRepresentation


def test_component_error_useful():
    """
    Check that a data-less frame gives useful error messages about not having
    data when the attributes asked for are possible coordinate components
    """
    from astropy.coordinates.builtin_frames import ICRS

    i = ICRS()

    with pytest.raises(ValueError) as excinfo:
        i.ra
    assert 'does not have associated data' in str(excinfo.value)

    with pytest.raises(AttributeError) as excinfo1:
        i.foobar
    with pytest.raises(AttributeError) as excinfo2:
        i.lon  # lon is *not* the component name despite being the underlying representation's name
    assert "object has no attribute 'foobar'" in str(excinfo1.value)
    assert "object has no attribute 'lon'" in str(excinfo2.value)


def test_cache_clear():
    from astropy.coordinates.builtin_frames import ICRS

    i = ICRS(1*u.deg, 2*u.deg)

    # Add an in frame units version of the rep to the cache.
    repr(i)

    assert len(i.cache['representation']) == 2

    i.cache.clear()

    assert len(i.cache['representation']) == 0


def test_inplace_array():
    from astropy.coordinates.builtin_frames import ICRS

    i = ICRS([[1, 2], [3, 4]]*u.deg, [[10, 20], [30, 40]]*u.deg)

    # Add an in frame units version of the rep to the cache.
    repr(i)

    # Check that repr() has added a rep to the cache
    assert len(i.cache['representation']) == 2

    # Modify the data
    i.data.lon[:, 0] = [100, 200]*u.deg

    # Clear the cache
    i.cache.clear()

    # This will use a second (potentially cached rep)
    assert_allclose(i.ra, [[100, 2], [200, 4]]*u.deg)
    assert_allclose(i.dec, [[10, 20], [30, 40]]*u.deg)


def test_inplace_change():
    from astropy.coordinates.builtin_frames import ICRS

    i = ICRS(1*u.deg, 2*u.deg)

    # Add an in frame units version of the rep to the cache.
    repr(i)

    # Check that repr() has added a rep to the cache
    assert len(i.cache['representation']) == 2

    # Modify the data
    i.data.lon[()] = 10*u.deg

    # Clear the cache
    i.cache.clear()

    # This will use a second (potentially cached rep)
    assert i.ra == 10 * u.deg
    assert i.dec == 2 * u.deg


def test_representation_with_multiple_differentials():
    from astropy.coordinates.builtin_frames import ICRS

    dif1 = r.CartesianDifferential([1, 2, 3]*u.km/u.s)
    dif2 = r.CartesianDifferential([1, 2, 3]*u.km/u.s**2)
    rep = r.CartesianRepresentation([1, 2, 3]*u.pc,
                                    differentials={'s': dif1, 's2': dif2})

    # check warning is raised for a scalar
    with pytest.raises(ValueError):
        ICRS(rep)


def test_representation_arg_backwards_compatibility():
    # TODO: this test can be removed when the `representation` argument is
    # removed from the BaseCoordinateFrame initializer.
    from astropy.coordinates.builtin_frames import ICRS

    c1 = ICRS(x=1*u.pc, y=2*u.pc, z=3*u.pc,
              representation_type=r.CartesianRepresentation)

    c2 = ICRS(x=1*u.pc, y=2*u.pc, z=3*u.pc,
              representation_type=r.CartesianRepresentation)

    c3 = ICRS(x=1*u.pc, y=2*u.pc, z=3*u.pc,
              representation_type='cartesian')

    assert c1.x == c2.x
    assert c1.y == c2.y
    assert c1.z == c2.z

    assert c1.x == c3.x
    assert c1.y == c3.y
    assert c1.z == c3.z

    assert c1.representation_type == c1.representation_type

    with pytest.raises(ValueError):
        ICRS(x=1*u.pc, y=2*u.pc, z=3*u.pc,
             representation_type='cartesian',
             representation='cartesian')


def test_missing_component_error_names():
    """
    This test checks that the component names are frame component names, not
    representation or differential names, when referenced in an exception raised
    when not passing in enough data. For example:

    ICRS(ra=10*u.deg)

    should state:

    TypeError: __init__() missing 1 required positional argument: 'dec'
    """
    from astropy.coordinates.builtin_frames import ICRS

    with pytest.raises(TypeError) as e:
        ICRS(ra=150 * u.deg)
    assert "missing 1 required positional argument: 'dec'" in str(e)

    with pytest.raises(TypeError) as e:
        ICRS(ra=150*u.deg, dec=-11*u.deg,
             pm_ra=100*u.mas/u.yr, pm_dec=10*u.mas/u.yr)
    assert "pm_ra_cosdec" in str(e)


def test_non_spherical_representation_unit_creation(unitphysics):
    from astropy.coordinates.builtin_frames import ICRS

    class PhysicsICRS(ICRS):
        default_representation = r.PhysicsSphericalRepresentation

    pic = PhysicsICRS(phi=1*u.deg, theta=25*u.deg, r=1*u.kpc)
    assert isinstance(pic.data, r.PhysicsSphericalRepresentation)

    picu = PhysicsICRS(phi=1*u.deg, theta=25*u.deg)
    assert isinstance(picu.data, unitphysics)


def test_attribute_repr():
    from astropy.coordinates.attributes import Attribute
    from astropy.coordinates.baseframe import BaseCoordinateFrame

    class Spam:
        def _astropy_repr_in_frame(self):
            return "TEST REPR"

    class TestFrame(BaseCoordinateFrame):
        attrtest = Attribute(default=Spam())

    assert "TEST REPR" in repr(TestFrame())
