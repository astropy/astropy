# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import pytest
import numpy as np

from ... import units as u
from ..builtin_frames import ICRS, Galactic, Galactocentric
from .. import builtin_frames as bf
from ...units import allclose as quantity_allclose
from ..errors import ConvertError
from .. import representation as r

def test_api():
    # transform observed Barycentric velocities to full-space Galactocentric
    gc_frame = Galactocentric()
    icrs = ICRS(ra=151.*u.deg, dec=-16*u.deg, distance=101*u.pc,
                pm_ra_cosdec=21*u.mas/u.yr, pm_dec=-71*u.mas/u.yr,
                radial_velocity=71*u.km/u.s)
    icrs.transform_to(gc_frame)

    # transform a set of ICRS proper motions to Galactic
    icrs = ICRS(ra=151.*u.deg, dec=-16*u.deg,
                pm_ra_cosdec=21*u.mas/u.yr, pm_dec=-71*u.mas/u.yr)
    icrs.transform_to(Galactic)

    # transform a Barycentric RV to a GSR RV
    icrs = ICRS(ra=151.*u.deg, dec=-16*u.deg, distance=1.*u.pc,
                pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
                radial_velocity=71*u.km/u.s)
    icrs.transform_to(Galactocentric)

all_kwargs = [
    dict(ra=37.4*u.deg, dec=-55.8*u.deg),
    dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc),
    dict(ra=37.4*u.deg, dec=-55.8*u.deg,
         pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr),
    dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
         pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr),
    dict(ra=37.4*u.deg, dec=-55.8*u.deg,
         radial_velocity=105.7*u.km/u.s),
    dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
         radial_velocity=105.7*u.km/u.s),
    dict(ra=37.4*u.deg, dec=-55.8*u.deg,
         radial_velocity=105.7*u.km/u.s,
         pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr),
    dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
         pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr,
         radial_velocity=105.7*u.km/u.s),
    # Now test other representation/differential types:
    dict(x=100.*u.pc, y=200*u.pc, z=300*u.pc,
         representation_type='cartesian'),
    dict(x=100.*u.pc, y=200*u.pc, z=300*u.pc,
         representation_type=r.CartesianRepresentation),
    dict(x=100.*u.pc, y=200*u.pc, z=300*u.pc,
         v_x=100.*u.km/u.s, v_y=200*u.km/u.s, v_z=300*u.km/u.s,
         representation_type=r.CartesianRepresentation,
         differential_type=r.CartesianDifferential),
    dict(x=100.*u.pc, y=200*u.pc, z=300*u.pc,
         v_x=100.*u.km/u.s, v_y=200*u.km/u.s, v_z=300*u.km/u.s,
         representation_type=r.CartesianRepresentation,
         differential_type='cartesian'),
]

@pytest.mark.parametrize('kwargs', all_kwargs)
def test_all_arg_options(kwargs):
    # Above is a list of all possible valid combinations of arguments.
    # Here we do a simple thing and just verify that passing them in, we have
    # access to the relevant attributes from the resulting object
    icrs = ICRS(**kwargs)
    gal = icrs.transform_to(Galactic)
    repr_gal = repr(gal)

    for k in kwargs:
        if k == 'differential_type':
            continue
        getattr(icrs, k)

    if 'pm_ra_cosdec' in kwargs: # should have both
        assert 'pm_l_cosb' in repr_gal
        assert 'pm_b' in repr_gal
        assert 'mas / yr' in repr_gal

        if 'radial_velocity' not in kwargs:
            assert 'radial_velocity' not in repr_gal

    if 'radial_velocity' in kwargs:
        assert 'radial_velocity' in repr_gal
        assert 'km / s' in repr_gal

        if 'pm_ra_cosdec' not in kwargs:
            assert 'pm_l_cosb' not in repr_gal
            assert 'pm_b' not in repr_gal

@pytest.mark.parametrize('cls,lon,lat', [
    [bf.ICRS, 'ra', 'dec'], [bf.FK4, 'ra', 'dec'], [bf.FK4NoETerms, 'ra', 'dec'],
    [bf.FK5, 'ra', 'dec'], [bf.GCRS, 'ra', 'dec'], [bf.HCRS, 'ra', 'dec'],
    [bf.LSR, 'ra', 'dec'], [bf.CIRS, 'ra', 'dec'], [bf.Galactic, 'l', 'b'],
    [bf.AltAz, 'az', 'alt'], [bf.Supergalactic, 'sgl', 'sgb'],
    [bf.GalacticLSR, 'l', 'b'], [bf.HeliocentricTrueEcliptic, 'lon', 'lat'],
    [bf.GeocentricTrueEcliptic, 'lon', 'lat'],
    [bf.BarycentricTrueEcliptic, 'lon', 'lat'],
    [bf.PrecessedGeocentric, 'ra', 'dec']
])
def test_expected_arg_names(cls, lon, lat):
    kwargs = {lon: 37.4*u.deg, lat: -55.8*u.deg, 'distance': 150*u.pc,
              'pm_{0}_cos{1}'.format(lon, lat): -21.2*u.mas/u.yr,
              'pm_{0}'.format(lat): 17.1*u.mas/u.yr,
              'radial_velocity': 105.7*u.km/u.s}
    frame = cls(**kwargs)


# these data are extracted from the vizier copy of XHIP:
# http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=+V/137A/XHIP
_xhip_head = """
------ ------------ ------------ -------- -------- ------------ ------------ ------- -------- -------- ------- ------ ------ ------
       R            D            pmRA     pmDE                               Di      pmGLon   pmGLat   RV      U      V      W
HIP    AJ2000 (deg) EJ2000 (deg) (mas/yr) (mas/yr) GLon (deg)   GLat (deg)   st (pc) (mas/yr) (mas/yr) (km/s)  (km/s) (km/s) (km/s)
------ ------------ ------------ -------- -------- ------------ ------------ ------- -------- -------- ------- ------ ------ ------
"""[1:-1]
_xhip_data = """
    19 000.05331690 +38.30408633    -3.17   -15.37 112.00026470 -23.47789171  247.12    -6.40   -14.33    6.30    7.3    2.0  -17.9
    20 000.06295067 +23.52928427    36.11   -22.48 108.02779304 -37.85659811   95.90    29.35   -30.78   37.80  -19.3   16.1  -34.2
    21 000.06623581 +08.00723430    61.48    -0.23 101.69697120 -52.74179515  183.68    58.06   -20.23  -11.72  -45.2  -30.9   -1.3
 24917 080.09698238 -33.39874984    -4.30    13.40 236.92324669 -32.58047131  107.38   -14.03    -1.15   36.10  -22.4  -21.3  -19.9
 59207 182.13915108 +65.34963517    18.17     5.49 130.04157185  51.18258601   56.00   -18.98    -0.49    5.70    1.5    6.1    4.4
 87992 269.60730667 +36.87462906   -89.58    72.46  62.98053142  25.90148234  129.60    45.64   105.79   -4.00  -39.5  -15.8   56.7
115110 349.72322473 -28.74087144    48.86    -9.25  23.00447250 -69.52799804  116.87    -8.37   -49.02   15.00  -16.8  -12.2  -23.6
"""[1:-1]

# in principal we could parse the above as a table, but doing it "manually"
# makes this test less tied to Table working correctly


@pytest.mark.parametrize('hip,ra,dec,pmra,pmdec,glon,glat,dist,pmglon,pmglat,rv,U,V,W',
                         [[float(val) for val in row.split()] for row in _xhip_data.split('\n')])
def test_xhip_galactic(hip, ra, dec, pmra, pmdec, glon, glat, dist, pmglon, pmglat, rv, U, V, W):
    i = ICRS(ra*u.deg, dec*u.deg, dist*u.pc,
             pm_ra_cosdec=pmra*u.marcsec/u.yr, pm_dec=pmdec*u.marcsec/u.yr,
             radial_velocity=rv*u.km/u.s)
    g = i.transform_to(Galactic)

    # precision is limited by 2-deciimal digit string representation of pms
    assert quantity_allclose(g.pm_l_cosb, pmglon*u.marcsec/u.yr, atol=.01*u.marcsec/u.yr)
    assert quantity_allclose(g.pm_b, pmglat*u.marcsec/u.yr, atol=.01*u.marcsec/u.yr)

    # make sure UVW also makes sense
    uvwg = g.cartesian.differentials['s']
    # precision is limited by 1-decimal digit string representation of vels
    assert quantity_allclose(uvwg.d_x, U*u.km/u.s, atol=.1*u.km/u.s)
    assert quantity_allclose(uvwg.d_y, V*u.km/u.s, atol=.1*u.km/u.s)
    assert quantity_allclose(uvwg.d_z, W*u.km/u.s, atol=.1*u.km/u.s)

@pytest.mark.parametrize('kwargs,expect_success', [
    [dict(ra=37.4*u.deg, dec=-55.8*u.deg), False],
    [dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc), True],
    [dict(ra=37.4*u.deg, dec=-55.8*u.deg,
          pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr), False],
    [dict(ra=37.4*u.deg, dec=-55.8*u.deg, radial_velocity=105.7*u.km/u.s), False],
    [dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
          radial_velocity=105.7*u.km/u.s), False],
    [dict(ra=37.4*u.deg, dec=-55.8*u.deg,
          radial_velocity=105.7*u.km/u.s,
          pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr), False],
    [dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
          pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr,
          radial_velocity=105.7*u.km/u.s), True]

])
def test_frame_affinetransform(kwargs, expect_success):
    """There are already tests in test_transformations.py that check that
    an AffineTransform fails without full-space data, but this just checks that
    things work as expected at the frame level as well.
    """

    icrs = ICRS(**kwargs)

    if expect_success:
        gc = icrs.transform_to(Galactocentric)

    else:
        with pytest.raises(ConvertError):
            icrs.transform_to(Galactocentric)

def test_differential_type_arg():
    """
    Test passing in an explicit differential class to the initializer or
    changing the differential class via set_representation_cls
    """
    from ..builtin_frames import ICRS

    icrs = ICRS(ra=1*u.deg, dec=60*u.deg,
                pm_ra=10*u.mas/u.yr, pm_dec=-11*u.mas/u.yr,
                differential_type=r.UnitSphericalDifferential)
    assert icrs.pm_ra == 10*u.mas/u.yr

    icrs = ICRS(ra=1*u.deg, dec=60*u.deg,
                pm_ra=10*u.mas/u.yr, pm_dec=-11*u.mas/u.yr,
                differential_type={'s': r.UnitSphericalDifferential})
    assert icrs.pm_ra == 10*u.mas/u.yr

    icrs = ICRS(ra=1*u.deg, dec=60*u.deg,
                pm_ra_cosdec=10*u.mas/u.yr, pm_dec=-11*u.mas/u.yr)
    icrs.set_representation_cls(s=r.UnitSphericalDifferential)
    assert quantity_allclose(icrs.pm_ra, 20*u.mas/u.yr)

    # incompatible representation and differential
    with pytest.raises(TypeError):
        ICRS(ra=1*u.deg, dec=60*u.deg,
             v_x=1*u.km/u.s, v_y=-2*u.km/u.s, v_z=-2*u.km/u.s,
             differential_type=r.CartesianDifferential)

    # specify both
    icrs = ICRS(x=1*u.pc, y=2*u.pc, z=3*u.pc,
                v_x=1*u.km/u.s, v_y=2*u.km/u.s, v_z=3*u.km/u.s,
                representation_type=r.CartesianRepresentation,
                differential_type=r.CartesianDifferential)
    assert icrs.x == 1*u.pc
    assert icrs.y == 2*u.pc
    assert icrs.z == 3*u.pc
    assert icrs.v_x == 1*u.km/u.s
    assert icrs.v_y == 2*u.km/u.s
    assert icrs.v_z == 3*u.km/u.s


def test_slicing_preserves_differential():
    icrs = ICRS(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
                pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr,
                radial_velocity=105.7*u.km/u.s)
    icrs2 = icrs.reshape(1,1)[:1,0]

    for name in icrs.representation_component_names.keys():
        assert getattr(icrs, name) == getattr(icrs2, name)[0]

    for name in icrs.get_representation_component_names('s').keys():
        assert getattr(icrs, name) == getattr(icrs2, name)[0]


def test_shorthand_attributes():
    # Check that attribute access works

    # for array data:
    n = 4
    icrs1 = ICRS(ra=np.random.uniform(0, 360, n)*u.deg,
                 dec=np.random.uniform(-90, 90, n)*u.deg,
                 distance=100*u.pc,
                 pm_ra_cosdec=np.random.normal(0, 100, n)*u.mas/u.yr,
                 pm_dec=np.random.normal(0, 100, n)*u.mas/u.yr,
                 radial_velocity=np.random.normal(0, 100, n)*u.km/u.s)
    v = icrs1.velocity
    pm = icrs1.proper_motion
    assert quantity_allclose(pm[0], icrs1.pm_ra_cosdec)
    assert quantity_allclose(pm[1], icrs1.pm_dec)

    # for scalar data:
    icrs2 = ICRS(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
                 pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr,
                 radial_velocity=105.7*u.km/u.s)
    v = icrs2.velocity
    pm = icrs2.proper_motion
    assert quantity_allclose(pm[0], icrs2.pm_ra_cosdec)
    assert quantity_allclose(pm[1], icrs2.pm_dec)

    # check that it fails where we expect:

    # no distance
    rv = 105.7*u.km/u.s
    icrs3 = ICRS(ra=37.4*u.deg, dec=-55.8*u.deg,
                 pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr,
                 radial_velocity=rv)
    with pytest.raises(ValueError):
        icrs3.velocity

    icrs3.set_representation_cls('cartesian')
    assert hasattr(icrs3, 'radial_velocity')
    assert quantity_allclose(icrs3.radial_velocity, rv)

    icrs4 = ICRS(x=30*u.pc, y=20*u.pc, z=11*u.pc,
                 v_x=10*u.km/u.s, v_y=10*u.km/u.s, v_z=10*u.km/u.s,
                 representation_type=r.CartesianRepresentation,
                 differential_type=r.CartesianDifferential)
    icrs4.radial_velocity


def test_negative_distance():
    """ Regression test: #7408
    Make sure that negative parallaxes turned into distances are handled right
    """

    RA = 150 * u.deg
    DEC = -11*u.deg
    c = ICRS(ra=RA, dec=DEC,
             distance=(-10*u.mas).to(u.pc, u.parallax()),
             pm_ra_cosdec=10*u.mas/u.yr,
             pm_dec=10*u.mas/u.yr)
    assert quantity_allclose(c.ra, RA)
    assert quantity_allclose(c.dec, DEC)

    c = ICRS(ra=RA, dec=DEC,
             distance=(-10*u.mas).to(u.pc, u.parallax()))
    assert quantity_allclose(c.ra, RA)
    assert quantity_allclose(c.dec, DEC)
