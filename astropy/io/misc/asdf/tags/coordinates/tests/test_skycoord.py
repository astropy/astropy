# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import numpy as np

import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic, FK4, FK5, Longitude

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree


# These tests are cribbed directly from the Examples section of
# https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html


def test_scalar_skycoord(tmpdir):

    c = SkyCoord(10, 20, unit="deg")  # defaults to ICRS frame
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


def test_vector_skycoord(tmpdir):

    c = SkyCoord([1, 2, 3], [-30, 45, 8], frame="icrs", unit="deg")  # 3 coords
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


def test_skycoord_fk4(tmpdir):

    coords = ["1:12:43.2 +1:12:43", "1 12 43.2 +1 12 43"]
    c = SkyCoord(coords, frame=FK4, unit=(u.deg, u.hourangle), obstime="J1992.21")
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.parametrize('coord', [
    SkyCoord("1h12m43.2s +1d12m43s", frame=Galactic),  # Units from string
    SkyCoord(frame="galactic", l="1h12m43.2s", b="+1d12m43s")
])
def test_skycoord_galactic(coord, tmpdir):

    tree = dict(coord=coord)
    assert_roundtrip_tree(tree, tmpdir)


def test_skycoord_ra_dec(tmpdir):

    ra = Longitude([1, 2, 3], unit=u.deg)  # Could also use Angle
    dec = np.array([4.5, 5.2, 6.3]) * u.deg  # Astropy Quantity
    c = SkyCoord(ra, dec, frame='icrs')
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)

    c = SkyCoord(frame=ICRS, ra=ra, dec=dec, obstime='2001-01-02T12:34:56')
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


def test_skycoord_override_defaults(tmpdir):

    c = FK4(1 * u.deg, 2 * u.deg)  # Uses defaults for obstime, equinox
    c = SkyCoord(c, obstime='J2010.11', equinox='B1965')  # Override defaults
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


def test_skycoord_cartesian(tmpdir):

    c = SkyCoord(w=0, u=1, v=2, unit='kpc', frame='galactic',
                   representation_type='cartesian')
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


def test_skycoord_vector_frames(tmpdir):

    c = SkyCoord([ICRS(ra=1*u.deg, dec=2*u.deg), ICRS(ra=3*u.deg, dec=4*u.deg)])
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.xfail(reason='Velocities are not properly serialized yet')
def test_skycoord_radial_velocity(tmpdir):

    c = SkyCoord(ra=1*u.deg, dec=2*u.deg, radial_velocity=10*u.km/u.s)
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.xfail(reason='Velocities are not properly serialized yet')
def test_skycoord_proper_motion(tmpdir):

    c = SkyCoord(ra=1*u.deg, dec=2*u.deg, pm_ra_cosdec=2*u.mas/u.yr,
                 pm_dec=1*u.mas/u.yr)
    tree = dict(coord=c)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.skip(reason='Apparent loss of precision during serialization')
def test_skycoord_extra_attribute(tmpdir):

    sc = SkyCoord(10*u.deg, 20*u.deg, equinox="2011-01-01T00:00", frame="fk4")
    tree = dict(coord=sc.transform_to("icrs"))

    def check_asdf(asdffile):
        assert hasattr(asdffile['coord'], 'equinox')

    assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check_asdf)


def test_skycoord_2d_obstime(tmpdir):

    sc = SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,m', frame='fk4',
                    obstime=['J1990.5', 'J1991.5']),
    tree = dict(coord=sc)
    assert_roundtrip_tree(tree, tmpdir)
