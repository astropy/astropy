# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree

import astropy.units as u
from astropy import units
from astropy.coordinates import ICRS, HCRS, Longitude, Latitude, Angle

from ....extension import AstropyExtension


def test_hcrs_basic(tmpdir):
    ra = Longitude(25, unit=units.deg)
    dec = Latitude(45, unit=units.deg)

    tree = {'coord': HCRS(ra=ra, dec=dec)}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())


def test_icrs_basic(tmpdir):
    wrap_angle = Angle(1.5, unit=units.rad)
    ra = Longitude(25, unit=units.deg, wrap_angle=wrap_angle)
    dec = Latitude(45, unit=units.deg)

    tree = {'coord': ICRS(ra=ra, dec=dec)}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())


def test_icrs_nodata(tmpdir):
    tree = {'coord': ICRS()}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())


def test_hcrs_compound(tmpdir):

    icrs = ICRS(ra=[0, 1, 2]*units.deg, dec=[3, 4, 5]*units.deg)

    tree = {'coord': icrs}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())


def test_sunpy(tmpdir):

    scoord = pytest.importorskip('sunpy.coordinates')

    hpc = scoord.Helioprojective(100*units.arcsec, 100*units.arcsec,
                                 obstime="2011-01-01")

    tree = {'coord': hpc}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())


def test_sunpy_hgs(tmpdir):

    scoord = pytest.importorskip('sunpy.coordinates')

    hpc = scoord.HeliographicStonyhurst(1*u.deg, 0*u.deg, 0.98*u.AU)

    tree = {'coord': hpc}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())
