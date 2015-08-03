# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Accuracy tests for GCRS coordinate transformations, primarily to/from AltAz.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ....tests.helper import pytest
from .... import units as u
from ....time import Time
from ...builtin_frames import AltAz
from ... import EarthLocation, get_sun, GCRS, CIRS, CartesianRepresentation


# shared by parametrized tests
std_altazs = [AltAz(location=EarthLocation(-90*u.deg, 65*u.deg),
                    obstime=Time('J2000')), #J2000 is often a default so this might work when others don't
              AltAz(location=EarthLocation(-90*u.deg, 65*u.deg),
                    obstime=Time('2014-01-01 00:00:00')),
              AltAz(location=EarthLocation(-90*u.deg, 65*u.deg),
                    obstime=Time('2014-08-01 08:00:00')),
              AltAz(location=EarthLocation(120*u.deg, -35*u.deg),
                    obstime=Time('2014-01-01 00:00:00'))
              ]


@pytest.mark.parametrize('altaz', std_altazs)
def test_gcrs_altaz_sun(altaz):
    """
    Sanity-check that the sun is at a reasonable distance from any altaz
    """
    earthecc = 0.017  # roughly earth orbital eccentricity

    sun = get_sun(altaz.obstime)
    assert sun.frame.name == 'gcrs'

    # the .to(u.au) is not necessary, it just makes the asserts on failure more readable
    assert sun.distance.to(u.au) < (earthecc + 1)*u.au
    assert sun.distance.to(u.au) > (earthecc - 1)*u.au

    sunaa = sun.transform_to(altaz)
    assert sunaa.distance.to(u.au) < (earthecc + 1)*u.au
    assert sunaa.distance.to(u.au) > (earthecc - 1)*u.au


@pytest.mark.parametrize('altaz', std_altazs)
def test_gcrs_altaz_moonish(altaz):
    """
    Sanity-check that an object resembling the moon goes to the right place
    """
    moondist = 385000*u.km  # approximate moon semi-major orbit axis of moon
    cart = CartesianRepresentation(3**-0.5*moondist, 3**-0.5*moondist, 3**-0.5*moondist)
    moon = GCRS(cart, obstime=altaz.obstime)

    moonaa = moon.transform_to(altaz)
    assert np.abs(moonaa.distance - moon.distance).to(u.km) > 1000*u.km
    assert np.abs(moonaa.distance - moon.distance).to(u.km) < 7000*u.km

    #also should add checks that the alt/az are different for different earth locations

@pytest.mark.parametrize('altaz', std_altazs)
def test_cirs_altaz_moonish(altaz):
    """
    Sanity-check that an object resembling the moon goes to the right place if starting from CIRS
    """
    moondist = 385000*u.km  # approximate moon semi-major orbit axis of moon
    cart = CartesianRepresentation(3**-0.5*moondist, 3**-0.5*moondist, 3**-0.5*moondist)
    moon = CIRS(cart, obstime=altaz.obstime)

    moonaa = moon.transform_to(altaz)
    assert np.abs(moonaa.distance - moon.distance).to(u.km) > 1000*u.km
    assert np.abs(moonaa.distance - moon.distance).to(u.km) < 7000*u.km
