# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Check `astropy.coordinates` against the `pytpm` example results from here:
http://phn.github.io/pytpm/conversions.html
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .... import units as u
from ....tests.helper import pytest
from ....time import Time
from ....coordinates import SkyCoord, Angle, EarthLocation, AltAz


def test_fk5_equinox_and_epoch_j2000_0_to_fk4_equinox_and_epoch_b1950_0():
    """
    http://phn.github.io/pytpm/conversions.html#fk5-equinox-and-epoch-j2000-0-to-fk4-equinox-and-epoch-b1950-0
    """
    pos = SkyCoord('12h22m54.899s 15d49m20.57s', frame='fk5')
    actual = pos.transform_to('fk4')
    expected = SkyCoord('12h20m22.935s 16d05m58.024s', frame='fk4')
    distance = actual.separation(expected)
    # TODO: what assertion precision makes sense here?
    assert distance < 1 * u.arcsec


def test_fk5_equinox_and_epoch_j2000_0_to_iau_1958_galactic_system():
    """
    http://phn.github.io/pytpm/conversions.html#fk5-equinox-and-epoch-j2000-0-to-iau-1958-galactic-system
    """
    # TODO: what assertion precision makes sense here?
    assert True

@pytest.mark.xfail
def test_fk5_equinox_and_epoch_j2000_0_to_topocentric_observed():
    """
    http://phn.github.io/pytpm/conversions.html#fk5-equinox-and-epoch-j2000-0-to-topocentric-observed
    """
    obstime = Time('2010-01-01 12:00:00', scale='utc')
    # Observatory position for `kpno` from here:
    # http://idlastro.gsfc.nasa.gov/ftp/pro/astro/observatory.pro
    location = EarthLocation(lon=Angle('111d36.0m'), lat=Angle('31d57.8m'),height=2120.*u.m)
    temperature = 0 * u.deg_C
    pressure = 0.781 * u.bar
    # relative_humidity = ?
    # obswl = ?
    aaframe = AltAz(obstime=obstime,location=location,
                    temperature=temperature, pressure=pressure)

    eq = SkyCoord('12h22m54.899s 15d49m20.57s', frame='fk5')
    hor_actual = eq.transform_to(aaframe)
    hor_expected = SkyCoord('264d55m06s 37d54m41s', frame='altaz')
    distance = hor_actual.separation(hor_expected)

    assert distance < 1 * u.arcsec
