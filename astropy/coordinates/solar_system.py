# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains convenience functions for retrieving solar system
ephemerides from jplephem.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from .sky_coordinate import SkyCoord
from ..utils.data import download_file
from .. import units as u
from ..constants import c as speed_of_light
from .representation import CartesianRepresentation
from .builtin_frames import GCRS
from .earth import EarthLocation


__all__ = ["get_planet", "get_moon"]

KERNEL = None


def _download_spk_file(show_progress=True):
    """
    Get the Satellite Planet Kernel (SPK) file `de430.bsp` from NASA JPL.

    Download the file from the JPL webpage once and subsequently access a
    cached copy. This file is ~120 MB, and covers years ~1550-2650 CE [1]_.

    .. [1] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de430-de431.txt
    """
    de430_url = ('http://naif.jpl.nasa.gov/pub/naif/'
                 'generic_kernels/spk/planets/de430.bsp')
    return download_file(de430_url, cache=True, show_progress=show_progress)


def _get_kernel(*args, **kwargs):
    """
    Try importing jplephem, download/retrieve from cahce the Satellite Planet
    Kernel.
    """
    global KERNEL

    try:
        from jplephem.spk import SPK
    except ImportError:
        raise ImportError("Solar system ephemeris calculations depend on"
                          "jplephem")

    if KERNEL is None:
        KERNEL = SPK.open(_download_spk_file())
    return KERNEL


def _get_absolute_planet_position(time, planet_index):
    """
    Calculate the position of planet ``planet_index`` in cartesian coordinates.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    cartesian_position : `~astropy.units.Quantity`
        GRCS position of the planet in cartesian coordinates

    earth_distance : `~astropy.units.Quantity`
        Distance between Earth and planet.

    Notes
    -----

    """
    kernel = _get_kernel()

    # Mercury and Venus have separate barycenter and planet vectors
    if planet_index < 3:
        barycenter_to_planet_ind = 100*planet_index + 99
        cartesian_position_planet = (kernel[0, planet_index].compute(time.jd) +
                                     kernel[planet_index,
                                            barycenter_to_planet_ind].compute(time.jd))
    # planet_index is 3 gets the moon
    elif planet_index == 3:
        barycenter_to_planet_ind = 301
        cartesian_position_planet = (kernel[0, planet_index].compute(time.jd) +
                                     kernel[planet_index,
                                            barycenter_to_planet_ind].compute(time.jd))
    else:
        cartesian_position_planet = kernel[0, planet_index].compute(time.jd)

    cartesian_position_earth = (kernel[0, 3].compute(time.jd) +
                                kernel[3, 399].compute(time.jd))

    # Quadrature sum of the difference of the position vectors
    earth_distance = np.sqrt(np.sum((np.array(cartesian_position_planet) -
                                     np.array(cartesian_position_earth))**2))*u.km

    earth_to_planet_vector = u.Quantity(cartesian_position_planet -
                                        cartesian_position_earth,
                                        unit=u.km)

    return earth_to_planet_vector, earth_distance


def _get_apparent_planet_position(time, planet_index):
    """
    Calculate the apparent position of planet ``planet_index`` in cartesian
    coordinates, given the approximate light travel time to the object.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Position of the planet defined as a vector from the Earth.
    """
    # Get distance of planet at `time`
    earth_to_planet_vector, earth_distance = _get_absolute_planet_position(time, planet_index)

    # The apparent position depends on the time that the light was emitted from
    # the distant planet, so subtract off the light travel time
    light_travel_time = earth_distance/speed_of_light
    emitted_time = time - light_travel_time

    # Calculate position given approximate light travel time.
    # TODO: this should be solved iteratively to converge on precise positions
    earth_to_planet_vector, earth_distance = _get_absolute_planet_position(emitted_time, planet_index)
    x, y, z = earth_to_planet_vector

    return CartesianRepresentation(x=x, y=y, z=z)


def get_planet(time, planet_index, location=None):
    """
    Get a `~astropy.coordinates.SkyCoord` for a planet as observed from a
    location on Earth.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto), excluding the
        special value 3, which corresponds to the Earth's Moon.

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, set to
        Greenwich.

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the planet
    """
    if location is None:
        location = EarthLocation.of_site('greenwich')

    cartrep = _get_apparent_planet_position(time, planet_index)

    return SkyCoord(GCRS(cartrep, obstime=time,
                         obsgeoloc=u.Quantity(location.geocentric, copy=False)))


def get_moon(time, location=None):
    """
    Get a `~astropy.coordinates.SkyCoord` for the Earth's moon as observed
    from a location on Earth.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, set to
        Greenwich.

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the planet
    """
    return get_planet(time, planet_index=3, location=location)
