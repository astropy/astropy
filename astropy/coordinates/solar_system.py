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
from .builtin_frames import GCRS, ICRS
from .builtin_frames.utils import get_jd12
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


def _get_barycentric_planet_position(time, planet_index):
    """
    Calculate the barycentric position of planet ``planet_index`` in cartesian coordinates.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto. 10 for Moon)

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) position of the planet in cartesian coordinates

    Notes
    -----

    """
    kernel = _get_kernel()

    jd1, jd2 = get_jd12(time,'tdb')
    
    # Mercury and Venus have separate barycenter and planet vectors
    if planet_index < 3:
        barycenter_to_planet_ind = 100*planet_index + 99
        cartesian_position_planet = (kernel[0, planet_index].compute(jd1, jd2) +
                                     kernel[planet_index,
                                            barycenter_to_planet_ind].compute(jd1, jd2))
    # planet_index 3 gets the Earth itself, not the barycenter of the Earth-Moon system
    elif planet_index == 3:
        barycenter_to_planet_ind = 399
        cartesian_position_planet = (kernel[0, planet_index].compute(jd1, jd2) +
                                     kernel[planet_index,
                                            barycenter_to_planet_ind].compute(jd1, jd2))
    # planet_index 10 gets the moon
    elif planet_index == 10:
        barycenter_to_planet_ind = 301
        cartesian_position_planet = (kernel[0, 3].compute(jd1, jd2) +
                                     kernel[3,
                                            barycenter_to_planet_ind].compute(jd1, jd2))
    else:
        cartesian_position_planet = kernel[0, planet_index].compute(jd1, jd2)

    barycen_to_planet_vector = u.Quantity(cartesian_position_planet, unit=u.km)    
    return CartesianRepresentation(barycen_to_planet_vector)

def _get_earth_planet_vector(time, planet_index):
    """
    Calculate the vector between the Geocenter and planet ``planet_index``. 

    This routine calculates the vector between the Earth's Geocenter and the planet
    specified by ``planet_index``. 

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto. 10 for Moon)

    Returns
    -------
    earth_planet_vector : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) vector from Geocenter to the planet in cartesian coordinates
        
    earth_distance : `~astropy.units.Quantity`
        Distance between Earth and planet.

    Notes
    -----

    """    
    earth_loc = _get_barycentric_planet_position(time, 3)    
    planet_loc = _get_barycentric_planet_position(time, planet_index)
    earth_planet_vector = planet_loc.xyz - earth_loc.xyz
    earth_distance = np.sqrt(earth_planet_vector.dot(earth_planet_vector))
    return earth_planet_vector, earth_distance
    
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
        Index of the planet (1-9 for Mercury through Pluto, 10 for Moon)

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) apparent position of the planet in cartesian coordinates
    """
    # Get distance of planet at `time`
    earth_to_planet_vector, earth_distance = _get_earth_planet_vector(time, planet_index)

    # The apparent position depends on the time that the light was emitted from
    # the distant planet, so subtract off the light travel time
    light_travel_time = earth_distance/speed_of_light
    emitted_time = time - light_travel_time

    # Calculate position given approximate light travel time.
    delta_light_travel_time = 20*u.s
    while (np.fabs(delta_light_travel_time)) > 1.0e-8*u.s:
        earth_to_planet_vector, earth_distance = _get_earth_planet_vector(emitted_time, planet_index)
        delta_light_travel_time = light_travel_time - earth_distance/speed_of_light
        light_travel_time = earth_distance/speed_of_light
        emitted_time = time - light_travel_time

    return _get_barycentric_planet_position(emitted_time, planet_index)

def get_planet(time, planet_index, location=None):
    """
    Get a `~astropy.coordinates.SkyCoord` for a planet as observed from a
    location on Earth.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto). The special value 10
        corresponds to the Earth's Moon.

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, set to
        a Geocentric observer

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the planet
    """
    cartrep = _get_apparent_planet_position(time, planet_index)
    icrs = ICRS(cartrep)
    if location is not None:
        gcrs = icrs.transform_to(GCRS(obstime=time, 
                                 obsgeoloc=u.Quantity(location.geocentric, 
                                 copy=False)))
    else:
        gcrs = icrs.transform_to(GCRS(obstime=time))
    return SkyCoord(gcrs)


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
        a Geocentric observer.

    planet_index : int
        Index of the planet (1-9 for Mercury through Pluto)

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the planet
    """
    return get_planet(time, planet_index=10, location=location)
