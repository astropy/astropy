# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains convenience functions for retrieving solar system
ephemerides from jplephem.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import six
from .sky_coordinate import SkyCoord
from ..utils.data import download_file
from .. import units as u
from ..constants import c as speed_of_light
from .representation import CartesianRepresentation
from .builtin_frames import GCRS, ICRS
from .builtin_frames.utils import get_jd12

__all__ = ["get_body", "get_moon"]

KERNEL = None

""" 
each value in the BODIES dictionary a list of kernel pairs needed
to find the barycentric position of that object from the JPL kernel.
"""
BODIES = {'sun': [(0,10)],
          'mercury': [(0, 1), (1, 199)], 
          'venus': [(0,2), (2, 299)], 
          'earth-moon-barycenter': [(0,3)],
          'earth':  [(0, 3), (3, 399)],
          'moon': [(0, 3), (3, 301)],
          'mars': [(0, 4)],
          'jupiter': [(0,5)],
          'saturn': [(0,6)],
          'uranus': [(0,7)],
          'neptune': [(0,8)],
          'pluto': [(0,9)]}          

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


def _get_barycentric_body_position(time, body_key):
    """
    Calculate the barycentric position of solar system body ``body_key`` in cartesian coordinates.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    body_key : int, str
        The solar system body to calculate. Can be a string, e.g. 'moon',
        or a valid integer index.

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) position of the body in cartesian coordinates

    Notes
    -----

    """
    # if we have an int as a body index, then create chain to that body
    if isinstance(body_key, int):
        # chain consists of single step, from 0 (solar system barycenter) to body_key
        chain = [(0, body_key)]
    elif isinstance(body_key, six.string_types):
        try:
            chain = BODIES[body_key]
        except KeyError:
            raise ValueError("Unknown body index '" + body_key + """', valid entries are
                              contained within solar_system.BODIES""")
    else:
        raise ValueError("'body_key' must be a string or integer") 
        
    kernel = _get_kernel()
    
    # are all the kernel pairs in this kernel?
    valid_chain = np.all([c in kernel.pairs.keys() for c in chain])
    if not valid_chain:
        raise ValueError("Postion of this body cannot be calculated using jpl kernel")
        
    jd1, jd2 = get_jd12(time, 'tdb')

    cartesian_position_body = np.sum([kernel[pair].compute(jd1, jd2) for pair in chain],
                                       axis=0)

    barycen_to_body_vector = u.Quantity(cartesian_position_body, unit=u.km)
    return CartesianRepresentation(barycen_to_body_vector)


def _get_earth_body_vector(time, body_key):
    """
    Calculate the vector between the Geocenter and body with ``body_key``.

    This routine calculates the vector between the Earth's Geocenter and the body
    specified by ``body_key``.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    body_key : int, str
        The solar system body to calculate. Can be a string, e.g. 'moon',
        or a valid integer index.

    Returns
    -------
    earth_body_vector : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) vector from Geocenter to the body in cartesian coordinates

    earth_distance : `~astropy.units.Quantity`
        Distance between Earth and body.

    Notes
    -----

    """
    earth_loc = _get_barycentric_body_position(time, 'earth')
    body_loc = _get_barycentric_body_position(time, body_key)
    earth_body_vector = body_loc.xyz - earth_loc.xyz
    earth_distance = np.sqrt(np.sum(earth_body_vector**2,axis=0))
    return earth_body_vector, earth_distance

def _get_apparent_body_position(time, body_key):
    """
    Calculate the apparent position of body ``body_key`` in cartesian
    coordinates, given the approximate light travel time to the object.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    body_key : int, str
        The solar system body to calculate. Can be a string, e.g. 'moon',
        or a valid integer index.

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) apparent position of the body in cartesian coordinates
    """
    # Get distance of body at `time`
    earth_to_body_vector, earth_distance = _get_earth_body_vector(time, body_key)

    # The apparent position depends on the time that the light was emitted from
    # the distant body, so subtract off the light travel time
    light_travel_time = earth_distance/speed_of_light
    emitted_time = time - light_travel_time

    # Calculate position given approximate light travel time.
    delta_light_travel_time = 20*u.s
    while (np.fabs(delta_light_travel_time)) > 1.0e-8*u.s:
        earth_to_body_vector, earth_distance = _get_earth_body_vector(emitted_time, body_key)
        delta_light_travel_time = light_travel_time - earth_distance/speed_of_light
        light_travel_time = earth_distance/speed_of_light
        emitted_time = time - light_travel_time

    return _get_barycentric_body_position(emitted_time, body_key)


def get_body(time, body_key, location=None):
    """
    Get a `~astropy.coordinates.SkyCoord` for a body as observed from a
    location on Earth.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    body_key : int
        Index of the body (1-9 for Mercury through Pluto). The special value 10
        corresponds to the Earth's Moon.

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, set to
        a Geocentric observer

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the body
    """
    cartrep = _get_apparent_body_position(time, body_key)
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
    Get a `~astropy.coordinates.SkyCoord` for the Earth's Moon as observed
    from a location on Earth.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, set to
        a Geocentric observer.

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the Moon
    """
    return get_body(time, body_key='moon', location=location)
