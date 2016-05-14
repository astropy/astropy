# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains convenience functions for retrieving solar system
ephemerides from jplephem.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import OrderedDict
import numpy as np
from .sky_coordinate import SkyCoord
from ..utils.data import download_file
from .. import units as u
from ..constants import c as speed_of_light
from .representation import CartesianRepresentation
from .builtin_frames import GCRS, ICRS
from .builtin_frames.utils import get_jd12, cartrepr_from_matmul
from .. import _erfa

__all__ = ["get_body", "get_moon", "get_body_barycentric", "SOLAR_SYSTEM_BODIES"]

KERNEL = None

"""
each value in the BODIES dictionary a list of kernel pairs needed
to find the barycentric position of that object from the JPL kernel.
"""
BODY_NAME_TO_KERNEL_SPEC = OrderedDict(
                                      (('sun', [(0, 10)]),
                                       ('mercury', [(0, 1), (1, 199)]),
                                       ('venus', [(0, 2), (2, 299)]),
                                       ('earth-moon-barycenter', [(0, 3)]),
                                       ('earth',  [(0, 3), (3, 399)]),
                                       ('moon', [(0, 3), (3, 301)]),
                                       ('mars', [(0, 4)]),
                                       ('jupiter', [(0, 5)]),
                                       ('saturn', [(0, 6)]),
                                       ('uranus', [(0, 7)]),
                                       ('neptune', [(0, 8)]),
                                       ('pluto', [(0, 9)]))
                                      )
SOLAR_SYSTEM_BODIES = tuple(BODY_NAME_TO_KERNEL_SPEC.keys())


def _download_spk_file(url=('http://naif.jpl.nasa.gov/pub/naif/'
                            'generic_kernels/spk/planets/de430.bsp'),
                       show_progress=True):
    """
    Get the Satellite Planet Kernel (SPK) file from NASA JPL.

    Download the file from the JPL webpage once and subsequently access a
    cached copy. The default is the file de430.bsp.

    This file is ~120 MB, and covers years ~1550-2650 CE [1]_.

    .. [1] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de430-de431.txt
    """
    return download_file(url, cache=True, show_progress=show_progress)


def _get_kernel(*args, **kwargs):
    """
    Try importing jplephem, download/retrieve from cache the Satellite Planet
    Kernel.
    """
    global KERNEL

    try:
        from jplephem.spk import SPK
    except ImportError:
        raise ImportError("Solar system ephemeris calculations require the "
                          "jplephem package "
                          "(https://pypi.python.org/pypi/jplephem)")

    if KERNEL is None:
        KERNEL = SPK.open(_download_spk_file())
    return KERNEL


def get_body_barycentric(time, body):
    """
    Calculate the barycentric position of the solar system body ``body``
    in cartesian coordinates.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    body : str
        The solar system body to calculate.

        The allowed values for ``body`` can be found in
        ``astropy.coordinates.SOLAR_SYSTEM_BODIES``.

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) position of the body in cartesian coordinates

    Notes
    -----

    """
    # lookup chain
    chain = BODY_NAME_TO_KERNEL_SPEC[body.lower()]

    kernel = _get_kernel()

    jd1, jd2 = get_jd12(time, 'tdb')

    cartesian_position_body = sum([kernel[pair].compute(jd1, jd2) for pair in chain])

    barycen_to_body_vector = u.Quantity(cartesian_position_body, unit=u.km)
    return CartesianRepresentation(barycen_to_body_vector)


def _get_earth_body_vector(time, body, earth_time=None):
    """
    Calculate the vector between the Geocenter and body with ``body``.

    This routine calculates the vector between the Earth's Geocenter and the body
    specified by ``body``.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation.

    body : str
        The solar system body to calculate.

        The allowed values for ``body`` can be found in
        ``astropy.coordinates.SOLAR_SYSTEM_BODIES``.

    earth_time : `~astropy.time.Time`
        Time used for position of Earth. When correcting for light travel time,
        one wants to use different times for the body in question and Earth.
        If this is set to ```None```, the same time is used for both.

    Returns
    -------
    earth_body_vector : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) vector from Geocenter to the body in cartesian coordinates

    earth_distance : `~astropy.units.Quantity`
        Distance between Earth and body.

    Notes
    -----

    """
    earth_time = earth_time if earth_time is not None else time
    earth_loc = get_body_barycentric(earth_time, 'earth')
    body_loc = get_body_barycentric(time, body)

    earth_body_vector = body_loc.xyz - earth_loc.xyz

    earth_distance = np.sqrt(np.sum(earth_body_vector**2, axis=0))
    return earth_body_vector, earth_distance


def _get_apparent_body_position(time, body):
    """
    Calculate the apparent position of body ``body`` in cartesian
    coordinates, given the approximate light travel time to the object.

    Uses ``jplephem`` with the DE430 kernel.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    body : str
        The solar system body to calculate.

        The allowed values for ``body`` can be found in
        ``astropy.coordinates.SOLAR_SYSTEM_BODIES``.

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) apparent position of the body in cartesian coordinates
    """
    # Calculate position given approximate light travel time.
    delta_light_travel_time = 20*u.s
    emitted_time = time
    light_travel_time = 0*u.s
    while np.any(np.fabs(delta_light_travel_time) > 1.0e-8*u.s):
        earth_to_body_vector, earth_distance = _get_earth_body_vector(emitted_time,
                                                                      body, time)
        delta_light_travel_time = light_travel_time - earth_distance/speed_of_light
        light_travel_time = earth_distance/speed_of_light
        emitted_time = time - light_travel_time

    return get_body_barycentric(emitted_time, body)


def get_body(time, body, location=None):
    """
    Get a `~astropy.coordinates.SkyCoord` for a body as observed from a
    location on Earth.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation

    body : str
        The solar system body to calculate.

        The allowed values for ``body`` can be found in
        ``astropy.coordinates.SOLAR_SYSTEM_BODIES``.

    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, set to
        a Geocentric observer

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        Coordinate for the body
    """
    cartrep = _get_apparent_body_position(time, body)
    icrs = ICRS(cartrep)
    if location is not None:
        obsgeoloc, obsgeovel = location.get_gcrs_posvel(time)
        gcrs = icrs.transform_to(GCRS(obstime=time,
                                      obsgeoloc=obsgeoloc,
                                      obsgeovel=obsgeovel))
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
    return get_body(time, body='moon', location=location)


def _apparent_position_in_true_coordinates(skycoord):
    """
    Convert Skycoord in GCRS frame into one in which RA and Dec
    are defined w.r.t to the true equinox and poles of the Earth
    """
    jd1, jd2 = get_jd12(skycoord.obstime, 'tt')
    _, _, _, _, _, _, _, rbpn = _erfa.pn00a(jd1, jd2)
    return SkyCoord(skycoord.frame.realize_frame(cartrepr_from_matmul(rbpn, skycoord)))
