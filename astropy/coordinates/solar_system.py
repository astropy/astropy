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
from ..utils.decorators import classproperty
from ..utils.state import ScienceState
from ..utils import indent
from .. import units as u
from .. import _erfa as erfa
from ..constants import c as speed_of_light
from .representation import CartesianRepresentation
from .orbital_elements import calc_moon
from .builtin_frames import GCRS, ICRS
from .builtin_frames.utils import get_jd12
from ..extern import six

__all__ = ["get_body", "get_moon", "get_body_barycentric",
           "get_body_barycentric_posvel", "solar_system_ephemeris"]


DEFAULT_JPL_EPHEMERIS = 'de430'

"""List of kernel pairs needed to calculate positions of a given object."""
BODY_NAME_TO_KERNEL_SPEC = OrderedDict(
                                      (('sun', [(0, 10)]),
                                       ('mercury', [(0, 1), (1, 199)]),
                                       ('venus', [(0, 2), (2, 299)]),
                                       ('earth-moon-barycenter', [(0, 3)]),
                                       ('earth', [(0, 3), (3, 399)]),
                                       ('moon', [(0, 3), (3, 301)]),
                                       ('mars', [(0, 4)]),
                                       ('jupiter', [(0, 5)]),
                                       ('saturn', [(0, 6)]),
                                       ('uranus', [(0, 7)]),
                                       ('neptune', [(0, 8)]),
                                       ('pluto', [(0, 9)]))
                                      )

"""Indices to the plan94 routine for the given object."""
PLAN94_BODY_NAME_TO_PLANET_INDEX = OrderedDict(
    (('mercury', 1),
     ('venus', 2),
     ('earth-moon-barycenter', 3),
     ('mars', 4),
     ('jupiter', 5),
     ('saturn', 6),
     ('uranus', 7),
     ('neptune', 8)))

_EPHEMERIS_NOTE = """
You can either give an explicit ephemeris or use a default, which is normally
a built-in ephemeris that does not require ephemeris files.  To change
the default to be the JPL ephemeris::

    >>> from astropy.coordinates import solar_system_ephemeris
    >>> solar_system_ephemeris.set('jpl')  # doctest: +SKIP

Use of any JPL ephemeris requires the jplephem package
(https://pypi.python.org/pypi/jplephem).
If needed, the ephemeris file will be downloaded (and cached).

One can check which bodies are covered by a given ephemeris using::
    >>> solar_system_ephemeris.bodies
    ('earth', 'sun', 'moon', 'mercury', 'venus', 'earth-moon-barycenter', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune')
"""[1:-1]


class solar_system_ephemeris(ScienceState):
    """Default ephemerides for calculating positions of Solar-System bodies.

    This can be one of the following::

    - 'builtin': polynomial approximations to the orbital elements.
    - 'de430' or 'de432s': short-cuts for recent JPL dynamical models.
    - 'jpl': Alias for the default JPL ephemeris (currently, 'de430').
    - URL: (str) The url to a SPK ephemeris in SPICE binary (.bsp) format.
    - `None`: Ensure an Exception is raised without an explicit ephemeris.

    The default is 'builtin', which uses the ``epv00`` and ``plan94``
    routines from the ``erfa`` implementation of the Standards Of Fundamental
    Astronomy library.

    Notes
    -----
    Any file required will be downloaded (and cached) when the state is set.
    The default Satellite Planet Kernel (SPK) file from NASA JPL (de430) is
    ~120MB, and covers years ~1550-2650 CE [1]_.  The smaller de432s file is
    ~10MB, and covers years 1950-2050 [2]_.  Older versions of the JPL
    ephemerides (such as the widely used de200) can be used via their URL [3]_.

    .. [1] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de430-de431.txt
    .. [2] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de432s.txt
    .. [3] http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/
    """
    _value = 'builtin'
    _kernel = None

    @classmethod
    def validate(cls, value):
        # make no changes if value is None
        if value is None:
            return cls._value
        # Set up Kernel; if the file is not in cache, this will download it.
        cls.get_kernel(value)
        return value

    @classmethod
    def get_kernel(cls, value):
        # ScienceState only ensures the `_value` attribute is up to date,
        # so we need to be sure any kernel returned is consistent.
        if cls._kernel is None or cls._kernel.origin != value:
            if cls._kernel is not None:
                cls._kernel.daf.file.close()
                cls._kernel = None
            kernel = _get_kernel(value)
            if kernel is not None:
                kernel.origin = value
            cls._kernel = kernel
        return cls._kernel

    @classproperty
    def kernel(cls):
        return cls.get_kernel(cls._value)

    @classproperty
    def bodies(cls):
        if cls._value is None:
            return None
        if cls._value.lower() == 'builtin':
            return (('earth', 'sun', 'moon') +
                    tuple(PLAN94_BODY_NAME_TO_PLANET_INDEX.keys()))
        else:
            return tuple(BODY_NAME_TO_KERNEL_SPEC.keys())


def _get_kernel(value):
    """
    Try importing jplephem, download/retrieve from cache the Satellite Planet
    Kernel corresponding to the given ephemeris.
    """
    if value is None or value.lower() == 'builtin':
        return None

    if value.lower() == 'jpl':
        value = DEFAULT_JPL_EPHEMERIS

    if value.lower() in ('de430', 'de432s'):
        value = ('http://naif.jpl.nasa.gov/pub/naif/generic_kernels'
                 '/spk/planets/{:s}.bsp'.format(value.lower()))
    else:
        try:
            six.moves.urllib.parse.urlparse(value)
        except Exception:
            raise ValueError('{} was not one of the standard strings and '
                             'could not be parsed as a URL'.format(value))

    try:
        from jplephem.spk import SPK
    except ImportError:
        raise ImportError("Solar system JPL ephemeris calculations require "
                          "the jplephem package "
                          "(https://pypi.python.org/pypi/jplephem)")

    return SPK.open(download_file(value, cache=True))


def _get_body_barycentric_posvel(body, time, ephemeris=None,
                                 get_velocity=True):
    """Calculate the barycentric position (and velocity) of a solar system body.

    Parameters
    ----------
    body : str or other
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    ephemeris : str, optional
        Ephemeris to use.  By default, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set``
    get_velocity : bool, optional
        Whether or not to calculate the velocity as well as the position.

    Returns
    -------
    position : `~astropy.coordinates.CartesianRepresentation` or tuple
        Barycentric (ICRS) position or tuple of position and velocity.

    Notes
    -----
    No velocity can be calculated with the built-in ephemeris for the Moon.

    Whether or not velocities are calculated makes little difference for the
    built-in ephemerides, but for most JPL ephemeris files, the execution time
    roughly doubles.
    """

    if ephemeris is None:
        ephemeris = solar_system_ephemeris.get()
        if ephemeris is None:
            raise ValueError(_EPHEMERIS_NOTE)
        kernel = solar_system_ephemeris.kernel
    else:
        kernel = _get_kernel(ephemeris)

    jd1, jd2 = get_jd12(time, 'tdb')
    if kernel is None:
        body = body.lower()
        earth_pv_helio, earth_pv_bary = erfa.epv00(jd1, jd2)
        if body == 'earth':
            body_pv_bary = earth_pv_bary

        elif body == 'moon':
            if get_velocity:
                raise KeyError("the Moon's velocity cannot be calculated with "
                               "the '{0}' ephemeris.".format(ephemeris))
            return calc_moon(time).cartesian

        else:
            sun_pv_bary = earth_pv_bary - earth_pv_helio
            if body == 'sun':
                body_pv_bary = sun_pv_bary
            else:
                try:
                    body_index = PLAN94_BODY_NAME_TO_PLANET_INDEX[body]
                except KeyError:
                    raise KeyError("{0}'s position and velocity cannot be "
                                   "calculated with the '{1}' ephemeris."
                                   .format(body, ephemeris))
                body_pv_helio = erfa.plan94(jd1, jd2, body_index)
                body_pv_bary = body_pv_helio + sun_pv_bary

        body_pos_bary = CartesianRepresentation(
            body_pv_bary[..., 0, :], unit=u.au, xyz_axis=-1, copy=False)
        if get_velocity:
            body_vel_bary = CartesianRepresentation(
                body_pv_bary[..., 1, :], unit=u.au/u.day, xyz_axis=-1,
                copy=False)

    else:
        if isinstance(body, six.string_types):
            # Look up kernel chain for JPL ephemeris, based on name
            try:
                kernel_spec = BODY_NAME_TO_KERNEL_SPEC[body.lower()]
            except KeyError:
                raise KeyError("{0}'s position cannot be calculated with "
                               "the {1} ephemeris.".format(body, ephemeris))
        else:
            # otherwise, assume the user knows what their doing and intentionally
            # passed in a kernel chain
            kernel_spec = body

        # jplephem cannot handle multi-D arrays, so convert to 1D here.
        jd1_shape = getattr(jd1, 'shape', ())
        if len(jd1_shape) > 1:
            jd1, jd2 = jd1.ravel(), jd2.ravel()
        # Note that we use the new jd1.shape here to create a 1D result array.
        # It is reshaped below.
        body_posvel_bary = np.zeros((2 if get_velocity else 1, 3) +
                                     getattr(jd1, 'shape', ()))
        for pair in kernel_spec:
            spk = kernel[pair]
            if spk.data_type == 3:
                # Type 3 kernels contain both position and velocity.
                posvel = spk.compute(jd1, jd2)
                if get_velocity:
                    body_posvel_bary += posvel.reshape(body_posvel_bary.shape)
                else:
                    body_posvel_bary[0] += posvel[:4]
            else:
                # spk.generate first yields the position and then the
                # derivative. If no velocities are desired, body_posvel_bary
                # has only one element and thus the loop ends after a single
                # iteration, avoiding the velocity calculation.
                for body_p_or_v, p_or_v in zip(body_posvel_bary,
                                               spk.generate(jd1, jd2)):
                    body_p_or_v += p_or_v

        body_posvel_bary.shape = body_posvel_bary.shape[:2] + jd1_shape
        body_pos_bary = CartesianRepresentation(body_posvel_bary[0],
                                                unit=u.km, copy=False)
        if get_velocity:
            body_vel_bary = CartesianRepresentation(body_posvel_bary[1],
                                                    unit=u.km/u.day, copy=False)

    return (body_pos_bary, body_vel_bary) if get_velocity else body_pos_bary


def get_body_barycentric_posvel(body, time, ephemeris=None):
    """Calculate the barycentric position and velocity of a solar system body.

    Parameters
    ----------
    body : str or other
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    ephemeris : str, optional
        Ephemeris to use.  By default, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set``

    Returns
    -------
    position, velocity : tuple of `~astropy.coordinates.CartesianRepresentation`
        Tuple of barycentric (ICRS) position and velocity.

    See also
    --------
    get_body_barycentric : to calculate position only.
        This is faster by about a factor two for JPL kernels, but has no
        speed advantage for the built-in ephemeris.

    Notes
    -----
    The velocity cannot be calculated for the Moon.  To just get the position,
    use :func:`~astropy.coordinates.get_body_barycentric`.

    """
    return _get_body_barycentric_posvel(body, time, ephemeris)


get_body_barycentric_posvel.__doc__ += indent(_EPHEMERIS_NOTE)[4:]


def get_body_barycentric(body, time, ephemeris=None):
    """Calculate the barycentric position of a solar system body.

    Parameters
    ----------
    body : str or other
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    ephemeris : str, optional
        Ephemeris to use.  By default, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set``

    Returns
    -------
    position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) position of the body in cartesian coordinates

    See also
    --------
    get_body_barycentric_posvel : to calculate both position and velocity.

    Notes
    -----
    """
    return _get_body_barycentric_posvel(body, time, ephemeris,
                                        get_velocity=False)


get_body_barycentric.__doc__ += indent(_EPHEMERIS_NOTE)[4:]


def _get_apparent_body_position(body, time, ephemeris):
    """Calculate the apparent position of body ``body`` relative to Earth.

    This corrects for the light-travel time to the object.

    Parameters
    ----------
    body : str or other
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    ephemeris : str, optional
        Ephemeris to use.  By default, use the one set with
        ``~astropy.coordinates.solar_system_ephemeris.set``

    Returns
    -------
    cartesian_position : `~astropy.coordinates.CartesianRepresentation`
        Barycentric (ICRS) apparent position of the body in cartesian coordinates
    """
    if ephemeris is None:
        ephemeris = solar_system_ephemeris.get()
    # builtin ephemeris and moon is a special case, with no need to account for
    # light travel time, since this is already included in the Meeus algorithm
    # used.
    if ephemeris == 'builtin' and body.lower() == 'moon':
        return get_body_barycentric(body, time, ephemeris)

    # Calculate position given approximate light travel time.
    delta_light_travel_time = 20. * u.s
    emitted_time = time
    light_travel_time = 0. * u.s
    earth_loc = get_body_barycentric('earth', time, ephemeris)
    while np.any(np.fabs(delta_light_travel_time) > 1.0e-8*u.s):
        body_loc = get_body_barycentric(body, emitted_time, ephemeris)
        earth_distance = (body_loc - earth_loc).norm()
        delta_light_travel_time = (light_travel_time -
                                   earth_distance/speed_of_light)
        light_travel_time = earth_distance/speed_of_light
        emitted_time = time - light_travel_time

    return get_body_barycentric(body, emitted_time, ephemeris)


_get_apparent_body_position.__doc__ += indent(_EPHEMERIS_NOTE)[4:]


def get_body(body, time, location=None, ephemeris=None):
    """
    Get a `~astropy.coordinates.SkyCoord` for a solar system body as observed
    from a location on Earth in the `~astropy.coordinates.GCRS` reference
    system.

    Parameters
    ----------
    body : str or other
        The solar system body for which to calculate positions.  Can also be a
        kernel specifier (list of 2-tuples) if the ``ephemeris`` is a JPL
        kernel.
    time : `~astropy.time.Time`
        Time of observation.
    location : `~astropy.coordinates.EarthLocation`, optional
        Location of observer on the Earth.  If not given, will be taken from
        ``time`` (if not present, a geocentric observer will be assumed).
    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set`` (which is
        set to 'builtin' by default).

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        GCRS Coordinate for the body

    Notes
    -----
    """
    if location is None:
        location = time.location

    cartrep = _get_apparent_body_position(body, time, ephemeris)
    icrs = ICRS(cartrep)
    if location is not None:
        obsgeoloc, obsgeovel = location.get_gcrs_posvel(time)
        gcrs = icrs.transform_to(GCRS(obstime=time,
                                      obsgeoloc=obsgeoloc,
                                      obsgeovel=obsgeovel))
    else:
        gcrs = icrs.transform_to(GCRS(obstime=time))
    return SkyCoord(gcrs)


get_body.__doc__ += indent(_EPHEMERIS_NOTE)[4:]


def get_moon(time, location=None, ephemeris=None):
    """
    Get a `~astropy.coordinates.SkyCoord` for the Earth's Moon as observed
    from a location on Earth in the `~astropy.coordinates.GCRS` reference
    system.

    Parameters
    ----------
    time : `~astropy.time.Time`
        Time of observation
    location : `~astropy.coordinates.EarthLocation`
        Location of observer on the Earth. If none is supplied, taken from
        ``time`` (if not present, a geocentric observer will be assumed).
    ephemeris : str, optional
        Ephemeris to use.  If not given, use the one set with
        ``astropy.coordinates.solar_system_ephemeris.set`` (which is
        set to 'builtin' by default).

    Returns
    -------
    skycoord : `~astropy.coordinates.SkyCoord`
        GCRS Coordinate for the Moon

    Notes
    -----
    """

    return get_body('moon', time, location=location, ephemeris=ephemeris)


get_moon.__doc__ += indent(_EPHEMERIS_NOTE)[4:]


def _apparent_position_in_true_coordinates(skycoord):
    """
    Convert Skycoord in GCRS frame into one in which RA and Dec
    are defined w.r.t to the true equinox and poles of the Earth
    """
    jd1, jd2 = get_jd12(skycoord.obstime, 'tt')
    _, _, _, _, _, _, _, rbpn = erfa.pn00a(jd1, jd2)
    return SkyCoord(skycoord.frame.realize_frame(
        skycoord.cartesian.transform(rbpn)))
