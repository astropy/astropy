# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.coordinates import (EarthLocation, SkyCoord, AltAz, get_sun,
                                 Angle)

import astropy.units as u
import datetime
from astropy.time import Time
import pytz
import numpy as np
################################################################################
# TODO: Temporary solution to IERS tables problems
from astropy.utils.data import download_file
from astropy.utils import iers

iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL,
                                                      cache=True))
################################################################################


from astropy.extern.six import string_types
from .exceptions import TargetNeverUpWarning, TargetAlwaysUpWarning
import warnings


from abc import ABCMeta, abstractmethod

import numpy as np

__all__ = ["Observer", "Target", "FixedTarget", "NonFixedTarget",
           "Constraint", "TimeWindow", "AltitudeRange",
           "AboveAirmass", "Observation"]

#__doctest_requires__ = {'*': ['scipy.integrate']}

def _generate_24hr_grid(t0, start, end, N, for_deriv=False):
    """
    Generate a nearly linearly spaced grid of time durations.

    The midpoints of these grid points will span times from ``t0``+``start``
    to ``t0``+``end``, including the end points, which is useful when taking
    numerical derivatives.

    Parameters
    ----------
    t0 : `~astropy.time.Time`
        Time queried for, grid will be built from or up to this time.

    start : float
        Number of days before/after ``t0`` to start the grid.

    end : float
        Number of days before/after ``t0`` to end the grid.

    N : int
        Number of grid points to generate

    for_deriv : bool
        Generate time series for taking numerical derivative (modify
        bounds)?

    Returns
    -------
    `~astropy.time.Time`
    """

    if for_deriv:
        time_grid = np.concatenate([[start - 1/(N-1)],
                                    np.linspace(start, end, N)[1:-1],
                                    [end + 1/(N-1)]])*u.day
    else:
        time_grid = np.linspace(start, end, N)*u.day

    return t0 + time_grid

class Observer(object):
    """
    A container class for information about an observer's location and
    environment.

    TODO: write this docstring
    """
    @u.quantity_input(elevation=u.m)
    def __init__(self, location=None, timezone='UTC', name=None, latitude=None,
                 longitude=None, elevation=0*u.m, pressure=None,
                 relative_humidity=None, temperature=None, description=None):
        """
        Parameters
        ----------
        name : str
            A short name for the telescope, observatory or location.

        location : `~astropy.coordinates.EarthLocation`
            The location (latitude, longitude, elevation) of the observatory.

        longitude : float, str, `~astropy.units.Quantity` (optional)
            The longitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Longitude` object.

        latitude : float, str, `~astropy.units.Quantity` (optional)
            The latitude of the observing location. Should be valid input for
            initializing a `~astropy.coordinates.Latitude` object.

        elevation : `~astropy.units.Quantity` (optional), default = 0 meters
            The elevation of the observing location, with respect to sea
            level. Defaults to sea level.

        pressure : `~astropy.units.Quantity` (optional)
            The ambient pressure. Defaults to zero (i.e. no atmosphere).

        relative_humidity : float (optional)
            The ambient relative humidity.

        temperature : `~astropy.units.Quantity` (optional)
            The ambient temperature.

        timezone : str or `datetime.tzinfo` (optional)
            The local timezone to assume. If a string, it will be passed through
            `pytz.timezone()` to produce the timezone object.

        description : str (optional)
            A short description of the telescope, observatory or observing
            location.
        """

        self.name = name
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity

        # If lat/long given instead of EarthLocation, convert them
        # to EarthLocation
        if location is None and (latitude is not None and
                                 longitude is not None):
            self.location = EarthLocation.from_geodetic(longitude, latitude,
                                                        elevation)

        elif isinstance(location, EarthLocation):
            self.location = location

        else:
            raise TypeError('Observatory location must be specified with '
                            'either (1) an instance of '
                            'astropy.coordinates.EarthLocation or (2) '
                            'latitude and longitude in degrees as '
                            'accepted by astropy.coordinates.Latitude and '
                            'astropy.coordinates.Latitude.')

        # Accept various timezone inputs, default to UTC
        if isinstance(timezone, datetime.tzinfo):
            self.timezone = timezone
        elif isinstance(timezone, string_types):
            self.timezone = pytz.timezone(timezone)
        else:
            raise TypeError('timezone keyword should be a string, or an '
                            'instance of datetime.tzinfo')

    def astropy_time_to_datetime(self, astropy_time):
        """
        Convert the `~astropy.time.Time` object ``astropy_time`` to a
        localized `~datetime.datetime` object.

        Timezones localized with `~pytz`.

        Parameters
        ----------
        astropy_time : `~astropy.time.Time`
            Scalar or list-like.

        Returns
        -------
        `~datetime.datetime`
            Localized datetime, where the timezone of the datetime is
            set by the ``timezone`` keyword argument of the
            `~astroplan.Observer` constructor.
        """

        if not astropy_time.isscalar:
            return [self.astropy_time_to_datetime(t) for t in astropy_time]

        # Convert astropy.time.Time to a UTC localized datetime (aware)
        utc_datetime = pytz.utc.localize(astropy_time.utc.datetime)

        # Convert UTC to local timezone
        return self.timezone.normalize(utc_datetime)

    def datetime_to_astropy_time(self, date_time):
        """
        Convert the `~datetime.datetime` object ``date_time`` to a
        `~astropy.time.Time` object.

        Timezones localized with `~pytz`. If the ``date_time`` is naive, the
        implied timezone is the ``timezone`` structure of ``self``.

        Parameters
        ----------
        date_time : `~datetime.datetime` or list-like

        Returns
        -------
        `~astropy.time.Time`
            Astropy time object (no timezone information preserved).
        """

        if hasattr(date_time, '__iter__'):
            return Time([self.datetime_to_astropy_time(t) for t in date_time])

        # For timezone-naive datetimes, assign local timezone
        if date_time.tzinfo is None:
            date_time = self.timezone.localize(date_time)

        return Time(date_time, location=self.location)

    def altaz(self, time, target=None, obswl=None):
        """
        Get an `~astropy.coordinates.AltAz` frame or coordinate.

        If ``target`` is None, generates an altitude/azimuth frame. Otherwise,
        calculates the transformation to that frame for the requested ``target``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            The time at which the observation is taking place. Will be used as
            the ``obstime`` attribute in the resulting frame or coordinate. This
            will be passed in as the first argument to the `~astropy.time.Time`
            initializer, so it can be anything that `~astropy.time.Time` will
            accept (including a `~astropy.time.Time` object)

        target : `~astroplan.FixedTarget`, `~astropy.coordinates.SkyCoord`, defaults to `None` (optional)
            Celestial object of interest. If ``target`` is `None`, return the
            `~astropy.coordinates.AltAz` frame without coordinates.

        obswl : `~astropy.units.Quantity` (optional)
            Wavelength of the observation used in the calculation.

        Returns
        -------
        `~astropy.coordinates.AltAz`
            If ``target`` is `None`, returns `~astropy.coordinates.AltAz` frame.
            If ``target`` is not `None`, returns the ``target`` transformed to
            the `~astropy.coordinates.AltAz` frame.

        """
        if not isinstance(time, Time):
            time = Time(time)

        altaz_frame = AltAz(location=self.location, obstime=time,
                            pressure=self.pressure, obswl=obswl,
                            temperature=self.temperature,
                            relative_humidity=self.relative_humidity)

        if target is None:
            return altaz_frame
        else:
            if not (hasattr(target, 'ra') or hasattr(target, 'dec') or
                    hasattr(target, 'coord')):
                raise TypeError('The target must be a coordinate (i.e. a '
                                'FixedTarget or SkyCoord).')

            if hasattr(target, 'coord'):
                coordinate = target.coord
            else:
                coordinate = target
            return coordinate.transform_to(altaz_frame)

    def parallactic_angle(self, time, target):
        '''
        Calculate the parallactic angle.

        Parameters
        ----------
        time : `~astropy.time.Time`
            Observation time.

        target : `~astroplan.FixedTarget` or `~astropy.coordinates.SkyCoord`
            Celestial object of interest.

        Returns
        -------
        `~astropy.coordinates.Angle`
            Parallactic angle.

        Notes
        -----
        The parallactic angle is the angle between the great circle that
        intersects a celestial object and the zenith, and the object's hour
        circle [1]_.

        .. [1] https://en.wikipedia.org/wiki/Parallactic_angle

        '''
        if not isinstance(time, Time):
            time = Time(time)

        if not (hasattr(target, 'ra') or hasattr(target, 'dec') or
                hasattr(target, 'coord')):
            raise TypeError('The target must be a coordinate (i.e. a '
                            'FixedTarget or SkyCoord).')

        if hasattr(target, 'coord'):
            coordinate = target.coord
        else:
            coordinate = target

        # Eqn (14.1) of Meeus' Astronomical Algorithms
        LST = time.sidereal_time('mean', longitude=self.location.longitude)
        H = (LST - coordinate.ra).radian
        q = np.arctan(np.sin(H) /
                      (np.tan(self.location.latitude.radian)*
                       np.cos(coordinate.dec.radian) -
                       np.sin(coordinate.dec.radian)*np.cos(H)))*u.rad

        return Angle(q)

    # Sun-related methods.
    @u.quantity_input(horizon=u.deg)
    def _horiz_cross(self, t, alt, rise_set, horizon=0*u.degree):
        """
        Find time ``t`` when values in array ``a`` go from
        negative to positive or positive to negative (exclude endpoints)

        ``return_limits`` will return nearest times to zero-crossing.

        Parameters
        ----------
        t : `~astropy.time.Time`
            Grid of times
        alt : `~astropy.units.Quantity`
            Grid of altitudes
        rise_set : {"rising",  "setting"}
            Calculate either rising or setting across the horizon
        horizon : float
            Number of degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)
        Returns
        -------
        Returns the lower and upper limits on the time and altitudes
        of the horizon crossing.
        """
        if rise_set == 'rising':
            # Find index where altitude goes from below to above horizon
            condition = (alt[:-1] < horizon) * (alt[1:] > horizon)
        elif rise_set == 'setting':
            # Find index where altitude goes from above to below horizon
            condition = (alt[:-1] > horizon) * (alt[1:] < horizon)

        if not np.any(condition):
            warnmsg = ('Target does not cross horizon={} within 24 '
                       'hours'.format(horizon))
            if (alt > horizon).all():
                warnings.warn(warnmsg, TargetAlwaysUpWarning)
            else:
                warnings.warn(warnmsg, TargetNeverUpWarning)
            return (None, None), (None, None)

        # Isolate horizon crossing
        nearest_index = np.argwhere(condition)[0][0]

        # Capture points on either side of horizon crossing
        lower_t, upper_t = t[nearest_index:nearest_index+2]
        lower_alt, upper_alt = alt[nearest_index:nearest_index+2]

        return (lower_t, upper_t), (lower_alt, upper_alt)

    @u.quantity_input(horizon=u.deg)
    def _two_point_interp(self, times, altitudes, horizon=0*u.deg):
        """
        Do linear interpolation between two ``altitudes`` at
        two ``times`` to determine the time where the altitude
        goes through zero.

        Parameters
        ----------
        times : `~astropy.time.Time`
            Two times for linear interpolation between

        altitudes : array of `~astropy.units.Quantity`
            Two altitudes for linear interpolation between

        horizon : `~astropy.units.Quantity`
            Solve for the time when the altitude is equal to
            reference_alt.

        Returns
        -------
        t : `~astropy.time.Time`
            Time when target crosses the horizon

        """
        if times[0] is None:
            return np.nan
        else:
            slope = (altitudes[1] - altitudes[0])/(times[1].jd - times[0].jd)
            return Time(times[1].jd - ((altitudes[1] - horizon)/slope).value,
                        format='jd')

    def _altitude_trig(self, LST, target):
        """
        Calculate the altitude of ``target`` at local sidereal times ``LST``.

        This method provides a factor of ~3 speed up over calling `altaz`, and
        inherently does *not* take the atmosphere into account.

        Parameters
        ----------
        LST : `~astropy.time.Time`
            Local sidereal times (array)

        target : {`~astropy.coordinates.SkyCoord`, `FixedTarget`} or similar
            Target celestial object's coordinates.

        Returns
        -------
        alt : `~astropy.unit.Quantity`
            Array of altitudes
        """
        alt = np.arcsin(np.sin(self.location.latitude.radian) *
                        np.sin(target.dec) +
                        np.cos(self.location.latitude.radian) *
                        np.cos(target.dec) *
                        np.cos(LST.radian - target.ra.radian))
        return alt

    def _calc_riseset(self, time, target, prev_next, rise_set, horizon, N=150):
        """
        Time at next rise/set of ``target``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : `~astropy.coordinates.SkyCoord`
            Position of target or multiple positions of that target
            at multiple times (if target moves, like the Sun)

        prev_next : str - either 'previous' or 'next'
            Test next rise/set or previous rise/set

        rise_set : str - either 'rising' or 'setting'
            Compute prev/next rise or prev/next set

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        horizon : `~astropy.units.Quantity`
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        N : int
            Number of altitudes to compute when searching for
            rise or set.

        Returns
        -------
        ret1 : `~astropy.time.Time`
            Time of rise/set
        """

        if not isinstance(time, Time):
            time = Time(time)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N)
        else:
            times = _generate_24hr_grid(time, -1, 0, N)

        altitudes = self.altaz(times, target).alt

        horizon_crossing_limits = self._horiz_cross(times, altitudes, rise_set,
                                                    horizon)
        return self._two_point_interp(*horizon_crossing_limits, horizon=horizon)

    def _calc_riseset_brentq(self, time, target, prev_next, rise_set, horizon,
                             N=150, **kwargs):
        """
        Time at next rise/set of ``target`` with `~scipy.optimize.brentq`.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : `~astropy.coordinates.SkyCoord`
            Position of target or multiple positions of that target
            at multiple times (if target moves, like the Sun)

        prev_next : str - either 'previous' or 'next'
            Test next rise/set or previous rise/set

        rise_set : str - either 'rising' or 'setting'
            Compute prev/next rise or prev/next set

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        horizon : `~astropy.units.Quantity`
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        N : int
            Number of altitudes to compute when searching for
            rise or set.

        Returns
        -------
        ret1 : `~astropy.time.Time`
            Time of rise/set
        """
        if not isinstance(time, Time):
            time = Time(time)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N)
        else:
            times = _generate_24hr_grid(time, -1, 0, N)

        altitudes = self.altaz(times, target).alt

        time_limits, alt_limits = self._horiz_cross(times, altitudes, rise_set,
                                                    horizon)

        calc_altitude = lambda t: (self.altaz(Time(t, format='jd'),
                                              target).alt - horizon).value
        from scipy.optimize import brentq
        root = brentq(calc_altitude, time_limits[0].jd, time_limits[1].jd,
                      **kwargs)
        return Time(root, format='jd')

    def _calc_transit(self, time, target, prev_next, antitransit=False, N=150):
        """
        Time at next transit of the meridian of `target`.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : `~astropy.coordinates.SkyCoord`
            Position of target or multiple positions of that target
            at multiple times (if target moves, like the Sun)

        prev_next : str - either 'previous' or 'next'
            Test next rise/set or previous rise/set

        antitransit : bool
            Toggle compute antitransit (below horizon, equivalent to midnight
            for the Sun)

        location : `~astropy.coordinates.EarthLocation`
            Location of observer

        N : int
            Number of altitudes to compute when searching for
            rise or set.

        Returns
        -------
        ret1 : `~astropy.time.Time`
            Time of transit/antitransit
        """
        if not isinstance(time, Time):
            time = Time(time)

        if prev_next == 'next':
            times = _generate_24hr_grid(time, 0, 1, N, for_deriv=True)
        else:
            times = _generate_24hr_grid(time, -1, 0, N, for_deriv=True)

        # The derivative of the altitude with respect to time is increasing
        # from negative to positive values at the anti-transit of the meridian
        if antitransit:
            rise_set = 'rising'
        else:
            rise_set = 'setting'

        altitudes = self.altaz(times, target).alt
        dt = Time((times.jd[1:] + times.jd[:-1])/2, format='jd')
        d_altitudes = altitudes.diff()

        horizon = 0*u.degree # Find when derivative passes through zero
        horizon_crossing_limits = self._horiz_cross(dt, d_altitudes, rise_set,
                                                    horizon)
        return self._two_point_interp(*horizon_crossing_limits, horizon=horizon)

    @u.quantity_input(horizon=u.deg)
    def calc_rise(self, time, target, which='nearest', horizon=0*u.degree):
        """
        Calculate rise time.

        Compute time of the next/previous/nearest rise of the ``target``
        object, where "rise" is defined as the time when the ``target``
        transitions from altitudes below the ``horizon`` to above the
        ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Rise time of target
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_rise = self._calc_riseset(time, target, 'next', 'rising',
                                           horizon)
            if which == 'next':
                return next_rise

        if which == 'previous' or which == 'nearest':
            previous_rise = self._calc_riseset(time, target, 'previous',
                                               'rising', horizon)
            if which == 'previous':
                return previous_rise

        if which == 'nearest':
            if abs(time - previous_rise) < abs(time - next_rise):
                return previous_rise
            else:
                return next_rise
        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    @u.quantity_input(horizon=u.deg)
    def calc_set(self, time, target, which='nearest', horizon=0*u.degree):
        """
        Calculate set time.

        Compute time of the next/previous/nearest set of ``target``, where
        "set" is defined as when the ``target`` transitions from altitudes
        above ``horizon`` to below ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunset relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Set time of target.
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_set = self._calc_riseset(time, target, 'next', 'setting',
                                          horizon)
            if which == 'next':
                return next_set

        if which == 'previous' or which == 'nearest':
            previous_set = self._calc_riseset(time, target, 'previous',
                                              'setting', horizon)
            if which == 'previous':
                return previous_set

        if which == 'nearest':
            if abs(time - previous_set) < abs(time - next_set):
                return previous_set
            else:
                return next_set

        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    def calc_meridian_transit(self, time, target, which='nearest'):
        """
        Calculate time at the transit of the meridian.

        Compute time of the next/previous/nearest transit of the ``target``
        object.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Transit time of target
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_transit = self._calc_transit(time, target, 'next')
            if which == 'next':
                return next_transit

        if which == 'previous' or which == 'nearest':
            previous_transit = self._calc_transit(time, target, 'previous')
            if which == 'previous':
                return previous_transit

        if which == 'nearest':
            if abs(time - previous_transit) < abs(time - next_transit):
                return previous_transit
            else:
                return next_transit
        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    def calc_meridian_antitransit(self, time, target, which='nearest'):
        """
        Calculate time at the antitransit of the meridian.

        Compute time of the next/previous/nearest antitransit of the ``target``
        object.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Antitransit time of target
        """
        if not isinstance(time, Time):
            time = Time(time)

        if which == 'next' or which == 'nearest':
            next_antitransit = self._calc_transit(time, target, 'next',
                                                  antitransit=True)
            if which == 'next':
                return next_antitransit

        if which == 'previous' or which == 'nearest':
            previous_antitransit = self._calc_transit(time, target, 'previous',
                                                      antitransit=True)
            if which == 'previous':
                return previous_antitransit

        if which == 'nearest':
            if abs(time - previous_antitransit) < abs(time - next_antitransit):
                return previous_antitransit
            else:
                return next_antitransit
        raise ValueError('"which" kwarg must be "next", "previous" or '
                         '"nearest".')

    @u.quantity_input(horizon=u.deg)
    def sunrise(self, time, which='nearest', horizon=0*u.degree):
        """
        Time of sunrise.

        Compute time of the next/previous/nearest sunrise, where
        sunrise is defined as when the Sun transitions from altitudes
        below ``horizon`` to above ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which sunrise relative to the present ``time`` would you
            like to calculate.

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Time of sunrise
        """
        return self.calc_rise(time, get_sun(time), which, horizon)

    @u.quantity_input(horizon=u.deg)
    def sunset(self, time, which='nearest', horizon=0*u.degree):
        """
        Time of sunset.

        Compute time of the next/previous/nearest sunset, where
        sunset is defined as when the Sun transitions from altitudes
        below ``horizon`` to above ``horizon``.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which sunset relative to the present ``time`` would you
            like to calculate

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        Returns
        -------
        `~astropy.time.Time`
            Time of sunset
        """
        return self.calc_set(time, get_sun(time), which, horizon)

    def noon(self, time, which='nearest'):
        """
        Time at solar noon.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which noon relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time at solar noon
        """
        return self.calc_meridian_transit(time, get_sun(time), which)

    def midnight(self, time, which='nearest'):
        """
        Time at solar midnight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which noon relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time at solar midnight
        """
        return self.calc_meridian_antitransit(time, get_sun(time), which)

    # Twilight convenience functions

    def evening_astronomical(self, time, which='nearest'):
        """
        Time at evening astronomical (-18 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sunset(time, which, horizon=-18*u.degree)

    def evening_nautical(self, time, which='nearest'):

        """
        Time at evening nautical (-12 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sunset(time, which, horizon=-12*u.degree)

    def evening_civil(self, time, which='nearest'):
        """
        Time at evening civil (-6 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sunset(time, which, horizon=-6*u.degree)

    def morning_astronomical(self, time, which='nearest'):
        """
        Time at morning astronomical (-18 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sunrise(time, which, horizon=-18*u.degree)

    def morning_nautical(self, time, which='nearest'):
        """
        Time at morning nautical (-12 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of twilight
        """
        return self.sunrise(time, which, horizon=-12*u.degree)

    def morning_civil(self, time, which='nearest'):
        """
        Time at morning civil (-6 degree) twilight.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observations. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        which : {'next', 'previous', 'nearest'}
            Choose which twilight relative to the present ``time`` would you
            like to calculate. Default is nearest.

        Returns
        -------
        `~astropy.time.Time`
            Time of sunset
        """
        return self.sunrise(time, which, horizon=-6*u.degree)

    # Moon-related methods.

    def moonrise(self, time, **kwargs):
        """
        Returns the local moonrise time.

        The default moonrise returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moonset(self, time, **kwargs):
        """
        Returns the local moonset time.

        The default moonset returned is the next one to occur.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).

        Keywords: str, optional
            previous
            next
        """
        raise NotImplementedError()

    def moon_illumination(self, time):
        """
        Returns a float giving the percent illumation.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).
        """
        raise NotImplementedError()

    def moon_position(self, time):
        """
        Returns the position of the moon in alt/az.

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object).
        """
        raise NotImplementedError()

    @u.quantity_input(horizon=u.deg)
    def can_see(self, time, target, horizon=0*u.degree, return_altaz=False):
        """
        Is ``target`` above ``horizon`` at this ``time``?

        Parameters
        ----------
        time : `~astropy.time.Time` or other (see below)
            Time of observation. This will be passed in as the first argument to
            the `~astropy.time.Time` initializer, so it can be anything that
            `~astropy.time.Time` will accept (including a `~astropy.time.Time`
            object)

        target : coordinate object (i.e. `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`)
            Target celestial object

        horizon : `~astropy.units.Quantity` (optional), default = zero degrees
            Degrees above/below actual horizon to use
            for calculating rise/set times (i.e.,
            -6 deg horizon = civil twilight, etc.)

        return_altaz : bool (optional)
            Also return the '~astropy.coordinates.AltAz' coordinate.
        """
        if not isinstance(time, Time):
            time = Time(time)

        altaz = self.altaz(time, target)
        observable = altaz.alt > horizon

        if not return_altaz:
            return observable
        else:
            return observable, altaz

class Target(object):
    """
    This is an abstract base class -- you can't instantiate
    examples of this class, but must work with one of its
    subclasses such as ``FixedTarget`` or ``NonFixedTarget``.

    Would need to import six, abc to make this a metaclass?
    """
    __metaclass__ = ABCMeta

    def __init__(self, name=None, ra=None, dec=None, marker=None):
        """
        Defines a single observation target.

        Parameters
        ----------
        name : str, optional

        ra : WHAT TYPE IS ra ?

        dec : WHAT TYPE IS dec ?

        marker : str, optional
            User-defined markers to differentiate between different types
            of targets (e.g., guides, high-priority, etc.).
        """
        raise NotImplementedError()

    @property
    def ra(self):
        """
        Right ascension.
        """
        if isinstance(self, FixedTarget):
            return self.coord.ra
        raise NotImplementedError()

    @property
    def dec(self):
        """
        Declination.
        """
        if isinstance(self, FixedTarget):
            return self.coord.dec
        raise NotImplementedError()


class FixedTarget(Target):
    """
    An object that is "fixed" with respect to the celestial sphere.
    """
    def __init__(self, coord, name=None, **kwargs):
        """
        TODO: Docstring.
        """
        if not (hasattr(coord, 'transform_to') and
                hasattr(coord, 'represent_as')):
            raise TypeError('`coord` must be a coordinate object.')

        self.name = name
        self.coord = coord

    @classmethod
    def from_name(cls, query_name, name=None, **kwargs):
        """
        Initialize a `FixedTarget` by querying for a name, using the machinery
        in `~astropy.coordinates.SkyCoord.from_name`.
        """
        # Allow manual override for name keyword so that the target name can
        # be different from the query name, otherwise assume name=queryname.
        if name is None:
            name = query_name
        return cls(SkyCoord.from_name(query_name), name=name, **kwargs)

class NonFixedTarget(Target):
    """
    Placeholder for future function.
    """


class Constraint(object):
    """
    An object containing observational constraints.

    A Constraints object is used in conjunction with a Target
    and an Observer object (via the apply_constraints method) to find out
    if a particular target is visible to the observer.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def apply_constraints(self, target, observer, constraint_list):
        """
        Returns information on a target's visibility.

        Finds out if a Target is observable by an Observer given a list
        of Constraint objects.  The list must contain at least one
        Constraint object.

        Parameters
        ----------
        target : WHAT TYPE IS Target OBJECT ?

        constraint_list : WHAT TYPE IS constraint_list ? `numpy.array` ??
        """
        raise NotImplementedError


class TimeWindow(Constraint):
    """
    An object containing start and end times for an observation.
    """

    def __init__(self, start, end):
        """
        Initializes a TimeWindow object.

        Parameters
        ----------
        start : STRING OR astropy.time OBJECT ?

        end : STRING OR astropy.time OBJECT ?
        """
        raise NotImplementedError


class AltitudeRange(Constraint):
    """
    An object containing upper and lower altitude limits.
    """

    def __init__(self, low, high):
        """
        Initializes an AltitudeRange object.

        Parameters
        ----------
        low : `~astropy.units.Quantity`

        high : `~astropy.units.Quantity`
        """
        raise NotImplementedError


class AboveAirmass(Constraint):
    """
    An object containing an airmass lower limit.
    """

    def __init__(self, low):
        """
        Initializes an AboveAirmass object.

        Parameters
        ----------
        low : float
        """
        raise NotImplementedError


class Observation(object):
    """
    Comments.
    """

    def __init__(self, target, time):
        """
        Initializes an Observation object.

        Parameters
        ----------
        target : WHAT TYPE IS Target OBJECT ?

        date : WHAT TYPE IS date OBJECT ?
        """
        raise NotImplementedError()

    # Observability properties.

    @property
    def alt(self):
        """
        Altitude at time of observation.
        """
        raise NotImplementedError()

    @property
    def az(self):
        """
        Azimuth at time of observation.
        """
        raise NotImplementedError()

    @property
    def airmass(self):
        """
        Airmass.
        """
        raise NotImplementedError()

    @property
    def pang(self):
        """
        Parallactic angle.
        """
        raise NotImplementedError()

    @property
    def ha(self):
        """
        Hour angle.
        """
        raise NotImplementedError()

    @property
    def moon_sep(self):
        """
        Separation between moon and object at time of observation.
        """
        raise NotImplementedError()

    # Time properties.

    @property
    def ut(self):
        """
        Time of observation in UST.
        """
        raise NotImplementedError()

    @property
    def lt(self):
        """
        Time of observation in local time.
        """
        raise NotImplementedError()

    @property
    def gmst(self):
        """
        Time of observation in GMST.
        """
        raise NotImplementedError()

    @property
    def lmst(self):
        """
        Time of observation in local mean sidereal time.
        """
        raise NotImplementedError()
