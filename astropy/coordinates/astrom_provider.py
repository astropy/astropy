# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains a helper function to fill erfa.astrom struct and a
ScienceState, which allows to speed up coordinate transformations at the
expense of accuracy.
"""
from abc import ABCMeta, abstractmethod
import warnings

import numpy as np
import erfa

from ..time import Time
from ..utils.state import ScienceState
# from ..utils import indent
from .. import units as u
from .builtin_frames.utils import (
    get_jd12, get_cip, prepare_earth_position_vel, get_polar_motion, get_dut1utc,
    pav2pv
)

from ..utils.exceptions import AstropyWarning

__all__ = []


class BaseAstromProvider(metaclass=ABCMeta):
    '''
    Baseclass for classes providing utility access to erfa astrometry functions.
    '''

    @abstractmethod
    def apci(self, frame):
        '''
        Call ``erfa.apci`` with the correct information extracted from ``frame``.
        '''
        pass

    @abstractmethod
    def apcs(self, frame):
        '''
        Call ``erfa.apcs`` with the correct information extracted from ``frame``.
        '''
        pass

    @abstractmethod
    def apio13(self, frame):
        '''
        Call ``erfa.apcs`` with the correct information extracted from ``frame``.
        '''
        pass


class AstromProvider(BaseAstromProvider):
    '''
    The default provider for astrometry values.
    A utility class to extract the necessary arguments for
    erfa functions from frame attributes.
    '''

    @staticmethod
    def apci(frame):
        jd1_tt, jd2_tt = get_jd12(frame.obstime, 'tt')
        cip = get_cip(jd1_tt, jd2_tt)
        earth_pv, earth_heliocentric = prepare_earth_position_vel(frame.obstime)
        return erfa.apci(jd1_tt, jd2_tt, earth_pv, earth_heliocentric, *cip)

    @staticmethod
    def apcs(frame):
        jd1_tt, jd2_tt = get_jd12(frame.obstime, 'tt')
        pv = pav2pv(
            frame.obsgeoloc.get_xyz(xyz_axis=-1).value,
            frame.obsgeovel.get_xyz(xyz_axis=-1).value
        )
        earth_pv, earth_heliocentric = prepare_earth_position_vel(frame.obstime)
        return erfa.apcs(jd1_tt, jd2_tt, pv, earth_pv, earth_heliocentric)

    @staticmethod
    def apio13(frame):
        lon, lat, height = frame.location.to_geodetic('WGS84')

        jd1_utc, jd2_utc = get_jd12(frame.obstime, 'utc')
        dut1utc = get_dut1utc(frame.obstime)

        return erfa.apio13(
            jd1_utc, jd2_utc, dut1utc,
            lon.to_value(u.radian),
            lat.to_value(u.radian),
            height.to_value(u.m),
            *get_polar_motion(frame.obstime),
            # all below are already in correct units because they are QuantityFrameAttribues
            frame.pressure.value,
            frame.temperature.value,
            frame.relative_humidity.value,
            frame.obswl.value,
        )


class InterpolatingAstromProvider(AstromProvider):
    '''
    A provider for astrometry values that does not call erfa
    for each individual timestamp but interpolates the linearly
    between support points.

    This can dramatically speed up coordinate transformations when
    obstime is an array of many values (factors of 10 to > 100 depending
    on the selected resolution and the time range of the values).
    '''

    @u.quantity_input(time_resolution=u.day)
    def __init__(self, time_resolution):
        if time_resolution.to_value(u.ms) < 1:
            warnings.warn(
                f'Using {self.__class__.__name__} with `time_resolution`'
                ' below 1ms might not improve performance',
                AstropyWarning
            )
        self.mjd_resolution = time_resolution.to_value(u.day)

    def _get_support_points(self, obstime):
        '''
        Calculate support points for the interpolation.
        This includes the minimum and maximum timestamp in obstime
        to improve precision when the time resolution is wider than
        the time range or only a single value is requested.
        '''
        mjd_scaled = np.ravel(obstime.mjd / self. mjd_resolution)

        # ndmin=1 needed in case only a single coordinate is transformed
        mjd_lower = np.array(np.floor(mjd_scaled), ndmin=1, copy=False)
        mjd_upper = np.array(np.ceil(mjd_scaled), ndmin=1, copy=False)

        # unique already does sorting
        mjd_u = np.unique(np.concatenate([
            [mjd_scaled.min(), mjd_scaled.max()],
            mjd_lower, mjd_upper,
        ]))

        return Time(
            mjd_u * self.mjd_resolution,
            format='mjd',
            scale=obstime.scale,
        )

    @staticmethod
    def _prepare_earth_position_vel(support, obstime):
        pv_support, heliocentric_support = prepare_earth_position_vel(support)

        # do interpolation
        earth_pv = np.empty(obstime.shape, dtype=erfa.dt_pv)
        earth_heliocentric = np.empty(obstime.shape + (3,))
        for dim in range(3):
            for key in 'pv':
                earth_pv[key][..., dim] = np.interp(
                    obstime.mjd,
                    support.mjd,
                    pv_support[key][..., dim]
                )
            earth_heliocentric[..., dim] = np.interp(
                obstime.mjd, support.mjd, heliocentric_support[..., dim]
            )

        return earth_pv, earth_heliocentric

    @staticmethod
    def _get_cip(support, obstime):
        jd1_tt_support, jd2_tt_support = get_jd12(support, 'tt')
        cip_support = get_cip(jd1_tt_support, jd2_tt_support)
        return tuple(
            np.interp(obstime.mjd, support.mjd, cip_component)
            for cip_component in cip_support
        )

    def apci(self, frame):
        obstime = frame.obstime
        support = self._get_support_points(obstime)

        jd1_tt, jd2_tt = get_jd12(obstime, 'tt')
        cip = self._get_cip(support, obstime)
        earth_pv, earth_heliocentric = self._prepare_earth_position_vel(support, obstime)

        astrom = erfa.apci(jd1_tt, jd2_tt, earth_pv, earth_heliocentric, *cip)
        return astrom

    def apcs(self, frame):
        obstime = frame.obstime
        support = self._get_support_points(obstime)

        jd1_tt, jd2_tt = get_jd12(obstime, 'tt')
        # get the position and velocity arrays for the observatory.  Need to
        # have xyz in last dimension, and pos/vel in one-but-last.
        earth_pv, earth_heliocentric = self._prepare_earth_position_vel(support, obstime)
        pv = pav2pv(
            frame.obsgeoloc.get_xyz(xyz_axis=-1).value,
            frame.obsgeovel.get_xyz(xyz_axis=-1).value
        )
        return erfa.apcs(jd1_tt, jd2_tt, pv, earth_pv, earth_heliocentric)


class astrom_provider(ScienceState):
    """
    ScienceState to select with astrom provider is used in
    coordinate transformations.
    """

    _value = AstromProvider()

    @classmethod
    def validate(cls, value):
        if not isinstance(value, BaseAstromProvider):
            raise TypeError(f'Must be an instance of {BaseAstromProvider!r}')
        return value
