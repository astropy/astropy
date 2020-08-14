# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains a helper function to fill erfa.astrom struct and a
ScienceState, which allows to speed up coordinate transformations at the
expense of accuracy.
"""
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

__all__ = ["astrom_interpolation_resolution", "get_astrom"]


class astrom_interpolation_resolution(ScienceState):
    """
    Time resolution for the support points used to interpolate
    slow changing astrom members.
    """

    _value = 0.0 * u.s

    @classmethod
    def validate(cls, value):
        value = u.Quantity(value, copy=False)
        return value.to(u.s)


@u.quantity_input(interpolation_resolution=u.s)
def get_astrom(frame, tcode, interpolation_resolution=None):
    """
    TBD
    """
    assert tcode in ['apio13', 'apci', 'apcs', 'apci13', 'apco13']

    if interpolation_resolution is None:
        interpolation_resolution = astrom_interpolation_resolution.get()

    if interpolation_resolution.to_value(u.ms) < 1:
        # below millisecond MJD resolution, one is probably better of
        # with no MJD binning
        mjd_resolution = None

    else:
        mjd_resolution = interpolation_resolution.to_value(u.day)

    obstime = frame.obstime
    jd1_tt, jd2_tt = get_jd12(obstime, 'tt')

    if mjd_resolution:
        # compute mjd support points for interpolation of Earth pv and cip
        mjd_scaled = obstime.mjd / mjd_resolution
        mjd_lower = np.array(np.floor(mjd_scaled), ndmin=1, copy=False)
        mjd_upper = np.array(np.ceil(mjd_scaled), ndmin=1, copy=False)

        mjd_u = np.unique(np.concatenate([
            [mjd_scaled.min(), mjd_scaled.max()],
            mjd_lower, mjd_upper,
        ]))  # does sorting

        obstime_support = Time(mjd_u * mjd_resolution, format='mjd', scale=obstime.scale)

    if tcode in ['apci', 'apcs']:
        # find the position and velocity of earth
        if mjd_resolution:
            (
                earth_pv_support, earth_heliocentric_support
                ) = prepare_earth_position_vel(obstime_support)
            # do interpolation
            earth_pv = np.empty(obstime.shape, dtype=erfa.dt_pv)
            earth_heliocentric = np.empty(obstime.shape + (3,))
            for dim in range(3):
                for key in 'pv':
                    earth_pv[key][..., dim] = np.interp(
                        obstime.mjd, obstime_support.mjd,
                        earth_pv_support[key][:, dim]
                    )
                earth_heliocentric[..., dim] = np.interp(
                    obstime.mjd, obstime_support.mjd, earth_heliocentric_support[:, dim]
                    )

        else:
            earth_pv, earth_heliocentric = prepare_earth_position_vel(obstime)

    if tcode == 'apio13':

        lon, lat, height = frame.location.to_geodetic('WGS84')

        xp, yp = get_polar_motion(obstime)
        jd1_utc, jd2_utc = get_jd12(obstime, 'utc')
        dut1utc = get_dut1utc(obstime)
        astrom = erfa.apio13(
            jd1_utc, jd2_utc, dut1utc,
            lon.to(u.radian).value, lat.to(u.radian).value,
            height.to(u.m).value,
            xp, yp,  # polar motion
            # all below are already in correct units because they are QuantityFrameAttribues
            frame.pressure.value,
            frame.temperature.value,
            frame.relative_humidity.value,
            frame.obswl.value,
        )

    elif tcode == 'apci':

        if mjd_resolution:
            jd1_tt_support, jd2_tt_support = get_jd12(obstime_support, 'tt')
            cip_support = get_cip(jd1_tt_support, jd2_tt_support)
            # do interpolation; cip is an ordinary array
            cip = tuple(
                np.interp(obstime.mjd, obstime_support.mjd, cip_support[i])
                for i in range(3)
                )

        else:
            cip = get_cip(jd1_tt, jd2_tt)

        astrom = erfa.apci(jd1_tt, jd2_tt, earth_pv, earth_heliocentric, *cip)

    elif tcode == 'apcs':

        # get the position and velocity arrays for the observatory.  Need to
        # have xyz in last dimension, and pos/vel in one-but-last.
        # (Note could use np.stack once our minimum numpy version is >=1.10.)
        pv = pav2pv(
            frame.obsgeoloc.get_xyz(xyz_axis=-1).value,
            frame.obsgeovel.get_xyz(xyz_axis=-1).value
        )
        astrom = erfa.apcs(jd1_tt, jd2_tt, pv, earth_pv, earth_heliocentric)

    elif tcode == 'apci13':
        pass
    elif tcode == 'apco13':
        pass

    return astrom
