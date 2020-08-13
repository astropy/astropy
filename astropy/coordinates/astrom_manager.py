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
    get_jd12, get_cip, prepare_earth_position_vel, get_polar_motion, get_dut1utc
)

__all__ = ["transform_precision"]


DEFAULT_PRECISION = None  # no MJD binning


class transform_precision(ScienceState):
    """
    TBD
    """

    _value = 0.  # value steers the transform precision (MJD_BINNING)


def get_astrom(frame, tcode, precision=None):
    """
    TBD
    """

    assert tcode in ['apio13', 'apci', 'apcs', 'apci13', 'apco13']

    if precision is None:
        precision = transform_precision.get()

    if precision < 1.e-3:
        # below millisecond MJD resolution, one is probably better of
        # with no MJD binning
        mjd_resolution = None

    else:
        mjd_resolution = precision / 86400.  # in days

    obstime = frame.obstime
    jd1_tt, jd2_tt = get_jd12(obstime, 'tt')

    if mjd_resolution:
        # compute mjd support points for interpolation of Earth pv and cip
        mjd_lower = np.int64(obstime.mjd / mjd_resolution)
        mjd_upper = mjd_lower + 1
        mjd_u = np.unique([mjd_lower, mjd_upper])  # does sorting

        obstime_support = Time(mjd_u * mjd_resolution, format='mjd', scale=obstime.scale)

    if tcode in ['apci', 'apcs']:
        # find the position and velocity of earth
        if mjd_resolution:
            (
                earth_pv_support, earth_heliocentric_support
                ) = prepare_earth_position_vel(obstime_support)
            # do interpolation
            earth_pv = np.empty(obstime.shape + (2, 3,))
            earth_heliocentric = np.empty(obstime.shape + (3,))
            for dim in range(3):
                for ipv in range(2):
                    earth_pv[..., ipv, dim] = np.interp(
                        obstime.mjd, obstime_support.mjd, earth_pv_support[:, ipv, dim]
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
            frame.relative_humidity,
            frame.obswl.value
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
        pv = np.concatenate(
                (frame.obsgeoloc.get_xyz(xyz_axis=-1).value[..., np.newaxis, :],
                 frame.obsgeovel.get_xyz(xyz_axis=-1).value[..., np.newaxis, :]),
                axis=-2)
        astrom = erfa.apcs(jd1_tt, jd2_tt, pv, earth_pv, earth_heliocentric)

    elif tcode == 'apci13':
        pass
    elif tcode == 'apco13':
        pass

    return astrom
