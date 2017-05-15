# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains a helper function to fill erfa.astrom struct and a
ScienceState, which allows to speed up coordinate transformations at the
expense of accuracy.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import OrderedDict

import numpy as np

from ..time import Time
from .sky_coordinate import SkyCoord
from ..utils.decorators import classproperty
from ..utils.state import ScienceState
# from ..utils import indent
from .. import units as u
from .. import _erfa as erfa
from ..extern import six
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

    if mjd_resolution:
        # apply binning to obstime
        mjd_binned = np.int64(obstime.mjd / mjd_resolution + 0.5)
        mjd_u, mjd_idx, mjd_uidx = np.unique(
            mjd_binned, return_index=True, return_inverse=True
            )

        # obstime_binned = obstime[mjd_idx]
        obstime_binned = Time(mjd_u * mjd_resolution, format='mjd', scale=obstime.scale)
    else:
        obstime_binned = obstime

    if tcode in ['apci', 'apcs']:
        # find the position and velocity of earth
        jd1_tt, jd2_tt = get_jd12(obstime_binned, 'tt')
        earth_pv, earth_heliocentric = prepare_earth_position_vel(obstime_binned)

    if tcode == 'apio13':

        lon, lat, height = frame.location.to_geodetic('WGS84')

        xp, yp = get_polar_motion(obstime_binned)
        jd1_utc, jd2_utc = get_jd12(obstime_binned, 'utc')
        dut1utc = get_dut1utc(obstime_binned)
        astrom_binned = erfa.apio13(
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

        x, y, s = get_cip(jd1_tt, jd2_tt)
        astrom_binned = erfa.apci(jd1_tt, jd2_tt, earth_pv, earth_heliocentric, x, y, s)

    elif tcode == 'apcs':

        # get the position and velocity arrays for the observatory.  Need to
        # have xyz in last dimension, and pos/vel in one-but-last.
        # (Note could use np.stack once our minimum numpy version is >=1.10.)
        pv = np.concatenate(
                (frame.obsgeoloc.get_xyz(xyz_axis=-1).value[..., np.newaxis, :],
                 frame.obsgeovel.get_xyz(xyz_axis=-1).value[..., np.newaxis, :]),
                axis=-2)
        astrom_binned = erfa.apcs(jd1_tt, jd2_tt, pv, earth_pv, earth_heliocentric)

    elif tcode == 'apci13':
        pass
    elif tcode == 'apco13':
        pass

    if mjd_resolution:
        astrom = astrom_binned[mjd_uidx]
        jd1_ut1, jd2_ut1 = get_jd12(obstime, 'ut1')
        astrom = erfa.aper13(jd1_ut1, jd2_ut1, astrom)
    else:
        astrom = astrom_binned

    return astrom

