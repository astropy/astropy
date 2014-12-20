# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains functions/values used repeatedly in different modules of
the ``builtin_frames`` package.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import warnings

import numpy as np

from ... import units as u
from ...time import Time
from ...utils import iers
from ...utils.exceptions import AstropyWarning

# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
EQUINOX_J2000 = Time('J2000', scale='utc')
EQUINOX_B1950 = Time('B1950', scale='tai')

# This is a time object that is the default "obstime" when such an attribute is
# necessary.  Currently, we use J2000.
DEFAULT_OBSTIME = Time('J2000', scale='utc')

PIOVER2 = np.pi / 2.

#comes from the mean of the 1962-2014 IERS B data
_DEFAULT_PM = ((0.035, 0.29)*u.arcsec).to(u.radian).value

_IERS_HINT = """
If you need enough precision such that this matters (~<10 arcsec), you can download the latest IERS predictions by doing:
from astropy.utils.data import download_file
from astropy.utils import iers
iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))"""


def get_polar_motion(time):
    """
    gets the two polar motion components in radians for use with apio13
    """
    #get the polar motion from the IERS table
    xp, yp, status = iers.IERS.open().pm_xy(time.jd1, time.jd2, return_status=True)

    wmsg = None
    if status == iers.TIME_BEFORE_IERS_RANGE:
        wmsg = ('Tried to get polar motions for a time before IERS data is '
                'valid. Defaulting to polar motion from the 50-yr mean.')
    elif status == iers.TIME_BEYOND_IERS_RANGE:
        wmsg = ('Tried to get polar motions for a time after IERS data is '
                'valid. Defaulting to polar motion from the 50-yr mean.' + _IERS_HINT)
    if wmsg is not None:
        warnings.warn(wmsg, AstropyWarning)
        pms = np.repeat(_DEFAULT_PM, np.prod(time.shape))
        return pms.reshape([2] + list(time.shape))
    else:
        return xp.to(u.radian).value, yp.to(u.radian).value


def get_dut1utc(time):
    """
    This function is used to get UT1-UTC in coordinates because normally it
    gives an error outside the IERS range, but in coordinates we want to allow
    it to go through but with a warning.
    """
    try:
        return time.delta_ut1_utc
    except IndexError as e:
        msg = e.args[0] + ' Assuming UT1-UTC=0 for coordinate transformations.' + _IERS_HINT
        warnings.warn(msg, AstropyWarning)
        return np.zeros(time.shape)
