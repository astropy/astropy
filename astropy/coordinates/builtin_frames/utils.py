# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains functions/values used repeatedly in different modules of
the ``builtin_frames`` package.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ...time import Time
from ...utils import iers

# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
EQUINOX_J2000 = Time('J2000', scale='utc')
EQUINOX_B1950 = Time('B1950', scale='tai')

# This is a time object that is the default "obstime" when such an attribute is
# necessary.  Currently, we use J2000.
DEFAULT_OBSTIME = Time('J2000', scale='utc')

PIOVER2 = np.pi / 2.


def get_polar_motion(time):
    """
    gets the two polar motion components in radians for use with apio13
    """
    #get the polar motion from the IERS table
    xp, yp = iers.IERS.open().pm_xy(time.jd1, time.jd2)
    return xp.to(u.radian).value, yp.to(u.radian).value
