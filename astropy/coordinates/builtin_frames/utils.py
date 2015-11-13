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
from ..representation import CartesianRepresentation

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
_DEFAULT_PM = (0.035, 0.29)*u.arcsec

_IERS_HINT = """
If you need enough precision such that this matters (~<10 arcsec), you can
use the latest IERS predictions by running:

    >>> from astropy.utils import iers
    >>> iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL)

"""


def cartrepr_from_matmul(pmat, coo, transpose=False):
    """
    Note that pmat should be an ndarray, *not* a matrix.
    """
    if pmat.shape[-2:] != (3, 3):
        raise ValueError("tried to do matrix multiplication with an array that "
                         "doesn't end in 3x3")
    if coo.isscalar:
        # a simpler path for scalar coordinates
        if transpose:
            pmat = pmat.T
        newxyz = np.sum(pmat * coo.cartesian.xyz, axis=-1)
    else:
        xyz = coo.cartesian.xyz.T
        # these expression are the same as iterating over the first dimension of
        # pmat and xyz and doing matrix multiplication on each in turn.  resulting
        # dimension is <coo shape> x 3
        pmat = pmat.reshape(pmat.size//9, 3, 3)
        if transpose:
            pmat = pmat.transpose(0, 2, 1)
        newxyz = np.sum(pmat * xyz.reshape(xyz.size//3, 1, 3), axis=-1).T

    return CartesianRepresentation(newxyz)


def get_polar_motion(time):
    """
    gets the two polar motion components in radians for use with apio13
    """
    #get the polar motion from the IERS table
    xp, yp, status = iers.IERS.open().pm_xy(time, return_status=True)

    wmsg = None
    if np.any(status == iers.TIME_BEFORE_IERS_RANGE):
        wmsg = ('Tried to get polar motions for times before IERS data is '
                'valid. Defaulting to polar motion from the 50-yr mean for those.')
        xp.ravel()[status.ravel()==iers.TIME_BEFORE_IERS_RANGE] = _DEFAULT_PM[0]
        yp.ravel()[status.ravel()==iers.TIME_BEFORE_IERS_RANGE] = _DEFAULT_PM[1]

        warnings.warn(wmsg, AstropyWarning)

    if np.any(status == iers.TIME_BEYOND_IERS_RANGE):
        wmsg = ('Tried to get polar motions for times after IERS data is '
                'valid. Defaulting to polar motion from the 50-yr mean for those.' + _IERS_HINT)

        xp.ravel()[status.ravel()==iers.TIME_BEYOND_IERS_RANGE] = _DEFAULT_PM[0]
        yp.ravel()[status.ravel()==iers.TIME_BEYOND_IERS_RANGE] = _DEFAULT_PM[1]

        warnings.warn(wmsg, AstropyWarning)

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

def get_jd12(time, scale):
    """
    Gets ``jd1`` and ``jd2`` from a time object in a particular scale.

    Parameters
    ----------
    time : `~astropy.time.Time`
        The time to get the jds for
    scale : str
        The time scale to get the jds for

    Returns
    -------
    jd1 : float
    jd2 : float
    """
    if time.scale == scale:
        newtime = time
    elif time.scale == 'ut1' or scale == 'ut1':
        olddt = time.delta_ut1_utc
        time.delta_ut1_utc = get_dut1utc(time)
        newtime = getattr(time, scale)
        time.delta_ut1_utc =  olddt  # ensures no changes to the input `time`
    else:
        newtime = getattr(time, scale)

    return newtime.jd1, newtime.jd2

