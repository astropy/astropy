# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

__doctest_skip__ = ['wcs_to_celestial_frame']


def add_stokes_axis_to_wcs(wcs, add_before_ind):
    """
    Add a new Stokes axis that is uncorrelated with any other axes.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        The WCS to add to
    add_before_ind : int
        Index of the WCS to insert the new Stokes axis in front of.
        To add at the end, do add_before_ind = wcs.wcs.naxis
        The beginning is at position 0.

    Returns
    -------
    A new `~astropy.wcs.WCS` instance with an additional axis
    """

    inds = [i + 1 for i in range(wcs.wcs.naxis)]
    inds.insert(add_before_ind, 0)
    newwcs = wcs.sub(inds)
    newwcs.wcs.ctype[add_before_ind] = 'STOKES'
    newwcs.wcs.cname[add_before_ind] = 'STOKES'
    return newwcs


def _wcs_to_celestial_frame_builtin(wcs):

    from ..coordinates import FK4, FK4NoETerms, FK5, ICRS, Galactic
    from ..time import Time
    from . import WCSSUB_CELESTIAL

    # Keep only the celestial part of the axes
    wcs = wcs.sub([WCSSUB_CELESTIAL])

    if wcs.wcs.lng == -1 or wcs.wcs.lat == -1:
        return None

    radesys = wcs.wcs.radesys

    if np.isnan(wcs.wcs.equinox):
        equinox = None
    else:
        equinox = wcs.wcs.equinox

    xcoord = wcs.wcs.ctype[0][:4]
    ycoord = wcs.wcs.ctype[1][:4]

    # Apply logic from FITS standard to determine the default radesys
    if radesys == '' and xcoord == 'RA--' and ycoord == 'DEC-':
        if equinox is None:
            radesys = "ICRS"
        elif equinox < 1984.:
            radesys = "FK4"
        else:
            radesys = "FK5"

    if radesys == 'FK4':
        if equinox is not None:
            equinox = Time(equinox, format='byear')
        frame = FK4(equinox=equinox)
    elif radesys == 'FK4-NO-E':
        if equinox is not None:
            equinox = Time(equinox, format='byear')
        frame = FK4NoETerms(equinox=equinox)
    elif radesys == 'FK5':
        if equinox is not None:
            equinox = Time(equinox, format='jyear')
        frame = FK5(equinox=equinox)
    elif radesys == 'ICRS':
        frame = ICRS()
    else:
        if xcoord == 'GLON' and ycoord == 'GLAT':
            frame = Galactic()
        else:
            frame = None

    return frame


WCS_FRAME_MAPPINGS = [[_wcs_to_celestial_frame_builtin]]


class custom_frame_mappings(object):
    def __init__(self, mappings=[]):
        if hasattr(mappings, '__call__'):
            mappings = [mappings]
        WCS_FRAME_MAPPINGS.append(mappings)

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        WCS_FRAME_MAPPINGS.pop()


def wcs_to_celestial_frame(wcs):
    """
    For a given WCS, return the coordinate frame that matches the celestial
    component of the WCS.

    Paramters
    ---------
    wcs : :class:`~astropy.wcs.WCS` instance
        The WCS to find the frame for

    Returns
    -------
    frame : :class:`~astropy.coordinates.base_frame.BaseFrame` subclass instance
        An instance of a :class:`~astropy.coordinates.base_frame.BaseFrame`
        subclass instance that best matches the specified WCS.

    Notes
    -----

    To extend this function to frames not defined in astropy.coordinates, you
    can write your own function which should take a :class:`~astropy.wcs.WCS`
    instance and should return either an instance of a frame, or `None` if no
    matching frame was found. You can register this function temporarily with::

    >>> from astropy.wcs.utils import wcs_to_celestial_frame, custom_frame_mappings
    >>> with custom_frame_mappings(my_function):
    ...     wcs_to_celestial_frame(...)
    """
    for mapping_set in WCS_FRAME_MAPPINGS:
        for func in mapping_set:
            frame = func(wcs)
            if frame is not None:
                return frame
    raise ValueError("Could not determine celestial frame corresponding to "
                     "the specified WCS object")
