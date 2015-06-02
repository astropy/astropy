# Adapted from Astropy core package until 1.0 is released
#
# Copyright (c) 2011-2014, Astropy Developers
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of the Astropy Team nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import numpy as np

__doctest_skip__ = ['wcs_to_celestial_frame']


def _wcs_to_celestial_frame_builtin(wcs):

    from astropy.coordinates import FK4, FK4NoETerms, FK5, ICRS, Galactic
    from astropy.time import Time
    from astropy.wcs import WCSSUB_CELESTIAL

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

    Parameters
    ----------
    wcs : :class:`~astropy.wcs.WCS` instance
        The WCS to find the frame for

    Returns
    -------
    frame : :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame` subclass instance
        An instance of a :class:`~astropy.coordinates.baseframe.BaseCoordinateFrame`
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
