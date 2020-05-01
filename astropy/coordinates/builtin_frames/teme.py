# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A coordinate or frame in the TEME (True Equator, Mean Equinox) system.

TEME is a True equator, Mean Equinox coordinate frame used in NORAD TLE
satellite files.
"""
# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.decorators import format_doc
from astropy.coordinates.representation import CartesianRepresentation, CartesianDifferential
from astropy.coordinates.baseframe import BaseCoordinateFrame, base_doc
from astropy.coordinates.attributes import TimeAttribute
from .utils import DEFAULT_OBSTIME

__all__ = ['TEME']

doc_footer_teme = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the frame is defined.  Used for determining the
        position of the Earth.
"""


@format_doc(base_doc, components="", footer=doc_footer_teme)
class TEME(BaseCoordinateFrame):
    """
    A coordinate or frame in the True Equator Mean Equinox frame (TEME).

    This frame is a geocentric system similar to CIRS or geocentric apparent place,
    except that the mean sidereal time is used to rotate from TIRS. TEME coordinates
    are most often used in combination with orbital data for satellites in the
    two-line-ephemeris format.

    Different implementations of the TEME frame exist. For clarity, this frame follows the
    conventions and relations to other frames that are set out in Vallado et al (2006).

    For more background on TEME, see the references provided in the
    :ref:`astropy-coordinates-seealso` section of the documentation.
    """

    default_representation = CartesianRepresentation
    default_differential = CartesianDifferential

    obstime = TimeAttribute()

# Transformation functions for getting to/from TEME and ITRS are in
# intermediate rotation transforms.py
