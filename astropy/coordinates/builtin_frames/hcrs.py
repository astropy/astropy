# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import SphericalRepresentation
from ..baseframe import (BaseCoordinateFrame, RepresentationMapping,
                         TimeFrameAttribute)
from .utils import DEFAULT_OBSTIME


class HCRS(BaseCoordinateFrame):
    """
    A coordinate or frame in a Heliocentric system, with axes aligned to ICRS.

    The ICRS has an origin at the Barycenter and axes which are fixed with
    respect to space.

    This coordinate system is distinct from ICRS mainly in that it is relative
    to the Sun's center-of-mass rather than the solar system Barycenter.
    In principle, therefore, this frame should include the effects of
    aberration (unlike ICRS), but this is not done, since they are very small,
    of the order of 8 milli-arcseconds.

    For more background on the ICRS and related coordinate transformations, see
    the references provided in the :ref:`astropy-coordinates-seealso` section of
    the documentation.

    This frame has these frame attributes:

    * ``obstime``
        The time at which the observation is taken.  Used for determining the
        position of the Sun.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation or `None` to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        Right ascension for this object (``dec`` must also be given and
        ``representation`` must be `None`).
    dec : `Angle`, optional, must be keyword
        Declination for this object (``ra`` must also be given and
        ``representation`` must be `None`).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be `None`).
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    obstime = TimeFrameAttribute(default=DEFAULT_OBSTIME)
