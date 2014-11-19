
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import inspect

# Dependencies
import numpy as np

# Project
from ...extern import six
from ...utils.compat.odict import OrderedDict
from ... import units as u
from ...time import Time
from ..angles import Angle
from ..representation import (SphericalRepresentation, CartesianRepresentation,
                             UnitSphericalRepresentation)
from ..baseframe import (BaseCoordinateFrame, frame_transform_graph, GenericFrame,
                        FrameAttribute, TimeFrameAttribute,
                        RepresentationMapping)
from ..transformations import FunctionTransform, DynamicMatrixTransform


# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
_EQUINOX_J2000 = Time('J2000', scale='utc')
_EQUINOX_B1950 = Time('B1950', scale='tai')



class AltAz(BaseCoordinateFrame):
    """
    A coordinate or frame in the Altitude-Azimuth system (i.e., Horizontal
    coordinates).

    .. warning::
        The AltAz class currently does not support any transformations. In a
        future version, it will support the standard IAU2000 AltAz<->ICRS
        transformations.  It is provided right now as a placeholder for storing
        as-observed horizontal coordinates.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    az : `Angle`, optional, must be keyword
        The Azimuth for this object (``alt`` must also be given and
        ``representation`` must be None).
    alt : `Angle`, optional, must be keyword
        The Altitude for this object (``az`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'az'),
                      RepresentationMapping('lat', 'alt')],
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=_EQUINOX_B1950)
    location = FrameAttribute(default=None)
    obstime = TimeFrameAttribute(default=None)

    def __init__(self, *args, **kwargs):
        from warnings import warn
        from astropy.utils.exceptions import AstropyWarning

        warn(AstropyWarning('The AltAz class currently does not support any '
                            'transformations.  In a future version, it will '
                            'support the standard IAU2000 AltAz<->ICRS '
                            'transformations.'))
        super(AltAz, self).__init__(*args, **kwargs)
