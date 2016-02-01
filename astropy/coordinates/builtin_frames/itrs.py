# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import CartesianRepresentation
from ..baseframe import BaseCoordinateFrame, TimeFrameAttribute
from .utils import DEFAULT_OBSTIME


class ITRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the International Terrestrial Reference System
    (ITRS).  This is approximately a geocentric system, although strictly it is
    defined by a series of reference locations near the surface of the Earth.
    For more background on the ITRS, see the references provided in the
    :ref:`astropy-coordinates-seealso` section of the documentation.
    """

    default_representation = CartesianRepresentation

    obstime = TimeFrameAttribute(default=DEFAULT_OBSTIME)

    @property
    def earth_location(self):
        """
        The data in this frame as an `~astropy.coordinates.EarthLocation` class.
        """
        from ..earth import EarthLocation

        cart = self.represent_as(CartesianRepresentation)
        return EarthLocation(x=cart.x, y=cart.y, z=cart.z)

# Self-transform is in intermediate_rotation_transforms.py with all the other
# ITRS transforms
