# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import CartesianRepresentation
from ..baseframe import BaseCoordinateFrame, TimeFrameAttribute, frame_transform_graph
from ..transformations import FunctionTransform
from .utils import DEFAULT_OBSTIME


class ITRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the International Terrestrial Reference System
    (ITRS).  This is approximately a geocentric system, although strictly it is
    defined by a series of reference locations near the surface of the Earth.
    For more background on the ITRS, see the references provided in the
    :ref:`astropy-coordinates-seealso` section of the documentation.

    Note that this frame does *not* have an ``obstime``, because it is
    rotating with the Earth.
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


@frame_transform_graph.transform(FunctionTransform, ITRS, ITRS)
def itrs_to_itrs(from_coo, to_frame):
    # this is a trivial self-transform because the frame rotates with the earth,
    # so obstime is only necessary for transforms
    return to_frame.realize_frame(from_coo.data)
