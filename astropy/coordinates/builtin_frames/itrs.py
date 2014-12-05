# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import CartesianRepresentation
from ..baseframe import BaseCoordinateFrame


class ITRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the Intertaional Terrestrial Reference System
    (ITRS).  This is approximately a geocentric system, although strictly it is
    defined by a series of reference locations.
    """

    default_representation = CartesianRepresentation
