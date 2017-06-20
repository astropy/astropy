# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..frame_attributes import TimeFrameAttribute

from ..baseframe import BaseRADecFrame
from .utils import DEFAULT_OBSTIME

class CIRS(BaseRADecFrame):
    """
    A coordinate or frame in the Celestial Intermediate Reference System (CIRS).

    This frame has one frame attribute:

    * ``obstime``
        The time at which the observation is taken.  Used for determining the
        position of the Earth and its precession.
    """

    obstime = TimeFrameAttribute(default=DEFAULT_OBSTIME)

CIRS.__doc__ += BaseRADecFrame.__doc__

# The "self-transform" is defined in icrs_cirs_transformations.py, because in
# the current implementation it goes through ICRS (like GCRS)
