# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..attributes import TimeAttribute

from .baseradec import _base_radec_docstring, BaseRADecFrame
from .utils import DEFAULT_OBSTIME


class CIRS(BaseRADecFrame):
    """
    A coordinate or frame in the Celestial Intermediate Reference System (CIRS).

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Earth and its precession.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


CIRS.__doc__ = CIRS.__doc__.format(params=_base_radec_docstring)

# The "self-transform" is defined in icrs_cirs_transformations.py, because in
# the current implementation it goes through ICRS (like GCRS)
