# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..attributes import TimeAttribute
from .utils import DEFAULT_OBSTIME
from .baseradec import _base_radec_docstring, BaseRADecFrame


class HCRS(BaseRADecFrame):
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

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Sun.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


HCRS.__doc__ = HCRS.__doc__.format(params=_base_radec_docstring)

# Transformations are defined in icrs_circ_transforms.py
