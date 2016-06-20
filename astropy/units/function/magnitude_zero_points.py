# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines magnitude zero points.  By default, they are used to
define corresponding magnitudes, but not enabled as regular physical units.
To enable them, do::

    >>> from astropy.units import magnitude_zero_points
    >>> magnitude_zero_points.enable()  # doctest: +SKIP
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as _numpy
from ..core import UnitBase, def_unit

# One cannot import L_bol0 directly, or the order of imports of units and
# constants starts to matter on python2. [#5121]
from ...constants import si as _si
from .. import si, astrophys


_ns = globals()

def_unit(['Bol', 'L_bol'], _si.L_bol0, namespace=_ns, prefixes=False,
         doc="Luminosity corresponding to absolute bolometric magnitude zero")
def_unit(['bol', 'f_bol'], _si.L_bol0 / (4 * _numpy.pi * (10.*astrophys.pc)**2),
         namespace=_ns, prefixes=False, doc="Irradiance corresponding to "
         "appparent bolometric magnitude zero")
def_unit(['AB'], 10.**(-0.4*48.6) * 1.e-3 * si.W / si.m**2 / si.Hz,
         namespace=_ns, prefixes=False,
         doc="AB magnitude zero flux density.")
def_unit(['ST'], 10.**(-0.4*21.1) * 1.e-3 * si.W / si.m**2 / si.AA,
         namespace=_ns, prefixes=False,
         doc="ST magnitude zero flux density.")

###########################################################################
# CLEANUP

del UnitBase
del def_unit
del si
del astrophys

###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from ..utils import generate_unit_summary as _generate_unit_summary
if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())


def enable():
    """
    Enable magnitude zero point units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable these
    units only temporarily.
    """
    # Local import to avoid cyclical import
    from ..core import add_enabled_units
    # Local import to avoid polluting namespace
    import inspect
    return add_enabled_units(inspect.getmodule(enable))
