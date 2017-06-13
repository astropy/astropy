# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines deprecated units.

These units are not available in the top-level `astropy.units`
namespace. To use these units, you must import the `astropy.units.deprecated`
module::

    >>> from astropy.units import deprecated
    >>> q = 10. * deprecated.emu  # doctest: +SKIP

To include them in `~astropy.units.UnitBase.compose` and the results of
`~astropy.units.UnitBase.find_equivalent_units`, do::

    >>> from astropy.units import deprecated
    >>> deprecated.enable()  # doctest: +SKIP

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

_ns = globals()


def _initialize_module():
    # Local imports to avoid polluting top-level namespace
    from . import cgs
    from . import astrophys
    from .core import def_unit, _add_prefixes

    def_unit(['emu'], cgs.Bi, namespace=_ns,
             doc='Biot: CGS (EMU) unit of current')

    # add only some *prefixes* as deprecated units
    _add_prefixes(astrophys.solMass, namespace=_ns, prefixes=True)


_initialize_module()


###########################################################################
# DOCSTRING

# This generates a docstring for this module that describes all of the
# standard units defined here.
from .utils import generate_unit_summary as _generate_unit_summary  # noqa
if __doc__ is not None:
    __doc__ += _generate_unit_summary(globals())


def enable():
    """
    Enable deprecated units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.

    This may be used with the ``with`` statement to enable deprecated
    units only temporarily.
    """
    # Local import to avoid cyclical import
    from .core import add_enabled_units
    # Local import to avoid polluting namespace
    import inspect
    return add_enabled_units(inspect.getmodule(enable))
