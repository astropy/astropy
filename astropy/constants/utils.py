# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utility functions for ``constants`` sub-package."""
import itertools

__all__ = []


def _get_c(codata, iaudata, module, not_in_module_only=True):
    from .constant import Constant

    for _nm, _c in itertools.chain(sorted(vars(codata).items()),
                                   sorted(vars(iaudata).items())):
        if not isinstance(_c, Constant):
            continue
        elif (not not_in_module_only) or (_c.abbrev not in module.__dict__):
            yield _c


def _set_c(codata, iaudata, module, not_in_module_only=True, doclines=None,
           set_class=False):

    for _c in _get_c(codata, iaudata, module,
                     not_in_module_only=not_in_module_only):
        if set_class:
            value = _c.__class__(_c.abbrev, _c.name, _c.value,
                                 _c._unit_string, _c.uncertainty,
                                 _c.reference)
        else:
            value = _c

        setattr(module, _c.abbrev, value)

        if doclines is not None:
            doclines.append('{0:^10} {1:^14.9g} {2:^16} {3}'.format(
                _c.abbrev, _c.value, _c._unit_string, _c.name))
