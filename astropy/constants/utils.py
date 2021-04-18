# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utility functions for ``constants`` sub-package."""
import itertools

__all__ = []


def _get_c(codata, iaudata, module, not_in_module_only=True):
    """
    Generator to return a Constant object.

    Parameters
    ----------
    codata, iaudata : obj
        Modules containing CODATA and IAU constants of interest.

    module : obj
        Namespace module of interest.

    not_in_module_only : bool
        If ``True``, ignore constants that are already in the
        namespace of ``module``.

    Returns
    -------
    _c : Constant
        Constant object to process.

    """
    from .constant import Constant

    for _nm, _c in itertools.chain(sorted(vars(codata).items()),
                                   sorted(vars(iaudata).items())):
        if not isinstance(_c, Constant):
            continue
        elif (not not_in_module_only) or (_c.abbrev not in module.__dict__):
            yield _c


def _set_c(codata, iaudata, module, not_in_module_only=True, doclines=None,
           set_class=False):
    """
    Set constants in a given module namespace.

    Parameters
    ----------
    codata, iaudata : obj
        Modules containing CODATA and IAU constants of interest.

    module : obj
        Namespace module to modify with the given ``codata`` and ``iaudata``.

    not_in_module_only : bool
        If ``True``, constants that are already in the namespace
        of ``module`` will not be modified.

    doclines : list or None
        If a list is given, this list will be modified in-place to include
        documentation of modified constants. This can be used to update
        docstring of ``module``.

    set_class : bool
        Namespace of ``module`` is populated with ``_c.__class__``
        instead of just ``_c`` from :func:`_get_c`.

    """
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
            doclines.append('{:^10} {:^14.9g} {:^16} {}'.format(
                _c.abbrev, _c.value, _c._unit_string, _c.name))
