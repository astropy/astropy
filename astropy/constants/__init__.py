# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains astronomical and physical constants for use in Astropy or other
places.

A typical use case might be::

    >>> from astropy.constants import c, m_e
    >>> # ... define the mass of something you want the rest energy of as m ...
    >>> m = m_e
    >>> E = m * c**2
    >>> E.to('MeV')  # doctest: +FLOAT_CMP
    <Quantity 0.510998927603161 MeV>

"""
import inspect
from contextlib import contextmanager

# Hack to make circular imports with units work
try:
    from .. import units
    del units
except ImportError:
    pass

from .constant import Constant, EMConstant  # noqa
from . import si  # noqa
from . import cgs  # noqa
from . import codata2014, iau2015  # noqa
from . import utils as _utils

# for updating the constants module docstring
_lines = [
    'The following constants are available:\n',
    '========== ============== ================ =========================',
    '   Name        Value            Unit       Description',
    '========== ============== ================ =========================',
]

# NOTE: Update this when default changes.
_utils._set_c(codata2014, iau2015, inspect.getmodule(inspect.currentframe()),
              not_in_module_only=True, doclines=_lines, set_class=True)

_lines.append(_lines[1])

if __doc__ is not None:
    __doc__ += '\n'.join(_lines)


# TODO: Re-implement in a way that is more consistent with astropy.units.
#       See https://github.com/astropy/astropy/pull/7008 discussions.
@contextmanager
def set_enabled_constants(modname):
    """
    Context manager to temporarily set values in the ``constants``
    namespace to an older version.
    See :ref:`astropy-constants-prior` for usage.

    Parameters
    ----------
    modname : {'astropyconst13'}
        Name of the module containing an older version.

    """

    # Re-import here because these were deleted from namespace on init.
    import inspect
    import warnings
    from . import utils as _utils

    # NOTE: Update this when default changes.
    if modname == 'astropyconst13':
        from .astropyconst13 import codata2010 as codata
        from .astropyconst13 import iau2012 as iaudata
    else:
        raise ValueError(
            'Context manager does not currently handle {}'.format(modname))

    module = inspect.getmodule(inspect.currentframe())

    # Ignore warnings about "Constant xxx already has a definition..."
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _utils._set_c(codata, iaudata, module,
                      not_in_module_only=False, set_class=True)

    try:
        yield
    finally:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # NOTE: Update this when default changes.
            _utils._set_c(codata2014, iau2015, module,
                          not_in_module_only=False, set_class=True)


# Clean up namespace
del inspect
del contextmanager
del _utils
del _lines
