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
import warnings
from contextlib import contextmanager

from astropy.utils import find_current_module
from astropy.utils.decorators import deprecated

# Hack to make circular imports with units work
from astropy import units
del units

# These lines import some namespaces into the top level
from .constant import Constant, EMConstant  # noqa
from . import si  # noqa
from . import cgs  # noqa
from .config import codata, iaudata  # noqa

from . import utils as _utils  # noqa

# for updating the constants module docstring
_lines = [
    'The following constants are available:\n',
    '========== ============== ================ =========================',
    '   Name        Value            Unit       Description',
    '========== ============== ================ =========================',
]

# Catch warnings about "already has a definition in the None system"
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'Constant .*already has a definition')
    _utils._set_c(codata, iaudata, find_current_module(),
                  not_in_module_only=True, doclines=_lines, set_class=True)

_lines.append(_lines[1])

if __doc__ is not None:
    __doc__ += '\n'.join(_lines)


@deprecated('4.0', alternative="Use 'astropy.physical_constants' and 'astropy.astronomical_constants'")  # noqa
@contextmanager
def set_enabled_constants(modname):
    """
    Context manager to temporarily set values in the ``constants``
    namespace to an older version.
    See :ref:`astropy-constants-prior` for usage.

    Parameters
    ----------
    modname : {'astropyconst13', 'astropyconst20'}
        Name of the module containing an older version.

    """

    # Re-import here because these were deleted from namespace on init.
    import importlib
    import warnings
    from astropy.utils import find_current_module
    from . import utils as _utils

    try:
        modmodule = importlib.import_module('.constants.' + modname, 'astropy')
        codata_context = modmodule.codata
        iaudata_context = modmodule.iaudata
    except ImportError as exc:
        exc.args += ('Context manager does not currently handle {}'
                     .format(modname),)
        raise

    module = find_current_module()

    # Ignore warnings about "Constant xxx already has a definition..."
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',
                                'Constant .*already has a definition')
        _utils._set_c(codata_context, iaudata_context, module,
                      not_in_module_only=False, set_class=True)

    try:
        yield
    finally:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',
                                    'Constant .*already has a definition')
            _utils._set_c(codata, iaudata, module,
                          not_in_module_only=False, set_class=True)


# Clean up namespace
del find_current_module
del deprecated
del warnings
del contextmanager
del _utils
del _lines
