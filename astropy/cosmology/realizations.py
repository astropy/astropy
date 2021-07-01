# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import warnings

from astropy import units as u
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.state import ScienceState

from . import parameters
from .core import Cosmology

__all__ = ["default_cosmology"] + parameters.available

__doctest_requires__ = {"*": ["scipy"]}


# Pre-defined cosmologies. This loops over the parameter sets in the
# parameters module and creates a cosmology instance with the same name as the
# parameter set in the current module's namespace.
for key in parameters.available:
    map = getattr(parameters, key)
    map["name"] = key
    cosmo = Cosmology.from_(map, format="mapping", move_to_meta=True)
    cosmo.__doc__ = (f"{key} instance of {cosmo.__class__.__name__} cosmology\n"
                     f"(from {cosmo.meta['reference']})")

    setattr(sys.modules[__name__], key, cosmo)

# don't leave these variables floating around in the namespace
del key, cosmo, map

#########################################################################
# The science state below contains the current cosmology.
#########################################################################


class default_cosmology(ScienceState):
    """
    The default cosmology to use.  To change it::

        >>> from astropy.cosmology import default_cosmology, WMAP7
        >>> with default_cosmology.set(WMAP7):
        ...     # WMAP7 cosmology in effect
        ...     pass

    Or, you may use a string::

        >>> with default_cosmology.set('WMAP7'):
        ...     # WMAP7 cosmology in effect
        ...     pass
    """

    _value = "Planck18"

    @staticmethod
    def get_cosmology_from_string(arg):
        """ Return a cosmology instance from a string.
        """
        if arg == "no_default":
            cosmo = None
        else:
            try:
                cosmo = getattr(sys.modules[__name__], arg)
            except AttributeError:
                s = "Unknown cosmology '{}'. Valid cosmologies:\n{}".format(
                    arg, parameters.available
                )
                raise ValueError(s)
        return cosmo

    @classmethod
    def validate(cls, value):
        if value is None:
            value = "Planck18"
        if isinstance(value, str):
            if value == "Planck18_arXiv_v2":
                warnings.warn(
                    f"{value} is deprecated in astropy 4.2, use Planck18 instead",
                    AstropyDeprecationWarning,
                )
            return cls.get_cosmology_from_string(value)
        elif isinstance(value, Cosmology):
            return value
        else:
            raise TypeError(
                "default_cosmology must be a string or Cosmology instance."
            )
