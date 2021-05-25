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


def _all_subclasses(cls):
    """Yield a (qualname, cls) of all subclasses (inclusive)."""
    # modified from https://stackoverflow.com/a/33607093
    yield cls.__qualname__, cls
    for subclass in cls.__subclasses__():
        yield from _all_subclasses(subclass)  # recursion in subclass


# Pre-defined cosmologies. This loops over the parameter sets in the
# parameters module and creates a cosmology instance with the same name as the
# parameter set in the current module's namespace.
_subclasses = dict(_all_subclasses(Cosmology))
for key in parameters.available:
    params = getattr(parameters, key)

    # TODO! this will need refactoring again when: parameter I/O is JSON/ECSSV
    cosmo_cls_name = params.pop("cosmology", None)

    if cosmo_cls_name in _subclasses:
        cosmo_cls = _subclasses[cosmo_cls_name]

        par = dict()
        meta = params.pop("meta", None) or {}
        for k, v in params.items():
            if k not in cosmo_cls._init_signature.parameters:
                meta.setdefault(k, v)  # merge into meta w/out overwriting
            elif k == "H0":
                par["H0"] = u.Quantity(v, u.km / u.s / u.Mpc)
            elif k == "m_nu":
                par["m_nu"] = u.Quantity(v, u.eV)
            else:
                par[k] = v

        ba = cosmo_cls._init_signature.bind_partial(**par, name=key, meta=meta)
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        cosmo.__doc__ = (f"{key} instance of {cosmo_cls.__name__} cosmology\n"
                         f"(from {meta['reference']})")

        setattr(sys.modules[__name__], key, cosmo)

# don't leave these variables floating around in the namespace
        del cosmo_cls, par, k, v, ba, cosmo
del key, params, cosmo_cls_name

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
