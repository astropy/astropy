# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import warnings

from astropy import units as u
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.utils.state import ScienceState

from . import parameters
from .core import Cosmology, FlatLambdaCDM, LambdaCDM

__all__ = ["default_cosmology"] + parameters.available

__doctest_requires__ = {"*": ["scipy"]}

# Pre-defined cosmologies. This loops over the parameter sets in the
# parameters module and creates a LambdaCDM or FlatLambdaCDM instance
# with the same name as the parameter set in the current module's namespace.

for key in parameters.available:
    par = getattr(parameters, key)

    if par["cosmology"] == "FlatLambdaCDM":
        cosmo = FlatLambdaCDM(
            par["H0"],
            par["Om0"],
            Tcmb0=par["Tcmb0"],
            Neff=par["Neff"],
            m_nu=u.Quantity(par["m_nu"], u.eV),
            name=key,
            Ob0=par["Ob0"],
        )
        docstr = "{} instance of FlatLambdaCDM cosmology\n\n(from {})"
        cosmo.__doc__ = docstr.format(key, par["reference"])
    else:
        warnings.warn("Please open a PR for your added cosmology realization.")
        # For a non-flat LCDM realization, the following is the code necessary
        # to create the realization, given the parameters. We comment this out
        # as there are not currently any such built-in realizations.
        # cosmo = LambdaCDM(
        #     par["H0"],
        #     par["Om0"],
        #     par["Ode0"],
        #     Tcmb0=par["Tcmb0"],
        #     Neff=par["Neff"],
        #     m_nu=u.Quantity(par["m_nu"], u.eV),
        #     name=key,
        #     Ob0=par["Ob0"],
        # )
        # docstr = "{} instance of LambdaCDM cosmology\n\n(from {})"
        # cosmo.__doc__ = docstr.format(key, par["reference"])

    setattr(sys.modules[__name__], key, cosmo)

# don't leave these variables floating around in the namespace
del key, par, cosmo

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
