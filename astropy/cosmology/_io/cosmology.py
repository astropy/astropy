# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""|Cosmology| I/O, using |Cosmology.to_format| and |Cosmology.from_format|.

This module provides functions to transform a |Cosmology| object to and from another
|Cosmology| object. The functions are registered with ``convert_registry`` under the
format name "astropy.cosmology". You probably won't need to use these functions as they
are present mainly for completeness and testing.

    >>> from astropy.cosmology import Cosmology, Planck18
    >>> Planck18.to_format("astropy.cosmology") is Planck18
    True
    >>> Cosmology.from_format(Planck18) is Planck18
    True
"""

from astropy.cosmology.connect import convert_registry
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Cosmology

__all__ = []  # nothing is publicly scoped


def from_cosmology(cosmo, /, **kwargs):
    """Return the |Cosmology| unchanged.

    Parameters
    ----------
    cosmo : `~astropy.cosmology.Cosmology`
        The cosmology to return.
    **kwargs
        This argument is required for compatibility with the standard set of
        keyword arguments in format `~astropy.cosmology.Cosmology.from_format`,
        e.g. "cosmology". If "cosmology" is included and is not `None`,
        ``cosmo`` is checked for correctness.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance
        Just ``cosmo`` passed through.

    Raises
    ------
    TypeError
        If the |Cosmology| object is not an instance of ``cosmo`` (and
        ``cosmology`` is not `None`).

    Examples
    --------
    >>> from astropy.cosmology import Cosmology, Planck18
    >>> print(Cosmology.from_format(Planck18))
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                  Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
    """
    # Check argument `cosmology`
    cosmology = kwargs.get("cosmology")
    if isinstance(cosmology, str):
        cosmology = _COSMOLOGY_CLASSES[cosmology]
    if cosmology is not None and not isinstance(cosmo, cosmology):
        raise TypeError(f"cosmology {cosmo} is not an {cosmology} instance.")

    return cosmo


def to_cosmology(cosmo, *args):
    """Return the |Cosmology| unchanged.

    Parameters
    ----------
    cosmo : `~astropy.cosmology.Cosmology`
        The cosmology to return.
    *args
        Not used.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance
        Just ``cosmo`` passed through.

    Examples
    --------
    >>> from astropy.cosmology import Planck18
    >>> Planck18.to_format("astropy.cosmology") is Planck18
    True
    """
    return cosmo


def cosmology_identify(origin, format, *args, **kwargs):
    """Identify if object is a `~astropy.cosmology.Cosmology`.

    Returns
    -------
    bool
    """
    itis = False
    if origin == "read":
        itis = isinstance(args[1], Cosmology) and (
            format in (None, "astropy.cosmology")
        )
    return itis


# ===================================================================
# Register

convert_registry.register_reader("astropy.cosmology", Cosmology, from_cosmology)
convert_registry.register_writer("astropy.cosmology", Cosmology, to_cosmology)
convert_registry.register_identifier("astropy.cosmology", Cosmology, cosmology_identify)
