# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The following are private functions. These functions are registered into
:meth:`~astropy.cosmology.Cosmology.to_format` and
:meth:`~astropy.cosmology.Cosmology.from_format` and should only be accessed
via these methods.
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
        itis = isinstance(args[1], Cosmology) and (format in (None, "astropy.cosmology"))
    return itis


# ===================================================================
# Register

convert_registry.register_reader("astropy.cosmology", Cosmology, from_cosmology)
convert_registry.register_writer("astropy.cosmology", Cosmology, to_cosmology)
convert_registry.register_identifier("astropy.cosmology", Cosmology, cosmology_identify)
