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

from __future__ import annotations

from typing import TYPE_CHECKING

from astropy.cosmology._src.core import _COSMOLOGY_CLASSES, Cosmology
from astropy.cosmology.io import convert_registry

if TYPE_CHECKING:
    from astropy.cosmology._src.typing import _CosmoT

__all__: list[str] = []  # nothing is publicly scoped


def from_cosmology(
    cosmo: _CosmoT, /, cosmology: type[_CosmoT] | str | None = None, **kwargs: object
) -> _CosmoT:
    """Return the |Cosmology| unchanged.

    Parameters
    ----------
    cosmo : `~astropy.cosmology.Cosmology`
        The cosmology to return.
    cosmology : type[`~astropy.cosmology.Cosmology`] | str | None, optional
        The |Cosmology| class to check against. If not `None`, ``cosmo`` is checked
        for correctness.
    **kwargs
        This argument is required for compatibility with the standard set of
        keyword arguments in format |Cosmology.from_format|.

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
    # hot path
    if cosmology is None:
        return cosmo

    cosmo_cls = (
        _COSMOLOGY_CLASSES[cosmology] if isinstance(cosmology, str) else cosmology
    )
    if not isinstance(cosmo, cosmo_cls):
        raise TypeError(f"cosmology {cosmo} is not an {cosmo_cls} instance.")

    return cosmo


def to_cosmology(cosmo: _CosmoT, *args: object) -> _CosmoT:
    """Return the |Cosmology| unchanged.

    Parameters
    ----------
    cosmo : `~astropy.cosmology.Cosmology`
        The cosmology to return.
    *args : object
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


def cosmology_identify(
    origin: str, format: str | None, *args: object, **kwargs: object
) -> bool:
    """Identify if object is a `~astropy.cosmology.Cosmology`.

    This checks if the 2nd argument is a |Cosmology| instance and the format is
    "astropy.cosmology" or `None`.

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
