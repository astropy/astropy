# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Built-in cosmologies.

See :attr:`~astropy.cosmology.realizations.available` for a full list.
"""

from __future__ import annotations

__all__ = [  # noqa: F822 (undefined name)
    "available",
    "default_cosmology",
    # Realizations (dynamic attribute, see __getattr__)
    "WMAP1",
    "WMAP3",
    "WMAP5",
    "WMAP7",
    "WMAP9",
    "Planck13",
    "Planck15",
    "Planck18",
]

import pathlib
import sys
from typing import TYPE_CHECKING

from astropy.utils.data import get_pkg_data_path
from astropy.utils.state import ScienceState

from . import _io  # Ensure IO methods are registered, to read realizations # noqa: F401
from .core import Cosmology

if TYPE_CHECKING:
    from typing import ClassVar

__doctest_requires__ = {"*": ["scipy"]}

_COSMOLOGY_DATA_DIR = pathlib.Path(
    get_pkg_data_path("cosmology", "data", package="astropy")
)
available = (
    "WMAP1",
    "WMAP3",
    "WMAP5",
    "WMAP7",
    "WMAP9",
    "Planck13",
    "Planck15",
    "Planck18",
)


def __getattr__(name: str) -> Cosmology:
    """Make specific realizations from data files with lazy import from ``PEP 562``.

    Raises
    ------
    AttributeError
        If "name" is not in :mod:`astropy.cosmology.realizations`
    """
    if name not in available:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}.")

    cosmo = Cosmology.read(
        str(_COSMOLOGY_DATA_DIR / name) + ".ecsv", format="ascii.ecsv"
    )
    object.__setattr__(
        cosmo,
        "__doc__",
        f"{name} instance of {cosmo.__class__.__qualname__} "
        f"cosmology\n(from {cosmo.meta['reference']})",
    )

    # Cache in this module so `__getattr__` is only called once per `name`.
    setattr(sys.modules[__name__], name, cosmo)

    return cosmo


def __dir__() -> list[str]:
    """Directory, including lazily-imported objects."""
    return __all__


#########################################################################
# The science state below contains the current cosmology.
#########################################################################


class default_cosmology(ScienceState):
    """The default cosmology to use.

    To change it::

        >>> from astropy.cosmology import default_cosmology, WMAP7
        >>> with default_cosmology.set(WMAP7):
        ...     # WMAP7 cosmology in effect
        ...     pass

    Or, you may use a string::

        >>> with default_cosmology.set('WMAP7'):
        ...     # WMAP7 cosmology in effect
        ...     pass

    To get the default cosmology:

        >>> default_cosmology.get()
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>,
                      Om0=0.30966, ...
    """

    _default_value: ClassVar[str] = "Planck18"
    _value: ClassVar[str | Cosmology] = "Planck18"

    @classmethod
    def validate(cls, value: Cosmology | str | None) -> Cosmology | None:
        """Return a Cosmology given a value.

        Parameters
        ----------
        value : None, str, or `~astropy.cosmology.Cosmology`

        Returns
        -------
        `~astropy.cosmology.Cosmology` instance

        Raises
        ------
        TypeError
            If ``value`` is not a string or |Cosmology|.
        """
        # None -> default
        if value is None:
            value = cls._default_value

        # Parse to Cosmology. Error if cannot.
        if isinstance(value, str):
            # special-case one string
            if value == "no_default":
                value = None
            else:
                value = cls._get_from_registry(value)
        elif not isinstance(value, Cosmology):
            raise TypeError(
                "default_cosmology must be a string or Cosmology instance, "
                f"not {value}."
            )

        return value

    @classmethod
    def _get_from_registry(cls, name: str) -> Cosmology:
        """Get a registered Cosmology realization.

        Parameters
        ----------
        name : str
            The built-in |Cosmology| realization to retrieve.

        Returns
        -------
        `astropy.cosmology.Cosmology`
            The cosmology realization of `name`.

        Raises
        ------
        ValueError
            If ``name`` is a str, but not for a built-in Cosmology.
        TypeError
            If ``name`` is for a non-Cosmology object.
        """
        try:
            value = getattr(sys.modules[__name__], name)
        except AttributeError:
            raise ValueError(
                f"Unknown cosmology {name!r}. Valid cosmologies:\n{available}"
            )

        if not isinstance(value, Cosmology):
            raise TypeError(f"cannot find a Cosmology realization called {name}.")

        return value
