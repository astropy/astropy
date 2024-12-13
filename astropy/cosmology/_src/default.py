"""Default Cosmology."""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__ = ["default_cosmology"]


from typing import TYPE_CHECKING

from astropy.utils.state import ScienceState

# isort: split
from astropy.cosmology._src.core import Cosmology

if TYPE_CHECKING:
    from typing import ClassVar


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
        from astropy.cosmology import realizations

        try:
            value = getattr(realizations, name)
        except AttributeError:
            msg = f"Unknown cosmology {name!r}. Valid cosmologies:\n{realizations.available}"
            raise ValueError(msg) from None

        if not isinstance(value, Cosmology):
            raise TypeError(f"cannot find a Cosmology realization called {name}.")

        return value
