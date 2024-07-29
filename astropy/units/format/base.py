# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

from typing import TYPE_CHECKING

from . import utils

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import ClassVar, Literal

    import numpy as np

    from astropy.units import NamedUnit, UnitBase
    from astropy.units.typing import Real


class Base:
    """
    The abstract base class of all unit formats.
    """

    registry: ClassVar[dict[str, type[Base]]] = {}
    _space: ClassVar[str] = " "
    _scale_unit_separator: ClassVar[str] = " "
    _times: ClassVar[str] = "*"
    name: ClassVar[str]  # Set by __init_subclass__ by the latest

    def __new__(cls, *args, **kwargs):
        # This __new__ is to make it clear that there is no reason to
        # instantiate a Formatter--if you try to you'll just get back the
        # class
        return cls

    def __init_subclass__(cls, **kwargs):
        # Keep a registry of all formats.  Key by the class name unless a name
        # is explicitly set (i.e., one *not* inherited from a superclass).
        if "name" not in cls.__dict__:
            cls.name = cls.__name__.lower()

        Base.registry[cls.name] = cls
        super().__init_subclass__(**kwargs)

    @classmethod
    def format_exponential_notation(
        cls, val: float | np.number, format_spec: str = ".8g"
    ) -> str:
        """
        Formats a value in exponential notation.

        Parameters
        ----------
        val : number
            The value to be formatted

        format_spec : str, optional
            Format used to split up mantissa and exponent

        Returns
        -------
        str
            The value in exponential notation in a this class's format.
        """
        x = format(val, format_spec).split("e")
        if len(x) != 2:
            return cls._format_mantissa(x[0])  # no exponent
        ex = x[1].lstrip("0+")
        if not ex:
            return cls._format_mantissa(x[0])  # exponent was zero
        if ex.startswith("-"):
            ex = "-" + ex[1:].lstrip("0")
        ex = f"10{cls._format_superscript(ex)}"
        m = cls._format_mantissa("" if x[0].rstrip("0") == "1." else x[0])
        return f"{m}{cls._times}{ex}" if m else ex

    @classmethod
    def _format_mantissa(cls, m: str) -> str:
        return m

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return f"({number})" if "/" in number or "." in number else number

    @classmethod
    def _format_unit_power(cls, unit: NamedUnit, power: Real = 1) -> str:
        """Format the unit for this format class raised to the given power.

        This is overridden in Latex where the name of the unit can depend on the power
        (e.g., for degrees).
        """
        name = unit.get_format_name(cls.name)
        if power != 1:
            name += cls._format_superscript(utils.format_power(power))
        return name

    @classmethod
    def _format_unit_list(cls, units: Iterable[tuple[NamedUnit, Real]]) -> str:
        return cls._space.join(
            cls._format_unit_power(base_, power) for base_, power in units
        )

    @classmethod
    def _format_fraction(
        cls,
        scale: str,
        numerator: str,
        denominator: str,
        *,
        fraction: Literal[True, "inline"] = "inline",
    ) -> str:
        if not (fraction is True or fraction == "inline"):
            raise ValueError(
                "format {cls.name!r} only supports inline fractions,"
                f"not fraction={fraction!r}."
            )

        if cls._space in denominator:
            denominator = f"({denominator})"
        if scale and numerator == "1":
            return f"{scale}/ {denominator}"
        return f"{scale}{numerator} / {denominator}"

    @classmethod
    def to_string(
        cls, unit: UnitBase, *, fraction: bool | Literal["inline", "multiline"] = True
    ) -> str:
        """Convert a unit to its string representation.

        Implementation for `~astropy.units.UnitBase.to_string`.

        Parameters
        ----------
        unit : |Unit|
            The unit to convert.
        fraction : {False|True|'inline'|'multiline'}, optional
            Options are as follows:

            - `False` : display unit bases with negative powers as they are
              (e.g., ``km s-1``);
            - 'inline' or `True` : use a single-line fraction (e.g., ``km / s``);
            - 'multiline' : use a multiline fraction (available for the
              ``latex``, ``console`` and ``unicode`` formats only; e.g.,
              ``$\\mathrm{\\frac{km}{s}}$``).

        Raises
        ------
        ValueError
            If ``fraction`` is not recognized.
        """
        # First the scale.  Normally unity, in which case we omit
        # it, but non-unity scale can happen, e.g., in decompositions
        # like u.Ry.decompose(), which gives "2.17987e-18 kg m2 / s2".
        if unit.scale == 1:
            s = ""
        else:
            s = cls.format_exponential_notation(unit.scale)

        # Now the unit baes, taking care that dimensionless does not have any
        # (but can have a scale; e.g., u.percent.decompose() gives "0.01").
        if len(unit.bases):
            if s:
                s += cls._scale_unit_separator
            if fraction:
                numerator, denominator = utils.get_grouped_by_powers(
                    unit.bases, unit.powers
                )
            else:
                numerator = list(unit.components())
                denominator = []
            if len(denominator):
                if len(numerator):
                    numerator = cls._format_unit_list(numerator)
                else:
                    numerator = "1"
                denominator = cls._format_unit_list(denominator)
                s = cls._format_fraction(s, numerator, denominator, fraction=fraction)
            else:
                s += cls._format_unit_list(numerator)

        return s

    @classmethod
    def parse(cls, s: str) -> UnitBase:
        """
        Convert a string to a unit object.
        """
        raise NotImplementedError(f"Can not parse with {cls.__name__} format")
