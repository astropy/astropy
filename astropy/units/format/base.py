# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from astropy.units.core import CompositeUnit, NamedUnit, Unit, get_current_unit_registry
from astropy.units.errors import UnitsWarning
from astropy.units.utils import maybe_simple_fraction
from astropy.utils.misc import did_you_mean

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import ClassVar, Literal

    import numpy as np

    from astropy.extern.ply.lex import LexToken
    from astropy.units import UnitBase
    from astropy.units.typing import UnitPower, UnitScale


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
        cls, val: UnitScale | np.number, format_spec: str = ".8g"
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
    def _format_unit_power(cls, unit: NamedUnit, power: UnitPower = 1) -> str:
        """Format the unit for this format class raised to the given power.

        This is overridden in Latex where the name of the unit can depend on the power
        (e.g., for degrees).
        """
        name = unit._get_format_name(cls.name)
        return name if power == 1 else name + cls._format_power(power)

    @classmethod
    def _format_power(cls, power: UnitPower) -> str:
        # If the denominator of `power` is a power of 2 then `power` is stored
        # as a `float` (see `units.utils.sanitize_power()`), but we still want
        # to display it as a fraction.
        return cls._format_superscript(
            str(maybe_simple_fraction(power) if isinstance(power, float) else power)
        )

    @classmethod
    def _format_unit_list(cls, units: Iterable[tuple[NamedUnit, UnitPower]]) -> str:
        return cls._space.join(
            cls._format_unit_power(base_, power) for base_, power in units
        )

    @classmethod
    def _format_inline_fraction(
        cls, scale: str, numerator: str, denominator: str
    ) -> str:
        if cls._space in denominator:
            denominator = f"({denominator})"
        if scale and numerator == "1":
            return f"{scale}/ {denominator}"
        return f"{scale}{numerator} / {denominator}"

    @classmethod
    def _format_multiline_fraction(
        cls, scale: str, numerator: str, denominator: str
    ) -> str:
        # By default, we just warn that we do not have a multiline format.
        warnings.warn(
            f"{cls.name!r} format does not support multiline "
            "fractions; using inline instead.",
            UnitsWarning,
        )
        return cls._format_inline_fraction(scale, numerator, denominator)

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
            - 'multiline' : use a multiline fraction if possible (available for
              the ``latex``, ``console`` and ``unicode`` formats; e.g.,
              ``$\\mathrm{\\frac{km}{s}}$``). If not possible, use 'inline'.

        Raises
        ------
        ValueError
            If ``fraction`` is not recognized.
        """
        # First the scale.  Normally unity, in which case we omit
        # it, but non-unity scale can happen, e.g., in decompositions
        # like u.Ry.decompose(), which gives "2.17987e-18 kg m2 / s2".
        s = "" if unit.scale == 1.0 else cls.format_exponential_notation(unit.scale)

        # dimensionless does not have any bases, but can have a scale;
        # e.g., u.percent.decompose() gives "0.01".
        if not unit.bases:
            return s

        if s:
            s += cls._scale_unit_separator
        # Unit powers are monotonically decreasing
        if not fraction or unit.powers[-1] > 0:
            return s + cls._format_unit_list(zip(unit.bases, unit.powers, strict=True))

        if fraction is True or fraction == "inline":
            formatter = cls._format_inline_fraction
        elif fraction == "multiline":
            formatter = cls._format_multiline_fraction
        else:
            raise ValueError(
                "fraction can only be False, 'inline', or 'multiline', "
                f"not {fraction!r}."
            )

        positive = []
        negative = []
        for base, power in zip(unit.bases, unit.powers, strict=True):
            if power > 0:
                positive.append((base, power))
            else:
                negative.append((base, -power))
        return formatter(
            s, cls._format_unit_list(positive) or "1", cls._format_unit_list(negative)
        )

    @classmethod
    def parse(cls, s: str) -> UnitBase:
        """
        Convert a string to a unit object.
        """
        raise NotImplementedError(f"Can not parse with {cls.__name__} format")


class _ParsingFormatMixin:
    """Provides private methods used in the formats that parse units."""

    _deprecated_units: ClassVar[frozenset[str]] = frozenset()

    @classmethod
    def _do_parse(cls, s: str, debug: bool = False) -> UnitBase:
        try:
            return cls._parser.parse(s, lexer=cls._lexer, debug=debug)
        except ValueError as e:
            if str(e):
                raise
            else:
                raise ValueError(f"Syntax error parsing unit '{s}'")

    @classmethod
    def _get_unit(cls, t: LexToken) -> UnitBase:
        try:
            return cls._validate_unit(t.value)
        except ValueError as e:
            registry = get_current_unit_registry()
            if t.value in registry.aliases:
                return registry.aliases[t.value]

            raise ValueError(f"At col {t.lexpos}, {str(e)}")

    @classmethod
    def _fix_deprecated(cls, x: str) -> list[str]:
        return [x + " (deprecated)" if x in cls._deprecated_units else x]

    @classmethod
    def _did_you_mean_units(cls, unit: str) -> str:
        """
        A wrapper around `astropy.utils.misc.did_you_mean` that deals with
        the display of deprecated units.

        Parameters
        ----------
        unit : str
            The invalid unit string

        Returns
        -------
        msg : str
            A message with alternatives, or the empty string.
        """
        return did_you_mean(unit, cls._units, fix=cls._fix_deprecated)

    @classmethod
    def _validate_unit(cls, unit: str, detailed_exception: bool = True) -> UnitBase:
        try:
            return cls._units[unit]
        except KeyError:
            if detailed_exception:
                raise ValueError(cls._invalid_unit_error_message(unit)) from None
            raise ValueError() from None

    @classmethod
    def _invalid_unit_error_message(cls, unit: str) -> str:
        return (
            f"Unit '{unit}' not supported by the {cls.__name__} standard. "
            + cls._did_you_mean_units(unit)
        )

    @classmethod
    def _decompose_to_known_units(cls, unit: CompositeUnit | NamedUnit) -> UnitBase:
        """
        Partially decomposes a unit so it is only composed of units that
        are "known" to a given format.
        """
        if isinstance(unit, CompositeUnit):
            return CompositeUnit(
                unit.scale,
                [cls._decompose_to_known_units(base) for base in unit.bases],
                unit.powers,
                _error_check=False,
            )
        if isinstance(unit, NamedUnit):
            try:
                return cls._validate_unit(unit._get_format_name(cls.name))
            except ValueError:
                if isinstance(unit, Unit):
                    return cls._decompose_to_known_units(unit._represents)
                raise
        raise TypeError(
            f"unit argument must be a 'NamedUnit' or 'CompositeUnit', not {type(unit)}"
        )
