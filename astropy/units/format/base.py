# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

from typing import TYPE_CHECKING

from astropy.units.utils import maybe_simple_fraction
from astropy.utils import classproperty

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from typing import ClassVar, Literal

    import numpy as np

    from astropy.units import NamedUnit, UnitBase
    from astropy.units.typing import Real, UnitPower


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

    @classproperty(lazy=True)
    def _fraction_formatters(cls) -> dict[bool | str, Callable[[str, str, str], str]]:
        return {
            True: cls._format_inline_fraction,
            "inline": cls._format_inline_fraction,
        }

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
    def _format_unit_list(cls, units: Iterable[tuple[NamedUnit, Real]]) -> str:
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
    def to_string(
        cls, unit: UnitBase, *, fraction: bool | Literal["inline"] = True
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
        # A separate `_to_string()` method allows subclasses (e.g. `FITS`) to implement
        # `to_string()` without needlessly interfering with the `to_string()`
        # implementations of their subclasses.
        return cls._to_string(unit, fraction=fraction)

    @classmethod
    def _to_string(cls, unit: UnitBase, *, fraction: bool | str) -> str:
        # First the scale.  Normally unity, in which case we omit
        # it, but non-unity scale can happen, e.g., in decompositions
        # like u.Ry.decompose(), which gives "2.17987e-18 kg m2 / s2".
        s = "" if unit.scale == 1 else cls.format_exponential_notation(unit.scale)

        # dimensionless does not have any bases, but can have a scale;
        # e.g., u.percent.decompose() gives "0.01".
        if not unit.bases:
            return s

        if s:
            s += cls._scale_unit_separator
        # Unit powers are monotonically decreasing
        if not fraction or unit.powers[-1] > 0:
            return s + cls._format_unit_list(zip(unit.bases, unit.powers, strict=True))

        positive = []
        negative = []
        for base, power in zip(unit.bases, unit.powers, strict=True):
            if power > 0:
                positive.append((base, power))
            else:
                negative.append((base, -power))
        try:
            return cls._fraction_formatters[fraction](
                s,
                cls._format_unit_list(positive) or "1",
                cls._format_unit_list(negative),
            )
        except KeyError:
            # We accept Booleans, but don't advertise them in the error message
            *all_but_last, last = (
                repr(key) for key in cls._fraction_formatters if isinstance(key, str)
            )
            supported_formats = (
                f"{', '.join(all_but_last)} or {last}" if all_but_last else last
            )
            raise ValueError(
                f"{cls.name!r} format only supports {supported_formats} "
                f"fractions, not {fraction=!r}."
            ) from None

    @classmethod
    def parse(cls, s: str) -> UnitBase:
        """
        Convert a string to a unit object.
        """
        raise NotImplementedError(f"Can not parse with {cls.__name__} format")
