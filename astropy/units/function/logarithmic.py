# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numbers

import numpy as np

from astropy.units import (
    CompositeUnit,
    Unit,
    UnitConversionError,
    UnitsError,
    UnitTypeError,
    dimensionless_unscaled,
)
from astropy.utils import lazyproperty
from astropy.utils.compat.numpycompat import NUMPY_LT_2_0

from .core import FunctionQuantity, FunctionUnitBase

__all__ = [
    "Decibel",
    "DecibelUnit",
    "Dex",
    "DexUnit",
    "LogQuantity",
    "LogUnit",
    "MagUnit",
    "Magnitude",
]


class LogUnit(FunctionUnitBase):
    """Logarithmic unit containing a physical one.

    Usually, logarithmic units are instantiated via specific subclasses
    such `~astropy.units.MagUnit`, `~astropy.units.DecibelUnit`, and
    `~astropy.units.DexUnit`.

    Parameters
    ----------
    physical_unit : `~astropy.units.Unit` or `string`
        Unit that is encapsulated within the logarithmic function unit.
        If not given, dimensionless.

    function_unit :  `~astropy.units.Unit` or `string`
        By default, the same as the logarithmic unit set by the subclass.

    """

    # the four essential overrides of FunctionUnitBase
    @lazyproperty
    def _default_function_unit(self):
        from .units import dex

        return dex

    @property
    def _quantity_class(self):
        return LogQuantity

    def from_physical(self, x):
        """Transformation from value in physical to value in logarithmic units.
        Used in equivalency.
        """
        # Local import to avoid circular dependency.
        from .units import dex

        return dex.to(self._function_unit, np.log10(x))

    def to_physical(self, x):
        """Transformation from value in logarithmic to value in physical units.
        Used in equivalency.
        """
        from .units import dex

        return 10 ** self._function_unit.to(dex, x)

    # ^^^^ the four essential overrides of FunctionUnitBase

    # add addition and subtraction, which imply multiplication/division of
    # the underlying physical units
    def _add_and_adjust_physical_unit(self, other, sign_self, sign_other):
        """Add/subtract LogUnit to/from another unit, and adjust physical unit.

        self and other are multiplied by sign_self and sign_other, resp.

        We wish to do:   ±lu_1 + ±lu_2  -> lu_f          (lu=logarithmic unit)
                  and     pu_1^(±1) * pu_2^(±1) -> pu_f  (pu=physical unit)

        Raises
        ------
        UnitsError
            If function units are not equivalent.
        """
        # First, insist on compatible logarithmic type. Here, plain u.mag,
        # u.dex, and u.dB are OK, i.e., other does not have to be LogUnit
        # (this will indirectly test whether other is a unit at all).
        try:
            getattr(other, "function_unit", other)._to(self._function_unit)
        except AttributeError:
            # if other is not a unit (i.e., does not have _to).
            return NotImplemented
        except UnitsError:
            raise UnitsError(
                "Can only add/subtract logarithmic units of compatible type."
            )

        other_physical_unit = getattr(other, "physical_unit", dimensionless_unscaled)
        physical_unit = CompositeUnit(
            1, [self._physical_unit, other_physical_unit], [sign_self, sign_other]
        )

        return self._copy(physical_unit)

    def __neg__(self):
        return self._copy(self.physical_unit ** (-1))

    def __add__(self, other):
        # Only know how to add to a logarithmic unit with compatible type,
        # be it a plain one (u.mag, etc.,) or another LogUnit
        return self._add_and_adjust_physical_unit(other, +1, +1)

    def __radd__(self, other):
        return self._add_and_adjust_physical_unit(other, +1, +1)

    def __sub__(self, other):
        return self._add_and_adjust_physical_unit(other, +1, -1)

    def __rsub__(self, other):
        # here, in normal usage other cannot be LogUnit; only equivalent one
        # would be u.mag,u.dB,u.dex.  But might as well use common routine.
        return self._add_and_adjust_physical_unit(other, -1, +1)


class MagUnit(LogUnit):
    """Logarithmic physical units expressed in magnitudes.

    Parameters
    ----------
    physical_unit : `~astropy.units.Unit` or `string`
        Unit that is encapsulated within the magnitude function unit.
        If not given, dimensionless.

    function_unit :  `~astropy.units.Unit` or `string`
        By default, this is ``mag``, but this allows one to use an equivalent
        unit such as ``2 mag``.
    """

    @lazyproperty
    def _default_function_unit(self):
        from .units import mag

        return mag

    @property
    def _quantity_class(self):
        return Magnitude


class DexUnit(LogUnit):
    """Logarithmic physical units expressed in magnitudes.

    Parameters
    ----------
    physical_unit : `~astropy.units.Unit` or `string`
        Unit that is encapsulated within the magnitude function unit.
        If not given, dimensionless.

    function_unit :  `~astropy.units.Unit` or `string`
        By default, this is ``dex``, but this allows one to use an equivalent
        unit such as ``0.5 dex``.
    """

    @lazyproperty
    def _default_function_unit(self):
        from .units import dex

        return dex

    @property
    def _quantity_class(self):
        return Dex

    def to_string(self, format="generic"):
        if format == "cds":
            if self.physical_unit == dimensionless_unscaled:
                return "[-]"  # by default, would get "[---]".
            else:
                return f"[{self.physical_unit.to_string(format=format)}]"
        else:
            return super().to_string()


class DecibelUnit(LogUnit):
    """Logarithmic physical units expressed in dB.

    Parameters
    ----------
    physical_unit : `~astropy.units.Unit` or `string`
        Unit that is encapsulated within the decibel function unit.
        If not given, dimensionless.

    function_unit :  `~astropy.units.Unit` or `string`
        By default, this is ``dB``, but this allows one to use an equivalent
        unit such as ``2 dB``.
    """

    @lazyproperty
    def _default_function_unit(self):
        from .units import dB

        return dB

    @property
    def _quantity_class(self):
        return Decibel


class LogQuantity(FunctionQuantity):
    """A representation of a (scaled) logarithm of a number with a unit.

    Parameters
    ----------
    value : number, `~astropy.units.Quantity`, `~astropy.units.LogQuantity`, or sequence of quantity-like.
        The numerical value of the logarithmic quantity. If a number or
        a `~astropy.units.Quantity` with a logarithmic unit, it will be
        converted to ``unit`` and the physical unit will be inferred from
        ``unit``.  If a `~astropy.units.Quantity` with just a physical unit,
        it will converted to the logarithmic unit, after, if necessary,
        converting it to the physical unit inferred from ``unit``.

    unit : str, `~astropy.units.UnitBase`, or `~astropy.units.FunctionUnitBase`, optional
        For an `~astropy.units.FunctionUnitBase` instance, the
        physical unit will be taken from it; for other input, it will be
        inferred from ``value``. By default, ``unit`` is set by the subclass.

    dtype : `~numpy.dtype`, optional
        The ``dtype`` of the resulting Numpy array or scalar that will
        hold the value.  If not provided, is is determined automatically
        from the input value.

    copy : bool, optional
        If `True` (default), then the value is copied.  Otherwise, a copy will
        only be made if ``__array__`` returns a copy, if value is a nested
        sequence, or if a copy is needed to satisfy an explicitly given
        ``dtype``.  (The `False` option is intended mostly for internal use,
        to speed up initialization where a copy is known to have been made.
        Use with care.)

    Examples
    --------
    Typically, use is made of an `~astropy.units.FunctionQuantity`
    subclass, as in::

        >>> import astropy.units as u
        >>> u.Magnitude(-2.5)
        <Magnitude -2.5 mag>
        >>> u.Magnitude(10.*u.count/u.second)
        <Magnitude -2.5 mag(ct / s)>
        >>> u.Decibel(1.*u.W, u.DecibelUnit(u.mW))  # doctest: +FLOAT_CMP
        <Decibel 30. dB(mW)>

    """

    # only override of FunctionQuantity
    _unit_class = LogUnit

    # additions that work just for logarithmic units
    def __add__(self, other):
        # Add function units, thus multiplying physical units. If no unit is
        # given, assume dimensionless_unscaled; this will give the appropriate
        # exception in LogUnit.__add__.
        new_unit = self.unit + getattr(other, "unit", dimensionless_unscaled)
        # Add actual logarithmic values, rescaling, e.g., dB -> dex.
        result = self._function_view + getattr(other, "_function_view", other)
        return self._new_view(result, new_unit)

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        new_unit = self.unit + getattr(other, "unit", dimensionless_unscaled)
        # Do calculation in-place using _function_view of array.
        function_view = self._function_view
        function_view += getattr(other, "_function_view", other)
        self._set_unit(new_unit)
        return self

    def __sub__(self, other):
        # Subtract function units, thus dividing physical units.
        new_unit = self.unit - getattr(other, "unit", dimensionless_unscaled)
        # Subtract actual logarithmic values, rescaling, e.g., dB -> dex.
        result = self._function_view - getattr(other, "_function_view", other)
        return self._new_view(result, new_unit)

    def __rsub__(self, other):
        new_unit = self.unit.__rsub__(getattr(other, "unit", dimensionless_unscaled))
        result = self._function_view.__rsub__(getattr(other, "_function_view", other))
        # Ensure the result is in right function unit scale
        # (with rsub, this does not have to be one's own).
        result = result.to(new_unit.function_unit)
        return self._new_view(result, new_unit)

    def __isub__(self, other):
        new_unit = self.unit - getattr(other, "unit", dimensionless_unscaled)
        # Do calculation in-place using _function_view of array.
        function_view = self._function_view
        function_view -= getattr(other, "_function_view", other)
        self._set_unit(new_unit)
        return self

    def __mul__(self, other):
        # Multiply by a float or a dimensionless quantity
        if isinstance(other, numbers.Number):
            # Multiplying a log means putting the factor into the exponent
            # of the unit
            new_physical_unit = self.unit.physical_unit**other
            result = self.view(np.ndarray) * other
            return self._new_view(result, self.unit._copy(new_physical_unit))
        else:
            return super().__mul__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
        if isinstance(other, numbers.Number):
            new_physical_unit = self.unit.physical_unit**other
            function_view = self._function_view
            function_view *= other
            self._set_unit(self.unit._copy(new_physical_unit))
            return self
        else:
            return super().__imul__(other)

    def __truediv__(self, other):
        # Divide by a float or a dimensionless quantity
        if isinstance(other, numbers.Number):
            # Dividing a log means putting the denominator into the exponent
            # of the unit
            new_physical_unit = self.unit.physical_unit ** (1 / other)
            result = self.view(np.ndarray) / other
            return self._new_view(result, self.unit._copy(new_physical_unit))
        else:
            return super().__truediv__(other)

    def __itruediv__(self, other):
        if isinstance(other, numbers.Number):
            new_physical_unit = self.unit.physical_unit ** (1 / other)
            function_view = self._function_view
            function_view /= other
            self._set_unit(self.unit._copy(new_physical_unit))
            return self
        else:
            return super().__itruediv__(other)

    def __pow__(self, other):
        # We check if this power is OK by applying it first to the unit.
        try:
            other = float(other)
        except TypeError:
            return NotImplemented
        new_unit = self.unit**other
        new_value = self.view(np.ndarray) ** other
        return self._new_view(new_value, new_unit)

    def __ilshift__(self, other):
        try:
            other = Unit(other)
        except UnitTypeError:
            return NotImplemented

        if not isinstance(other, self._unit_class):
            return NotImplemented

        try:
            factor = self.unit.physical_unit._to(other.physical_unit)
        except UnitConversionError:
            # Maybe via equivalencies?  Now we do make a temporary copy.
            try:
                value = self._to_value(other)
            except UnitConversionError:
                return NotImplemented

            self.view(np.ndarray)[...] = value
        else:
            self.view(np.ndarray)[...] += self.unit.from_physical(factor)

        self._set_unit(other)
        return self

    # Methods that do not work for function units generally but are OK for
    # logarithmic units as they imply differences and independence of
    # physical unit.
    def var(self, axis=None, dtype=None, out=None, ddof=0):
        unit = self.unit.function_unit**2
        return self._wrap_function(np.var, axis, dtype, out=out, ddof=ddof, unit=unit)

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        unit = self.unit._copy(dimensionless_unscaled)
        return self._wrap_function(np.std, axis, dtype, out=out, ddof=ddof, unit=unit)

    if NUMPY_LT_2_0:

        def ptp(self, axis=None, out=None):
            unit = self.unit._copy(dimensionless_unscaled)
            return self._wrap_function(np.ptp, axis, out=out, unit=unit)

    else:

        def __array_function__(self, function, types, args, kwargs):
            # TODO: generalize this to all supported functions!
            if function is np.ptp:
                unit = self.unit._copy(dimensionless_unscaled)
                return self._wrap_function(np.ptp, *args[1:], unit=unit, **kwargs)
            else:
                return super().__array_function__(function, types, args, kwargs)

    def diff(self, n=1, axis=-1):
        unit = self.unit._copy(dimensionless_unscaled)
        return self._wrap_function(np.diff, n, axis, unit=unit)

    def ediff1d(self, to_end=None, to_begin=None):
        unit = self.unit._copy(dimensionless_unscaled)
        return self._wrap_function(np.ediff1d, to_end, to_begin, unit=unit)

    _supported_functions = FunctionQuantity._supported_functions | {
        getattr(np, function) for function in ("var", "std", "ptp", "diff", "ediff1d")
    }


class Dex(LogQuantity):
    _unit_class = DexUnit


class Decibel(LogQuantity):
    _unit_class = DecibelUnit


class Magnitude(LogQuantity):
    _unit_class = MagUnit
