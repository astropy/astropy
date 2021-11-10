# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import abc
from numbers import Number
from collections.abc import ItemsView
from inspect import _empty
from types import MappingProxyType

import numpy as np

from astropy.utils.decorators import classproperty

from .core import Unit, UnitBase, UnitConversionError, UnitsError
from .physical import _physical_unit_mapping, _unit_physical_mapping, get_physical_type
from .quantity import Quantity

__all__ = ["UnitSpec", "CompoundUnitSpec", "NullUnitSpec", "UnitSpecBase"]

# registry of UnitSpec action options
_ACTION_USPEC_REGISTRY = {}
_USPEC_ACTION_REGISTRY = {}


class UnitSpecBase(metaclass=abc.ABCMeta):
    """Abstract base class for Unit Specification.

    Use this in ``isinstance`` checks.

    Parameters
    ----------
    target : Any
    qcls : `~astropy.units.Quantity` or subclass type
        The Quantity class. This can be used to enforce, e.g.
        `~astropy.coordinates.Angle` quantities.
    bound : bool, optional
        Whether the target type must be exactly 'qcls' (False) or any subtype
        (inclusive) of 'qcls' (True, default).
    **kw
        Not used.
    """

    def __init_subclass__(cls, action, **kwargs):
        # check registry
        if action in _ACTION_USPEC_REGISTRY:
            raise KeyError(f"`{action}` already in registry.")

        if action is not _empty:  # optionally register
            # register subclass and action
            _ACTION_USPEC_REGISTRY[action] = cls
            _USPEC_ACTION_REGISTRY[cls] = action

    @classproperty
    def action(cls):
        """Get action associated with unit specification."""
        return _USPEC_ACTION_REGISTRY[cls]

    def __init__(self, target, qcls=Quantity, bound=True, **kw):
        self.target = target
        self.qcls = qcls
        self.bound = bound

    @abc.abstractmethod
    def __call__(self, value, **kwargs):
        """The call method must be implemented by subclasses.

        Subclasses should NOT call this with ``super``.
        This raises a `NotImplementedError`.
        """
        raise NotImplementedError

    def _check_bound(self, value, bound=None):
        bound = self.bound if bound is None else bound

        if not bound and value.__class__ is not self.qcls:  # the strict check
            raise TypeError(f"{value!r} is not type {self.qcls!r}.")
        elif not isinstance(value, self.qcls):  # bound=True but wrong type
            raise TypeError(f"{value!r} is not qcls {self.qcls!r} (or subclass).")

        return True

    # =====================

    def __add__(self, other):
        # check if other is UnitSpecBase
        if not isinstance(other, UnitSpecBase):
            return NotImplemented

        return CompoundUnitSpec(self, other)

    def __repr__(self):
        return (f"{self.__class__.__qualname__}({self.target},"
                f" qcls={self.qcls.__name__}, bound={self.bound})")


class CompoundUnitSpec(UnitSpecBase, action=_empty):

    def __init__(self, *specs):
        # TODO! flatten specs to allow for CompoundUnitSpec(CompoundUnitSpec())
        # also check if there are conflicting specs
        self._specs = specs

    def __call__(self, value, **kwargs):
        # TODO! better iteration over specs
        for spec in self._specs:
            try:
                out = spec(value, **kwargs)
            except (UnitConversionError, TypeError):
                pass
            else:
                break
        else:
            # TODO! better error message
            raise TypeError("All unitspecs failed")  # TODO! better error message

        return out

    # =====================

    def __iadd__(self, other):
        # check if other is UnitSpecBase
        if not isinstance(other, UnitSpecBase):
            raise TypeError(f"{other!r} must be {UnitSpecBase!r}")

        # allow for other to be CompoundUnitSpec
        if isinstance(other, CompoundUnitSpec):
            other = other._specs
        else:
            other = (other, )

        self._specs = self._specs + other
        return self

    def __repr__(self):
        # TODO! nicely split long lines
        return f"{self.__class__.__qualname__}{self._specs!r}"

# =============================================================================


# TODO! this as a class? Make it the base of the non-compound chain?
# Subclasses return themselves, UnitSpec returns a subclass.
def UnitSpec(target, action=None, qcls=Quantity, bound=True, **kwargs):
    # First check if target is a UnitSpec. This class doesn't convert, so
    # those pass thru directly. Note: UnitSpec subclasses can do whatever.
    if isinstance(target, UnitSpecBase):
        return target

    if action not in _ACTION_USPEC_REGISTRY:
        raise KeyError(f"action {action!r} must be one of {_ACTION_USPEC_REGISTRY.keys()}")

    cls = _ACTION_USPEC_REGISTRY[action]
    uspec = cls(target, qcls=qcls, bound=bound, **kwargs)
    return uspec

# class UnitSpec(UnitSpecBase, action=_empty):
#     """Unit Specification.
#
#     Should not be subclassed.
#
#     Parameters
#     ----------
#     target : `UnitSpecBase` or None or Any
#     action : str or None, optional
#         The action to take on the target.
#     qcls : `~astropy.units.Quantity` or subclass type
#         The Quantity class. This can be used to enforce, e.g.
#         `~astropy.coordinates.Angle` quantities.
#     """
#
#     def __new__(cls, target=None, action=None, qcls=Quantity, bound=True, **kw):
#         if cls is not UnitSpec:
#             raise TypeError("UnitSpec cannot be subclassed.")
#
#         # First check if target is a UnitSpec. This class doesn't convert, so
#         # those pass thru directly. Note: UnitSpec subclasses can do whatever.
#         if isinstance(target, UnitSpecBase):
#             return target
#
#         return super().__new__(cls)
#
#     def __init__(self, target, action=None, qcls=Quantity, bound=True, **kw):
#         # second check it's a valid action...
#         if action not in _ACTION_USPEC_REGISTRY:
#             raise KeyError(f"action {action!r} must be one of {_ACTION_USPEC_REGISTRY.keys()}")
#
#         # UnitSpec wraps its subclasses so only UnitSpec need ever be imported.
#         # Need to determine which subclass to wrap
#         uspec = _ACTION_USPEC_REGISTRY[action](target, qcls=qcls, bound=bound, **kw)
#         object.__setattr__(self, "_uspec", uspec)
#         self.__doc__ = self._uspec.__doc__  # update docstring
#         # self.__call__.__doc__ = self._uspec.__call__.__doc__
#
#     def __call__(self, value, **kwargs):
#         """Calls the wrapped UnitSpec."""
#         # pass through to wrapped UnitSpec
#         return self._uspec(value, **kwargs)
#
#     @property
#     def action(self):
#         """Unit Specification."""
#         return self._uspec.action
#
#     def __getattr__(self, attr):
#         return getattr(self._uspec, attr)
#
#     # def __setattr__(self, attr, val):
#     #     setattr(self._uspec, attr, val)
#
#     def __dir__(self):
#         return sorted(set(object.__dir__(self)).union(set(dir(self._uspec))))
#
#     # =====================
#
#     def __repr__(self):
#         return "[UnitSpec]" + repr(self._uspec).replace("UnitSpec", "")


# -----------------------------------------------------------------------------


class NullUnitSpec(UnitSpecBase, action=None):
    """Null Unit Specification. Returns value on call.

    Examples
    --------
    >>> from astropy.units import NullUnitSpec
    >>> nulluspec = NullUnitSpec()
    >>> nulluspec
    NullUnitSpec()

    >>> print(nulluspec.action)
    None

    >>> nulluspec(1)
    1
    >>> nulluspec([1, 2, 3])
    [1, 2, 3]

    Alternatively to construct a `~astropy.units.unitspec.NullUnitSpec`
    with `~astropy.units.unitspec.UnitSpec`:

    >>> from astropy.units import UnitSpec
    >>> UnitSpec(None)
    [UnitSpec]NullUnitSpec()
    """

    def __init__(self, target=None, *args, **kwargs):
        if target is not None:
            raise ValidateError("NullUnitSpec can only have target=None.")
        pass

    def __call__(self, value):
        """Returns 'value' unchanged."""
        return value

    def __repr__(self):
        return f"{self.__class__.__qualname__}()"


# -----------------------------------------------------------------------------


class UnitSpecValidatePhysicalType(UnitSpecBase, action="validate"):
    """UnitSpec to validate physical type.

    Parameters
    ----------
    physical_type : `~astropy.units.PhysicalType`-like
        e.g. the string 'length'
        This is the ``target`` in `UnitSpecBase`.
    qcls : type
        The :class:`~astropy.units.Quantity` class.
        Can be used to check that an input is not only the correct physical
        type, but also the correct Quantity type.
    bound : bool, optional
        Whether the target type must be exactly 'qcls' (False) or any subtype
        (inclusive) of 'qcls' (True, default).

    strict_dimensionless : bool, optional
        Whether to be strict about stuff like numerics (`~number.Number` or
        `~numpy.ndarray`) being equivalent to dimensionless quantities.

    Examples
    --------
    >>> uspec = UnitSpec("length", action="validate")
    >>> uspec
    [UnitSpec]ValidatePhysicalType(PhysicalType('length'), qcls=Quantity,
                                   action='validate', bound=True)
    """

    def __init__(self, physical_type, qcls=Quantity, bound=True, strict_dimensionless=False, **kw):
        # make sure it's actually a physical type. this will catch any errors.
        ptype = get_physical_type(physical_type)
        super().__init__(ptype, qcls=qcls, bound=bound, **kw)

        self.strict_dimensionless = strict_dimensionless

    def __call__(self, value, bound=None, strict_dimensionless=None, **kw):
        """

        Parameters
        ----------
        value : Any
        strict_dimensionless : bool or None, optional
            None uses default value set at initialization.
        **kw
            Not used.
        """
        # what if 'value' was just a number? if not 'strict', numbers will be
        # treated as dimensionless quantities.
        strict = (self.strict_dimensionless if strict_dimensionless is None
                  else strict_dimensionless)
        bound = self.bound if bound is None else bound
        if (not strict
            and (self.target == "dimensionless")
            and (isinstance(value, (Number, np.generic))
                 or (isinstance(value, np.ndarray)
                     and np.issubdtype(value.dtype, np.number)))):
            pass  # He is the `~astropy.units.one`.

        # checks that 'value' is the the correct class type
        elif not self._check_bound(value, bound=bound):
            pass # error messages handled in _check_bound

        # check units is the correct physical type, including any enabled
        # equivalencies. TODO! a cleaner check
        if not value.unit.is_equivalent(self.target._unit):
            raise UnitConversionError(f"`{value.unit.physical_type}` is not equivalent to type '{self.target}'")
        else:
            pass  # passes all checks!

        return value

    def __repr__(self):
        return (f"{self.__class__.__qualname__}({self.target},"
                f" qcls={self.qcls.__name__}, bound={self.bound},"
                f" strict_dimensionless={self.strict_dimensionless})")


# -----------------------------------------------------------------------------


class UnitSpecConvertToUnit(UnitSpecBase, action="to unit"):
    """Convert input to target unit.

    Parameters
    ----------
    unit : unit-like
        Not Physical Type
    qcls : type
    bound : bool
    """

    def __init__(self, unit, qcls=Quantity, bound=True, **kw):
        unit = Unit(unit)  # make sure it's a Unit. catches any errors.
        super().__init__(unit, qcls=qcls, bound=bound, **kw)

    def __call__(self, value, bound=None, **kw):
        """Convert input to target units.

        Parameters
        ----------
        value
        bound : bool

        Returns
        -------
        `~astropy.units.Quantity`
        """
        self._check_bound(value, bound=bound)
        # convert value to desired quantity
        return self.qcls(value, unit=self.target, copy=False)


class UnitSpecConvertToValue(UnitSpecConvertToUnit, action="to value"):
    """Convert input to value in target units.

    Parameters
    ----------
    unit : unit-like
        Not Physical Type
    qcls : type
    bound : bool
    """

    def __call__(self, value, bound=None, **kw):
        """Convert input to value in target units.

        Parameters
        ----------
        value
        bound : bool

        Returns
        -------
        `~astropy.units.Quantity`
        """
        return super().__call__(value, bound=bound, **kw).to_value()


class UnitSpecAssignUnits(UnitSpecConvertToUnit, action="from value"):
    """Assign input target units.

    Equivalent to UnitSpec 'to unit', with strictness set to False.

    Parameters
    ----------
    unit : unit-like
        Not Physical Type
    qcls : type
    """

    def __init__(self, unit, qcls=Quantity, bound=True, **kw):
        super().__init__(unit, qcls=qcls, bound=bound, **kw)

    def __call__(self, value, bound=None, **kw):
        """Assign input target units.

        Parameters
        ----------
        value

        Returns
        -------
        `~astropy.units.Quantity`
        """
        # only apply type bound check if it's like qcls
        # numbers and 'raw' input are *assigned* units, not checked.
        if hasattr(value, "unit"):  # TODO? is this a sufficient check, or too general?
            self._check_bound(value, bound=bound)

        kw.setdefault("copy", False)
        return self.qcls(value, unit=self.target, **kw)
