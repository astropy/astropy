# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import abc
from numbers import Number
from collections.abc import ItemsView
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

    def __init__(self, target, qcls=Quantity, bound=True, **kw):
        self.target = target
        self.qcls = qcls
        self.bound = bound

    @classmethod
    def _get_action(cls, kls=None):
        # needed for `UnitSpec.__init_subclass__`, where `kls` is is passed by
        # the wrap of `_get_action` by `classproperty`
        return _USPEC_ACTION_REGISTRY[kls or cls]

    @classproperty
    def action(cls):
        """Get action associated with unit specification."""
        return cls._get_action()

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
        return f"{self.__class__.__qualname__}({self.target}, qcls={self.qcls.__name__})"


class CompoundUnitSpec(UnitSpecBase):

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
                return out

        raise TypeError("All unitspecs failed")  # TODO! better error message

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


class UnitSpec(UnitSpecBase):
    """Unit Specification.

    Parameters
    ----------
    target : `UnitSpecBase` or None or Any
    action : str or None, optional
        The action to take on the target.
    qcls : `~astropy.units.Quantity` or subclass type
        The Quantity class. This can be used to enforce, e.g.
        `~astropy.coordinates.Angle` quantities.
    """

    def __init_subclass__(cls, action, **kwargs):
        # check registry
        if action in _ACTION_USPEC_REGISTRY:
            raise KeyError(f"`{action}` already in registry.")

        # register subclass and action
        _ACTION_USPEC_REGISTRY[action] = cls
        _USPEC_ACTION_REGISTRY[cls] = action

        # reimplement UnitSpecBase action
        cls.action = classproperty(UnitSpecBase._get_action)

    def __new__(cls, target=None, action=None, qcls=Quantity, **kw):
        # First check if target is a UnitSpec. This class doesn't convert, so
        # those pass thru directly. Note: UnitSpec subclasses can do whatever.
        if isinstance(target, UnitSpecBase):
            return target

        # second check it's a valid action...
        elif action not in _ACTION_USPEC_REGISTRY:
            raise KeyError(f"action {action!r} must be one of {_ACTION_USPEC_REGISTRY.keys()}")

        # if so, create the UnitSpec for that action.
        self = super().__new__(cls)
        if cls is UnitSpec:
            # UnitSpec wraps its subclasses so only UnitSpec need ever be imported.
            # Need to determine which subclass to wrap
            self._uspec = _ACTION_USPEC_REGISTRY[action](target, qcls=qcls)
            self.__doc__ = self._uspec.__doc__  # update docstring
            # self.__call__.__doc__ = self._uspec.__call__.__doc__

        return self

    def __call__(self, value, **kwargs):
        """Calls the wrapped UnitSpec."""
        # pass through to wrapped UnitSpec
        return self._uspec(value, **kwargs)

    @property
    def action(self):
        """Unit Specification."""
        return self._uspec.action

    @action.setter
    def action(self, value):
        if value != self._uspec.action:
            self._uspec = self._ACTION_USPEC_REGISTRY[value](self.target, qcls=self.qcls)

    # =====================

    def __repr__(self):
        if self.__class__ is UnitSpec:
            return f"[UnitSpec]{self._uspec}"

        # can / should be overridden by subclasses. This is just a default.
        return (f"{self.__class__.__qualname__}({self.target!r}, qcls={self.qcls.__qualname__}, "
                f"action={self.action!r}, bound={self.bound})")


# -----------------------------------------------------------------------------


class NullUnitSpec(UnitSpec, action=None):
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
    
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, value):
        """Returns 'value' unchanged."""
        return value

    def __repr__(self):
        return f"{self.__class__.__qualname__}()"


# -----------------------------------------------------------------------------


class ValidatePhysicalType(UnitSpec, action="validate"):
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
                                   action='validate', strict=False)
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


# -----------------------------------------------------------------------------


class ConvertToUnit(UnitSpec, action="to unit"):
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


class ConvertToValue(ConvertToUnit, action="to value"):
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


class AssignUnits(ConvertToUnit, action="from value"):
    """Assign input target units.
    
    Equivalent to UnitSpec 'to unit', with strictness set to False.

    Parameters
    ----------
    unit : unit-like
        Not Physical Type
    qcls : type
    """

    def __init__(self, unit, qcls=Quantity, **kw):
        kw.pop("bound", None)  # make sure not present
        super().__init__(unit, qcls=qcls, bound=True, **kw)

    def __call__(self, value, **kw):
        """Assign input target units.
    
        Parameters
        ----------
        value
    
        Returns
        -------
        `~astropy.units.Quantity`
        """
        kw.pop("bound", None)  # make sure not present
        return self.qcls(value, unit=self.target, copy=False)
