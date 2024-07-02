# Licensed under a 3-clause BSD style license - see LICENSE.rst
# ruff: noqa: RUF009

from __future__ import annotations

import inspect
from abc import ABCMeta, abstractmethod
from dataclasses import KW_ONLY, dataclass, replace
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar

import numpy as np

from astropy.io.registry import UnifiedReadWriteMethod
from astropy.utils.decorators import classproperty
from astropy.utils.metadata import MetaData

from ._utils import all_parameters
from .connect import (
    CosmologyFromFormat,
    CosmologyRead,
    CosmologyToFormat,
    CosmologyWrite,
)
from .parameter import Parameter
from .parameter._descriptors import ParametersAttribute

if TYPE_CHECKING:
    from collections.abc import Mapping

    from typing_extensions import Self

    from astropy.cosmology.funcs.comparison import _FormatType

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com), Roban Kramer
# (robanhk@gmail.com), and Nathaniel Starkman (n.starkman@mail.utoronto.ca).

# Many of these adapted from Hogg 1999, astro-ph/9905116
# and Linder 2003, PRL 90, 91301

__all__ = ["Cosmology", "CosmologyError", "FlatCosmologyMixin"]


##############################################################################
# Parameters

# registry of cosmology classes with {key=name : value=class}
_COSMOLOGY_CLASSES = dict()

# typing
# NOTE: private b/c RTD error
_CosmoT = TypeVar("_CosmoT", bound="Cosmology")
_FlatCosmoT = TypeVar("_FlatCosmoT", bound="FlatCosmologyMixin")


# dataclass transformation
def _with_signature(cls: type[Cosmology]) -> type[Cosmology]:
    """Decorator to precompute the class' signature.

    This provides around a 20x speedup for future calls of ``inspect.signature(cls)``.
    `Cosmology` has a lot of I/O methods that use the signature, so this is a
    significant speedup for those methods.

    Note that CPython does not promise that this precomputation is a stable feature.
    If it is removed, the worst that will happen is that the signature will be
    computed on the fly, the speedup will be lost, and this decorator can be
    deprecated.
    """
    cls.__signature__ = None  # clear the signature cache
    cls.__signature__ = inspect.signature(cls)  # add the new signature to the class
    return cls


def dataclass_decorator(cls):
    """Decorator for the dataclass transform.

    Returns
    -------
    cls : type
        The `cls` transformed into a frozen `~dataclasses.dataclass`.
        The ``__eq__`` method is custom (``eq=False``).
        The signature is precomputed and added to the class.
    """
    return _with_signature(dataclass(frozen=True, repr=True, eq=False, init=True)(cls))


##############################################################################


class CosmologyError(Exception):
    pass


# TODO: replace with `field(converter=lambda x: None if x is None else str(x))` when
#       the `converter` argument is available in `field` (py3.13, maybe?).
#       See https://peps.python.org/pep-0712/
@dataclass(frozen=True, slots=True)
class _NameField:
    default: str | None = None

    def __get__(self, instance: Cosmology | None, owner: type) -> str:
        # Called from the class. `dataclass` uses this to create ``__init__``.
        if instance is None:
            return self.default
        # Called from the instance
        return instance._name

    def __set__(self, instance: Cosmology, value: str | None) -> None:
        object.__setattr__(instance, "_name", (None if value is None else str(value)))


@dataclass_decorator
class Cosmology(metaclass=ABCMeta):
    """Base-class for all Cosmologies.

    Parameters
    ----------
    *args
        Arguments into the cosmology; used by subclasses, not this base class.
    name : str or None (optional, keyword-only)
        The name of the cosmology.
    meta : dict or None (optional, keyword-only)
        Metadata for the cosmology, e.g., a reference.
    **kwargs
        Arguments into the cosmology; used by subclasses, not this base class.

    Notes
    -----
    Class instances are static -- you cannot (and should not) change the
    values of the parameters.  That is, all of the above attributes
    (except meta) are read only.

    For details on how to create performant custom subclasses, see the
    documentation on :ref:`astropy-cosmology-fast-integrals`.

    Cosmology subclasses are automatically registered in a global registry
    and with various I/O methods. To turn off or change this registration,
    override the ``_register_cls`` classmethod in the subclass.
    """

    _: KW_ONLY
    name: _NameField = _NameField()
    """The name of the cosmology realization, e.g. 'Planck2018' or `None`."""

    meta: MetaData = MetaData()

    # Unified I/O object interchange methods
    from_format: ClassVar = UnifiedReadWriteMethod(CosmologyFromFormat)
    to_format: ClassVar = UnifiedReadWriteMethod(CosmologyToFormat)

    # Unified I/O read and write methods
    read: ClassVar = UnifiedReadWriteMethod(CosmologyRead)
    write: ClassVar = UnifiedReadWriteMethod(CosmologyWrite)

    # Parameters
    parameters = ParametersAttribute(attr_name="_parameters")
    """Immutable mapping of the Parameters.

    If accessed from the class, this returns a mapping of the Parameter
    objects themselves.  If accessed from an instance, this returns a
    mapping of the values of the Parameters.
    """

    _derived_parameters = ParametersAttribute(attr_name="_parameters_derived")
    """Immutable mapping of the derived Parameters.

    If accessed from the class, this returns a mapping of the Parameter
    objects themselves.  If accessed from an instance, this returns a
    mapping of the values of the Parameters.
    """

    _parameters: ClassVar = MappingProxyType[str, Parameter]({})
    _parameters_derived: ClassVar = MappingProxyType[str, Parameter]({})
    _parameters_all: ClassVar = frozenset[str]()

    __signature__: ClassVar[inspect.Signature | None] = None

    # ---------------------------------------------------------------

    def __init_subclass__(cls):
        super().__init_subclass__()

        # -------------------
        # Parameters

        all_params = all_parameters(cls)
        cls._parameters = MappingProxyType(
            {k: v for k, v in all_params.items() if not v.derived}
        )
        cls._derived_parameters = MappingProxyType(
            {k: v for k, v in all_params.items() if v.derived}
        )
        cls._parameters_all = frozenset(all_params)

        # -------------------
        # Registration

        if not inspect.isabstract(cls):  # skip abstract classes
            cls._register_cls()

    @classmethod
    def _register_cls(cls):
        # register class
        _COSMOLOGY_CLASSES[cls.__qualname__] = cls

        # register to YAML
        from astropy.cosmology._io.yaml import register_cosmology_yaml

        register_cosmology_yaml(cls)

    # ---------------------------------------------------------------

    def __post_init__(self):  # noqa: B027
        """Post-initialization, for subclasses to override if they need."""

    @property
    @abstractmethod
    def is_flat(self):
        """Return bool; `True` if the cosmology is flat.

        This is abstract and must be defined in subclasses.
        """
        raise NotImplementedError("is_flat is not implemented")

    def clone(self, *, meta=None, **kwargs):
        """Returns a copy of this object with updated parameters, as specified.

        This cannot be used to change the type of the cosmology, so ``clone()``
        cannot be used to change between flat and non-flat cosmologies.

        Parameters
        ----------
        meta : mapping or None (optional, keyword-only)
            Metadata that will update the current metadata.
        **kwargs
            Cosmology parameter (and name) modifications. If any parameter is
            changed and a new name is not given, the name will be set to "[old
            name] (modified)".

        Returns
        -------
        newcosmo : `~astropy.cosmology.Cosmology` subclass instance
            A new instance of this class with updated parameters as specified.
            If no arguments are given, then a reference to this object is
            returned instead of copy.

        Examples
        --------
        To make a copy of the ``Planck13`` cosmology with a different matter
        density (``Om0``), and a new name:

            >>> from astropy.cosmology import Planck13
            >>> Planck13.clone(name="Modified Planck 2013", Om0=0.35)
            FlatLambdaCDM(name='Modified Planck 2013', H0=<Quantity 67.77 km / (Mpc s)>,
                          Om0=0.35, ...

        If no name is specified, the new name will note the modification.

            >>> Planck13.clone(Om0=0.35).name
            'Planck13 (modified)'
        """
        # Quick return check, taking advantage of the Cosmology immutability.
        if meta is None and not kwargs:
            return self

        # There are changed parameter or metadata values.
        # The name needs to be changed accordingly, if it wasn't already.
        _modname = self.name + " (modified)"
        kwargs.setdefault("name", (_modname if self.name is not None else None))

        # mix new meta into existing, preferring the former.
        meta = meta if meta is not None else {}
        cloned = replace(self, meta=self.meta | meta, **kwargs)

        # Check if nothing has changed.
        # TODO! or should return self?
        if (cloned.name == _modname) and not meta and cloned.is_equivalent(self):
            object.__setattr__(cloned, "name", self.name)

        return cloned

    @classproperty
    def _init_has_kwargs(cls):
        return (
            next(reversed(cls.__signature__.parameters.values())).kind
            == inspect.Parameter.VAR_KEYWORD
        )

    # ---------------------------------------------------------------
    # comparison methods

    def is_equivalent(self, other: Any, /, *, format: _FormatType = False) -> bool:
        r"""Check equivalence between Cosmologies.

        Two cosmologies may be equivalent even if not the same class.
        For example, an instance of ``LambdaCDM`` might have :math:`\Omega_0=1`
        and :math:`\Omega_k=0` and therefore be flat, like ``FlatLambdaCDM``.

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance, positional-only
            The object to which to compare.
        format : bool or None or str, optional keyword-only
            Whether to allow, before equivalence is checked, the object to be
            converted to a |Cosmology|. This allows, e.g. a |Table| to be
            equivalent to a Cosmology.
            `False` (default) will not allow conversion. `True` or `None` will,
            and will use the auto-identification to try to infer the correct
            format. A `str` is assumed to be the correct format to use when
            converting.
            ``format`` is broadcast to match the shape of ``other``.
            Note that the cosmology arguments are not broadcast against
            ``format``, so it cannot determine the output shape.

        Returns
        -------
        bool
            True if cosmologies are equivalent, False otherwise.

        Examples
        --------
        Two cosmologies may be equivalent even if not of the same class.
        In this examples the ``LambdaCDM`` has ``Ode0`` set to the same value
        calculated in ``FlatLambdaCDM``.

            >>> import astropy.units as u
            >>> from astropy.cosmology import LambdaCDM, FlatLambdaCDM
            >>> cosmo1 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
            >>> cosmo2 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3)
            >>> cosmo1.is_equivalent(cosmo2)
            True

        While in this example, the cosmologies are not equivalent.

            >>> cosmo3 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, Tcmb0=3 * u.K)
            >>> cosmo3.is_equivalent(cosmo2)
            False

        Also, using the keyword argument, the notion of equivalence is extended
        to any Python object that can be converted to a |Cosmology|.

            >>> from astropy.cosmology import Planck18
            >>> tbl = Planck18.to_format("astropy.table")
            >>> Planck18.is_equivalent(tbl, format=True)
            True

        The list of valid formats, e.g. the |Table| in this example, may be
        checked with ``Cosmology.from_format.list_formats()``.

        As can be seen in the list of formats, not all formats can be
        auto-identified by ``Cosmology.from_format.registry``. Objects of
        these kinds can still be checked for equivalence, but the correct
        format string must be used.

            >>> tbl = Planck18.to_format("yaml")
            >>> Planck18.is_equivalent(tbl, format="yaml")
            True
        """
        from .funcs import cosmology_equal

        try:
            return cosmology_equal(
                self, other, format=(None, format), allow_equivalent=True
            )
        except Exception:
            # `is_equivalent` allows `other` to be any object and returns False
            # if `other` cannot be converted to a Cosmology, rather than
            # raising an Exception.
            return False

    def __equiv__(self, other: Any, /) -> bool:
        """Cosmology equivalence. Use ``.is_equivalent()`` for actual check!

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance, positional-only
            The object in which to compare.

        Returns
        -------
        bool or `NotImplemented`
            `NotImplemented` if ``other`` is from a different class.
            `True` if ``other`` is of the same class and has matching parameters
            and parameter values.
            `False` otherwise.
        """
        if other.__class__ is not self.__class__:
            return NotImplemented  # allows other.__equiv__

        # Check all parameters in 'other' match those in 'self' and 'other' has
        # no extra parameters (latter part should never happen b/c same class)
        # We do not use `self.parameters == other.parameters` because it does not work
        # for aggregating the truthiness of arrays, e.g. `m_nu`.
        return self._parameters_all == other._parameters_all and all(
            np.all(getattr(self, k) == getattr(other, k)) for k in self._parameters_all
        )

    def __eq__(self, other: Any, /) -> bool:
        """Check equality between Cosmologies.

        Checks the Parameters and immutable fields (i.e. not "meta").

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance, positional-only
            The object in which to compare.

        Returns
        -------
        bool
            `True` if Parameters and names are the same, `False` otherwise.
        """
        if other.__class__ is not self.__class__:
            return NotImplemented  # allows other.__eq__

        eq = (
            # non-Parameter checks: name
            self.name == other.name
            # check all parameters in 'other' match those in 'self' and 'other'
            # has no extra parameters (latter part should never happen b/c same
            # class). We do not use `self.parameters == other.parameters` because
            # it does not work for aggregating the truthiness of arrays, e.g. `m_nu`.
            # TODO! element-wise when there are array cosmologies
            and self._parameters_all == other._parameters_all
            and all(
                np.all(getattr(self, k) == getattr(other, k))
                for k in self._parameters_all
            )
        )

        return eq

    # ---------------------------------------------------------------

    def __str__(self):
        """Return a string representation of the cosmology."""
        name_str = "" if self.name is None else f'name="{self.name}", '
        param_strs = (f"{k!s}={v!s}" for k, v in self.parameters.items())
        return f"{type(self).__name__}({name_str}{', '.join(param_strs)})"

    def __astropy_table__(self, cls, copy, **kwargs):
        """Return a `~astropy.table.Table` of type ``cls``.

        Parameters
        ----------
        cls : type
            Astropy ``Table`` class or subclass.
        copy : bool
            Ignored.
        **kwargs : dict, optional
            Additional keyword arguments. Passed to ``self.to_format()``.
            See ``Cosmology.to_format.help("astropy.table")`` for allowed kwargs.

        Returns
        -------
        `astropy.table.Table` or subclass instance
            Instance of type ``cls``.
        """
        return self.to_format("astropy.table", cls=cls, **kwargs)


# Manipulate the dataclass fields for fields that are not created by `dataclass.fields`
# and thus do not have a mechanism for setting "compare", "repr", etc.
Cosmology.__dataclass_fields__["meta"].compare = False
Cosmology.__dataclass_fields__["meta"].repr = False


@dataclass_decorator
class FlatCosmologyMixin(metaclass=ABCMeta):
    """Mixin class for flat cosmologies.

    Do NOT instantiate directly. Note that all instances of
    ``FlatCosmologyMixin`` are flat, but not all flat cosmologies are
    instances of ``FlatCosmologyMixin``. As example, ``LambdaCDM`` **may**
    be flat (for the a specific set of parameter values), but
    ``FlatLambdaCDM`` **will** be flat.
    """

    _parameters: ClassVar[MappingProxyType[str, Parameter]]
    _parameters_derived: ClassVar[MappingProxyType[str, Parameter]]

    def __init_subclass__(cls: type[_FlatCosmoT]) -> None:
        super().__init_subclass__()

        # Determine the non-flat class.
        # This will raise a TypeError if the MRO is inconsistent.
        cls.__nonflatclass__  # noqa: B018

    # ===============================================================

    @classmethod  # TODO! make metaclass-method
    def _get_nonflat_cls(
        cls, kls: type[_CosmoT] | None = None
    ) -> type[Cosmology] | None:
        """Find the corresponding non-flat class.

        The class' bases are searched recursively.

        Parameters
        ----------
        kls : :class:`astropy.cosmology.Cosmology` class or None, optional
            If `None` (default) this class is searched instead of `kls`.

        Raises
        ------
        TypeError
            If more than one non-flat class is found at the same level of the
            inheritance. This is similar to the error normally raised by Python
            for an inconsistent method resolution order.

        Returns
        -------
        type
            A :class:`Cosmology` subclass this class inherits from that is not a
            :class:`FlatCosmologyMixin` subclass.
        """
        _kls = cls if kls is None else kls

        # Find non-flat classes
        nonflat: set[type[Cosmology]]
        nonflat = {
            b
            for b in _kls.__bases__
            if issubclass(b, Cosmology) and not issubclass(b, FlatCosmologyMixin)
        }

        if not nonflat:  # e.g. subclassing FlatLambdaCDM
            nonflat = {
                k for b in _kls.__bases__ if (k := cls._get_nonflat_cls(b)) is not None
            }

        if len(nonflat) > 1:
            raise TypeError(
                "cannot create a consistent non-flat class resolution order "
                f"for {_kls} with bases {nonflat} at the same inheritance level."
            )
        if not nonflat:  # e.g. FlatFLRWMixin(FlatCosmologyMixin)
            return None

        return nonflat.pop()

    __nonflatclass__ = classproperty(
        _get_nonflat_cls, lazy=True, doc="Return the corresponding non-flat class."
    )

    # ===============================================================

    @property
    def is_flat(self):
        """Return `True`, the cosmology is flat."""
        return True

    @property
    @abstractmethod
    def nonflat(self: _FlatCosmoT) -> _CosmoT:
        """Return the equivalent non-flat-class instance of this cosmology."""

    def clone(
        self, *, meta: Mapping | None = None, to_nonflat: bool = False, **kwargs
    ) -> Self:
        """Returns a copy of this object with updated parameters, as specified.

        This cannot be used to change the type of the cosmology, except for
        changing to the non-flat version of this cosmology.

        Parameters
        ----------
        meta : mapping or None (optional, keyword-only)
            Metadata that will update the current metadata.
        to_nonflat : bool, optional keyword-only
            Whether to change to the non-flat version of this cosmology.
        **kwargs
            Cosmology parameter (and name) modifications. If any parameter is
            changed and a new name is not given, the name will be set to "[old
            name] (modified)".

        Returns
        -------
        newcosmo : `~astropy.cosmology.Cosmology` subclass instance
            A new instance of this class with updated parameters as specified.
            If no arguments are given, then a reference to this object is
            returned instead of copy.

        Examples
        --------
        To make a copy of the ``Planck13`` cosmology with a different matter
        density (``Om0``), and a new name:

            >>> from astropy.cosmology import Planck13
            >>> Planck13.clone(name="Modified Planck 2013", Om0=0.35)
            FlatLambdaCDM(name='Modified Planck 2013', H0=<Quantity 67.77 km / (Mpc s)>,
                          Om0=0.35, ...

        If no name is specified, the new name will note the modification.

            >>> Planck13.clone(Om0=0.35).name
            'Planck13 (modified)'

        The keyword 'to_nonflat' can be used to clone on the non-flat equivalent
        cosmology. For :class:`~astropy.cosmology.FLRW` cosmologies this means
        ``Ode0`` can be modified:

            >>> Planck13.clone(to_nonflat=True, Ode0=1)
            LambdaCDM(name='Planck13 (modified)', H0=<Quantity 67.77 km / (Mpc s)>,
                      Om0=0.30712, Ode0=1.0, ...
        """
        if to_nonflat:
            return self.nonflat.clone(meta=meta, **kwargs)
        return super().clone(meta=meta, **kwargs)

    # ===============================================================

    def __equiv__(self, other):
        """Flat-|Cosmology| equivalence.

        Use `astropy.cosmology.funcs.cosmology_equal` with
        ``allow_equivalent=True`` for actual checks!

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance
            The object to which to compare for equivalence.

        Returns
        -------
        bool or `NotImplemented`
            `True` if ``other`` is of the same class / non-flat class (e.g.
            |FlatLambdaCDM| and |LambdaCDM|) has matching parameters and
            parameter values.
            `False` if ``other`` is of the same class but has different
            parameters.
            `NotImplemented` otherwise.
        """
        if isinstance(other, FlatCosmologyMixin):
            return super().__equiv__(other)  # super gets from Cosmology

        # check if `other` is the non-flat version of this class this makes the
        # assumption that any further subclass of a flat cosmo keeps the same
        # physics.
        if not issubclass(other.__class__, self.__nonflatclass__):
            return NotImplemented

        # Check if have equivalent parameters and all parameters in `other`
        # match those in `self`` and `other`` has no extra parameters.
        params_eq = (
            # no extra parameters
            self._parameters_all == other._parameters_all
            # equal
            and all(
                np.all(getattr(self, k) == getattr(other, k)) for k in self.parameters
            )
            # flatness check
            and other.is_flat
        )

        return params_eq
