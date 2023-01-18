# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import abc
import inspect
from typing import TYPE_CHECKING, Any, TypeVar

import numpy as np

from astropy.io.registry import UnifiedReadWriteMethod
from astropy.utils.decorators import classproperty
from astropy.utils.metadata import MetaData

from .connect import (
    CosmologyFromFormat,
    CosmologyRead,
    CosmologyToFormat,
    CosmologyWrite,
)
from .parameter import Parameter

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Mapping

    from astropy.cosmology.funcs.comparison import _FormatType

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com), Roban Kramer
# (robanhk@gmail.com), and Nathaniel Starkman (n.starkman@mail.utoronto.ca).

# Many of these adapted from Hogg 1999, astro-ph/9905116
# and Linder 2003, PRL 90, 91301

__all__ = ["Cosmology", "CosmologyError", "FlatCosmologyMixin"]

__doctest_requires__ = {}  # needed until __getattr__ removed


##############################################################################
# Parameters

# registry of cosmology classes with {key=name : value=class}
_COSMOLOGY_CLASSES = dict()

# typing
_CosmoT = TypeVar("_CosmoT", bound="Cosmology")
_FlatCosmoT = TypeVar("_FlatCosmoT", bound="FlatCosmologyMixin")

##############################################################################


class CosmologyError(Exception):
    pass


class Cosmology(metaclass=abc.ABCMeta):
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
    Class instances are static -- you cannot (and should not) change the values
    of the parameters.  That is, all of the above attributes (except meta) are
    read only.

    For details on how to create performant custom subclasses, see the
    documentation on :ref:`astropy-cosmology-fast-integrals`.
    """

    meta = MetaData()

    # Unified I/O object interchange methods
    from_format = UnifiedReadWriteMethod(CosmologyFromFormat)
    to_format = UnifiedReadWriteMethod(CosmologyToFormat)

    # Unified I/O read and write methods
    read = UnifiedReadWriteMethod(CosmologyRead)
    write = UnifiedReadWriteMethod(CosmologyWrite)

    # Parameters
    __parameters__: tuple[str, ...] = ()
    __all_parameters__: tuple[str, ...] = ()

    # ---------------------------------------------------------------

    def __init_subclass__(cls):
        super().__init_subclass__()

        # -------------------
        # Parameters

        # Get parameters that are still Parameters, either in this class or above.
        parameters = []
        derived_parameters = []
        for n in cls.__parameters__:
            p = getattr(cls, n)
            if isinstance(p, Parameter):
                derived_parameters.append(n) if p.derived else parameters.append(n)

        # Add new parameter definitions
        for n, v in cls.__dict__.items():
            if n in parameters or n.startswith("_") or not isinstance(v, Parameter):
                continue
            derived_parameters.append(n) if v.derived else parameters.append(n)

        # reorder to match signature
        ordered = [
            parameters.pop(parameters.index(n))
            for n in cls._init_signature.parameters.keys()
            if n in parameters
        ]
        parameters = ordered + parameters  # place "unordered" at the end
        cls.__parameters__ = tuple(parameters)
        cls.__all_parameters__ = cls.__parameters__ + tuple(derived_parameters)

        # -------------------
        # register as a Cosmology subclass
        _COSMOLOGY_CLASSES[cls.__qualname__] = cls

    @classproperty(lazy=True)
    def _init_signature(cls):
        """Initialization signature (without 'self')."""
        # get signature, dropping "self" by taking arguments [1:]
        sig = inspect.signature(cls.__init__)
        sig = sig.replace(parameters=list(sig.parameters.values())[1:])
        return sig

    # ---------------------------------------------------------------

    def __init__(self, name=None, meta=None):
        self._name = str(name) if name is not None else name
        self.meta.update(meta or {})

    @property
    def name(self):
        """The name of the Cosmology instance."""
        return self._name

    @property
    @abc.abstractmethod
    def is_flat(self):
        """
        Return bool; `True` if the cosmology is flat.
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
            FlatLambdaCDM(name="Modified Planck 2013", H0=67.77 km / (Mpc s),
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
        new_meta = {**self.meta, **meta}
        # Mix kwargs into initial arguments, preferring the former.
        new_init = {**self._init_arguments, "meta": new_meta, **kwargs}
        # Create BoundArgument to handle args versus kwargs.
        # This also handles all errors from mismatched arguments
        ba = self._init_signature.bind_partial(**new_init)
        # Instantiate, respecting args vs kwargs
        cloned = type(self)(*ba.args, **ba.kwargs)

        # Check if nothing has changed.
        # TODO! or should return self?
        if (cloned.name == _modname) and not meta and cloned.is_equivalent(self):
            cloned._name = self.name

        return cloned

    @property
    def _init_arguments(self):
        # parameters
        kw = {n: getattr(self, n) for n in self.__parameters__}

        # other info
        kw["name"] = self.name
        kw["meta"] = self.meta

        return kw

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
        return set(self.__all_parameters__) == set(other.__all_parameters__) and all(
            np.all(getattr(self, k) == getattr(other, k))
            for k in self.__all_parameters__
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
            # class) TODO! element-wise when there are array cosmologies
            and set(self.__all_parameters__) == set(other.__all_parameters__)
            and all(
                np.all(getattr(self, k) == getattr(other, k))
                for k in self.__all_parameters__
            )
        )

        return eq

    # ---------------------------------------------------------------

    def __repr__(self):
        namelead = f"{self.__class__.__qualname__}("
        if self.name is not None:
            namelead += f'name="{self.name}", '
        # nicely formatted parameters
        fmtps = (f"{k}={getattr(self, k)}" for k in self.__parameters__)

        return namelead + ", ".join(fmtps) + ")"

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


class FlatCosmologyMixin(metaclass=abc.ABCMeta):
    """
    Mixin class for flat cosmologies. Do NOT instantiate directly.
    Note that all instances of ``FlatCosmologyMixin`` are flat, but not all
    flat cosmologies are instances of ``FlatCosmologyMixin``. As example,
    ``LambdaCDM`` **may** be flat (for the a specific set of parameter values),
    but ``FlatLambdaCDM`` **will** be flat.
    """

    __all_parameters__: tuple[str, ...]
    __parameters__: tuple[str, ...]

    def __init_subclass__(cls: type[_FlatCosmoT]) -> None:
        super().__init_subclass__()

        # Determine the non-flat class.
        # This will raise a TypeError if the MRO is inconsistent.
        cls.__nonflatclass__

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

    @abc.abstractmethod
    def nonflat(self: _FlatCosmoT) -> _CosmoT:
        """Return the equivalent non-flat-class instance of this cosmology."""

    def clone(self, *, meta: Mapping | None = None, to_nonflat: bool = False, **kwargs):
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
            FlatLambdaCDM(name="Modified Planck 2013", H0=67.77 km / (Mpc s),
                          Om0=0.35, ...

        If no name is specified, the new name will note the modification.

            >>> Planck13.clone(Om0=0.35).name
            'Planck13 (modified)'

        The keyword 'to_nonflat' can be used to clone on the non-flat equivalent
        cosmology.

            >>> Planck13.clone(to_nonflat=True)
            LambdaCDM(name="Planck13", ...

            >>> Planck13.clone(H0=70, to_nonflat=True)
            LambdaCDM(name="Planck13 (modified)", H0=70.0 km / (Mpc s), ...
        """
        if to_nonflat:
            return self.nonflat.clone(meta=meta, **kwargs)
        return super().clone(meta=meta, **kwargs)

    # ===============================================================

    def __equiv__(self, other):
        """flat-|Cosmology| equivalence.

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
            set(self.__all_parameters__) == set(other.__all_parameters__)  # no extra
            # equal
            and all(
                np.all(getattr(self, k) == getattr(other, k))
                for k in self.__parameters__
            )
            # flatness check
            and other.is_flat
        )

        return params_eq


# -----------------------------------------------------------------------------


def __getattr__(attr):
    from . import flrw

    if hasattr(flrw, attr) and attr not in ("__path__",):
        import warnings

        from astropy.utils.exceptions import AstropyDeprecationWarning

        warnings.warn(
            f"`astropy.cosmology.core.{attr}` has been moved (since v5.0) and "
            f"should be imported as ``from astropy.cosmology import {attr}``."
            " In future this will raise an exception.",
            AstropyDeprecationWarning,
        )

        return getattr(flrw, attr)

    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}.")
