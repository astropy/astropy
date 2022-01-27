# Licensed under a 3-clause BSD style license - see LICENSE.rst

import abc
import functools
import inspect
from types import FunctionType, MappingProxyType

import numpy as np

import astropy.units as u
from astropy.io.registry import UnifiedReadWriteMethod
from astropy.utils.decorators import classproperty
from astropy.utils.metadata import MetaData

from .connect import CosmologyFromFormat, CosmologyRead, CosmologyToFormat, CosmologyWrite
from .parameter import Parameter

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com), Roban Kramer
# (robanhk@gmail.com), and Nathaniel Starkman (n.starkman@mail.utoronto.ca).

# Many of these adapted from Hogg 1999, astro-ph/9905116
# and Linder 2003, PRL 90, 91301

__all__ = ["Cosmology", "CosmologyError", "FlatCosmologyMixin"]

__doctest_requires__ = {}  # needed until __getattr__ removed

# registry of cosmology classes with {key=name : value=class}
_COSMOLOGY_CLASSES = dict()


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
    __parameters__ = ()
    __all_parameters__ = ()

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
        ordered = [parameters.pop(parameters.index(n))
                   for n in cls._init_signature.parameters.keys()
                   if n in parameters]
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
            Cosmology parameter (and name) modifications.
            If any parameter is changed and a new name is not given, the name
            will be set to "[old name] (modified)".

        Returns
        -------
        newcosmo : `~astropy.cosmology.Cosmology` subclass instance
            A new instance of this class with updated parameters as specified.
            If no modifications are requested, then a reference to this object
            is returned instead of copy.

        Examples
        --------
        To make a copy of the ``Planck13`` cosmology with a different matter
        density (``Om0``), and a new name:

            >>> from astropy.cosmology import Planck13
            >>> newcosmo = Planck13.clone(name="Modified Planck 2013", Om0=0.35)

        If no name is specified, the new name will note the modification.

            >>> Planck13.clone(Om0=0.35).name
            'Planck13 (modified)'
        """
        # Quick return check, taking advantage of the Cosmology immutability.
        if meta is None and not kwargs:
            return self

        # There are changed parameter or metadata values.
        # The name needs to be changed accordingly, if it wasn't already.
        kwargs.setdefault("name", (self.name + " (modified)"
                                   if self.name is not None else None))

        # mix new meta into existing, preferring the former.
        new_meta = {**self.meta, **(meta or {})}
        # Mix kwargs into initial arguments, preferring the former.
        new_init = {**self._init_arguments, "meta": new_meta, **kwargs}
        # Create BoundArgument to handle args versus kwargs.
        # This also handles all errors from mismatched arguments
        ba = self._init_signature.bind_partial(**new_init)
        # Return new instance, respecting args vs kwargs
        return self.__class__(*ba.args, **ba.kwargs)

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

    def is_equivalent(self, other, *, format=False):
        r"""Check equivalence between Cosmologies.

        Two cosmologies may be equivalent even if not the same class.
        For example, an instance of ``LambdaCDM`` might have :math:`\Omega_0=1`
        and :math:`\Omega_k=0` and therefore be flat, like ``FlatLambdaCDM``.

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance or Any
            The object to check for equivalence.
        format : bool or None or str, optional keyword-only
            Whether to allow, before equivalence is checked, the object to be
            converted to a |Cosmology|. This allows, e.g. a |Table| to be
            equivalent to a Cosmology.
            `False` (default) will not allow conversion. `True` or `None` will,
            and will use the auto-identification to try to infer the correct
            format. A `str` is assumed to be the correct format to use when
            converting.

        Returns
        -------
        bool
            `True` if cosmologies are equivalent, `False` otherwise.

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

        Also, using the keyword argument ``format``, the notion of equivalence
        is extended to any Python object that can be converted to a |Cosmology|.

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
        # Allow for different formats to be considered equivalent.
        if format is not False:
            format = None if format is True else format  # str->str, None/True->None
            try:
                other = Cosmology.from_format(other, format=format)
            except Exception:  # TODO! should enforce only TypeError
                return False

        # The options are: 1) same class & parameters; 2) same class, different
        # parameters; 3) different classes, equivalent parameters; 4) different
        # classes, different parameters. (1) & (3) => True, (2) & (4) => False.
        equiv = self.__equiv__(other)
        if equiv is NotImplemented and hasattr(other, "__equiv__"):
            equiv = other.__equiv__(self)  # that failed, try from 'other'

        return equiv if equiv is not NotImplemented else False

    def __equiv__(self, other):
        """Cosmology equivalence. Use ``.is_equivalent()`` for actual check!

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance
            The object in which to compare.

        Returns
        -------
        bool or `NotImplemented`
            `NotImplemented` if 'other' is from a different class.
            `True` if 'other' is of the same class and has matching parameters
            and parameter values. `False` otherwise.
        """
        if other.__class__ is not self.__class__:
            return NotImplemented  # allows other.__equiv__

        # check all parameters in 'other' match those in 'self' and 'other' has
        # no extra parameters (latter part should never happen b/c same class)
        params_eq = (set(self.__all_parameters__) == set(other.__all_parameters__)
                     and all(np.all(getattr(self, k) == getattr(other, k))
                             for k in self.__all_parameters__))
        return params_eq

    def is_equal(self, other, *, check_meta=True, format=False):
        r"""Check equality between Cosmologies.

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance or Any
            The object to check for equality.
        check_meta : bool, optional keyword-only
            Whether to also check the metadata when determining equality.
        format : bool or None or str, optional keyword-only
            Whether to allow, before equality is checked, the object to be
            converted to a |Cosmology|. This allows, e.g. a |Table| to be
            equal to a Cosmology.
            `False` (default) will not allow conversion. `True` or `None` will,
            and will use the auto-identification to try to infer the correct
            format. A `str` is assumed to be the correct format to use when
            converting.

        Returns
        -------
        bool
            `True` if cosmologies are equal, `False` otherwise.

        See Also
        --------
        astropy.cosmology.Cosmology.is_equivalent
            Cosmologies may be equivalent, even if not the same class or name
            or metadata.

        Examples
        --------
        For a simple check that the two objects are |Cosmology|, with the same
        names and parameter values, use the standard python equality operator.

            >>> from astropy.cosmology import Planck13, Planck18
            >>> Planck18 == Planck18
            True

            >>> Planck13 != Planck18
            True

        ``is_equal`` also checks that the metadata (`astropy.Cosmology.meta`)
        are equal. To not check the metadata, set the keyword argument
        ``check_meta`` to `False`.

            >>> cosmo = Planck18.clone(name="Planck18", meta=dict(info="new"))
            >>> Planck18.is_equal(cosmo)
            False

            >>> Planck18.is_equal(cosmo, check_meta=False)
            True

        When the cosmologies have different names or parameter values they
        are never equal.

            >>> Planck18.is_equal(Planck13)
            False

        Using the keyword argument ``format``, the notion of equality is
        extended to any Python object that can be converted to a |Cosmology|.

            >>> tbl = Planck18.to_format("astropy.table")
            >>> Planck18.is_equal(tbl, format=True)
            True

        The list of valid formats, e.g. the |Table| in this example, may be
        checked with ``Cosmology.from_format.list_formats()``.

        As can be seen in the list of formats, not all formats can be
        auto-identified by ``Cosmology.from_format.registry``. Objects of
        these kinds can still be checked for equality, but the correct
        format string must be used.

            >>> tbl = Planck18.to_format("yaml")
            >>> Planck18.is_equal(tbl, format="yaml")
            True
        """
        # Allow for different formats to be considered equivalent.
        if format is not False:
            format = None if format is True else format  # str->str, None/True->None
            try:
                other = Cosmology.from_format(other, format=format)
            except Exception:  # TODO! should enforce only TypeError
                return False

        # Parameter equality
        eq = self.__eq__(other)
        if eq is NotImplemented and hasattr(other, "__eq__"):
            eq = other.__eq__(self)  # that failed, try from 'other'

        is_eq = eq if eq is not NotImplemented else False  # Ensure boolean

        # Metadata
        if check_meta and is_eq:  # Only check if required.
            is_eq &= (dict(self.meta) == dict(other.meta)) if hasattr(other, "meta") else False
            # using dict() to not care about ordering.

        return is_eq

    def __eq__(self, other):
        """Check equality between Cosmologies.

        Checks the Parameters and immutable fields (i.e. not "meta").

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance
            The object in which to compare.

        Returns
        -------
        bool
            `True` if Parameters and names are the same, `False` otherwise.
        """
        if other.__class__ is not self.__class__:
            return NotImplemented  # allows other.__eq__

        # check all parameters in 'other' match those in 'self'
        equivalent = self.__equiv__(other)
        # non-Parameter checks: name
        name_eq = (self.name == other.name)

        return equivalent and name_eq

    # ---------------------------------------------------------------

    def __repr__(self):
        ps = {k: getattr(self, k) for k in self.__parameters__}  # values
        cps = {k: getattr(self.__class__, k) for k in self.__parameters__}  # Parameter objects

        namelead = f"{self.__class__.__qualname__}("
        if self.name is not None:
            namelead += f"name=\"{self.name}\", "
        # nicely formatted parameters
        fmtps = (k + '=' + format(v, cps[k].format_spec if v is not None else '')
                 for k, v in ps.items())

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

    @property
    def is_flat(self):
        """Return `True`, the cosmology is flat."""
        return True


# -----------------------------------------------------------------------------


def __getattr__(attr):
    from . import flrw

    if hasattr(flrw, attr):
        import warnings

        from astropy.utils.exceptions import AstropyDeprecationWarning

        warnings.warn(
            f"`astropy.cosmology.core.{attr}` has been moved (since v5.0) and "
            f"should be imported as ``from astropy.cosmology import {attr}``."
            " In future this will raise an exception.",
            AstropyDeprecationWarning
        )

        return getattr(flrw, attr)

    raise AttributeError(f"module {__name__!r} has no attribute {attr!r}.")
