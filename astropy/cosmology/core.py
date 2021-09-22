# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
from abc import ABCMeta
from types import MappingProxyType

import numpy as np

import astropy.units as u
from astropy.io.registry import UnifiedReadWriteMethod
from astropy.utils import isiterable
from astropy.utils.decorators import classproperty
from astropy.utils.metadata import MetaData

from .connect import CosmologyFromFormat, CosmologyRead, CosmologyToFormat, CosmologyWrite

# Originally authored by Andrew Becker (becker@astro.washington.edu),
# and modified by Neil Crighton (neilcrighton@gmail.com) and Roban
# Kramer (robanhk@gmail.com).

# Many of these adapted from Hogg 1999, astro-ph/9905116
# and Linder 2003, PRL 90, 91301

__all__ = ["Cosmology", "CosmologyError", "FlatCosmologyMixin", "Parameter"]

__doctest_requires__ = {}  # needed until __getattr__ removed

# registry of cosmology classes with {key=name : value=class}
_COSMOLOGY_CLASSES = dict()


class CosmologyError(Exception):
    pass


class Parameter:
    """Cosmological parameter (descriptor).

    Should only be used with a :class:`~astropy.cosmology.Cosmology` subclass.

    Parameters
    ----------
    fget : callable or None, optional
        Function to get the value from instances of the cosmology class.
        If None (default) returns the corresponding private attribute.
        Often not set here, but as a decorator with ``getter``.
    doc : str or None, optional
        Parameter description. If 'doc' is None and 'fget' is not, then 'doc'
        is taken from ``fget.__doc__``.

    Examples
    --------
    The most common use case of ``Parameter`` is to access the corresponding
    private attribute.

        >>> from astropy.cosmology import LambdaCDM
        >>> from astropy.cosmology.core import Parameter
        >>> class Example1(LambdaCDM):
        ...     param = Parameter(doc="example parameter")
        ...     def __init__(self, param=15):
        ...         self._param = param
        >>> Example1.param
        <Parameter 'param' at ...
        >>> Example1.param.default
        15

        >>> ex = Example1(param=12357)
        >>> ex.param
        12357

    ``Parameter`` also supports custom ``getter`` methods.

        >>> import astropy.units as u
        >>> class Example2(LambdaCDM):
        ...     param = Parameter(doc="example parameter")
        ...     def __init__(self, param=15):
        ...         self._param = param
        ...     @param.getter
        ...     def param(self):
        ...         return self._param * u.m

        >>> ex2 = Example2(param=12357)
        >>> ex2.param
        <Quantity 12357. m>

    .. doctest::
       :hide:

       >>> from astropy.cosmology.core import _COSMOLOGY_CLASSES
       >>> _ = _COSMOLOGY_CLASSES.pop(Example1.__qualname__)
       >>> _ = _COSMOLOGY_CLASSES.pop(Example2.__qualname__)
    """

    def __init__(self, fget=None, doc=None):
        # modeled after https://docs.python.org/3/howto/descriptor.html#properties
        self.__doc__ = fget.__doc__ if (doc is None and fget is not None) else doc
        self.fget = fget if not hasattr(fget, "fget") else fget.__get__
        # TODO! better detection if descriptor.

        # misc
        self.__wrapped__ = fget  # so always have access to `fget`
        self.__name__ = getattr(fget, "__name__", None)  # compat with other descriptors

    def __set_name__(self, objcls, name):
        # attribute name
        self._attr_name = name
        self._attr_name_private = "_" + name

    @property
    def name(self):
        """Parameter name."""
        return self._attr_name

    # -------------------------------------------
    # descriptor

    def __get__(self, obj, objcls=None):
        # get from class
        if obj is None:
            return self
        # get from obj, allowing for custom ``getter``
        if self.fget is None:  # default to private attr (diff from `property`)
            return getattr(obj, self._attr_name_private)
        return self.fget(obj)

    def __set__(self, obj, value):
        raise AttributeError("can't set attribute")

    def __delete__(self, obj):
        raise AttributeError("can't delete attribute")

    # -------------------------------------------
    # from 'property'

    def getter(self, fget):
        """Make new Parameter with custom ``fget``.

        Parameters
        ----------
        fget : callable

        Returns
        -------
        `~astropy.cosmology.Parameter`
            Copy of this Parameter but with custom ``fget``.
        """
        return type(self)(fget=fget, doc=self.__doc__)

    # -------------------------------------------

    def __repr__(self):
        return f"<Parameter {self._attr_name!r} at {hex(id(self))}>"


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

    def __init_subclass__(cls):
        super().__init_subclass__()

        # names of Parameters
        cls.__parameters__ = frozenset({
            n for n in dir(cls)
            # select only <Parameter>, and guard from private attributes
            if (not n.startswith("_") and isinstance(getattr(cls, n, None), Parameter))
        })

        _COSMOLOGY_CLASSES[cls.__qualname__] = cls

    def __new__(cls, *args, **kwargs):
        self = super().__new__(cls)

        # bundle and store initialization arguments on the instance
        ba = cls._init_signature.bind_partial(*args, **kwargs)
        ba.apply_defaults()  # and fill in the defaults
        self._init_arguments = ba.arguments

        return self

    def __init__(self, *args, name=None, meta=None, **kwargs):
        self._name = name
        self.meta.update(meta or {})

    @classproperty(lazy=True)
    def _init_signature(cls):
        """Initialization signature (without 'self')."""
        # get signature, dropping "self" by taking arguments [1:]
        sig = inspect.signature(cls.__init__)
        sig = sig.replace(parameters=list(sig.parameters.values())[1:])
        return sig

    @property
    def name(self):
        """The name of the Cosmology instance."""
        return self._name

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

    # -----------------------------------------------------

    def __eq__(self, other):
        """Check equality on all immutable fields (i.e. not "meta").

        Parameters
        ----------
        other : `~astropy.cosmology.Cosmology` subclass instance
            The object in which to compare.

        Returns
        -------
        bool
            True if all immutable fields are equal, False otherwise.
        """
        if not isinstance(other, Cosmology):
            return False

        sias, oias = self._init_arguments, other._init_arguments

        # check if the cosmologies have identical signatures.
        # this protects against one cosmology having a superset of input
        # parameters to another cosmology.
        if (sias.keys() ^ oias.keys()) - {'meta'}:
            return False

        # are all the non-excluded immutable arguments equal?
        return all((np.all(oias[k] == v) for k, v in sias.items()
                    if k != "meta"))


class FlatCosmologyMixin(metaclass=ABCMeta):
    """
    Mixin class for flat cosmologies. Do NOT instantiate directly.
    Note that all instances of ``FlatCosmologyMixin`` are flat, but not all
    flat cosmologies are instances of ``FlatCosmologyMixin``. As example,
    ``LambdaCDM`` **may** be flat (for the a specific set of parameter values),
    but ``FlatLambdaCDM`` **will** be flat.
    """


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
