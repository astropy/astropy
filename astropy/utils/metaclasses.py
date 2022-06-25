# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy's ABCs."""

from __future__ import annotations
from ast import Raise
import inspect

import sys
import weakref
from abc import ABCMeta
from typing import Sequence, Dict, Any, Tuple, TypeVar, Union
from astropy.utils.compat.optional_deps import HAS_TYPING_EXTENSIONS

if sys.version_info >= (3, 11):
    from typing import Self
elif HAS_TYPING_EXTENSIONS:
    from typing_extensions import Self
else:
    Self = TypeVar("Self", bound="FactoryMeta")


__all__: list[str] = []  # FactoryMeta is in development


TypeArgs = Tuple[str, Sequence[type], Dict[str, Any]]


class FactoryMeta(ABCMeta):
    """Meta-class for generating classes of a data class.

    Classes derived from this metaclass.

    Parameters
    ----------
    name: str
        The name of the class. See `type` for details.
    bases: tuple[type, ...]
        The bases of the class. See `type` for details.
    namespace: dict[str, Any]
        The namespace of the class. See `type` for details.
    data_cls : type or None, optional keyword-only
        The un-mixed type when defining an intermixed class.
        Passing the ``data_cls`` means the generted intermixed class can be
        reused.

    base_cls : type or None, optional keyword-only
        The un-mixed type when defining an intermixed class that should be
        used as a base class for future generated intermixed classes.
        This is useful when a type is a base class for many derived classes,
        e.g. `numpy.ndarray` for the entire ``MaskedNDArray`` family.

    Warnings
    --------
    This is an experimental feature and subject to change. In particular, the
    private methods are not stable.

    Examples
    --------
    For real examples implemented in the code, see `astropy.utils.masked.Masked`
    or `astropy.uncertainty.Distribution`. Here we will work through a toy
    problem.

    Let's say we like how `numpy.ndarray` can cast elements of an array to
    different datatypes and we want to add this functionality to some builtin
    types. This starts by defining the inter-mix base class, using
    `FactoryMeta` as a metaclass. The ``_inmix_make_instance`` defines
    how new instances are constructed. In this case it's just calling the same
    args and kwargs. We also add a method ``astype`` for casting each element of
    self to a float.

        >>> class SupportsAsType(metaclass=FactoryMeta):
        ...     @classmethod
        ...     def _inmix_make_instance(cls, data, *args, **kwargs):
        ...         inmixcls = cls._inmix_get_subclass(type(data))
        ...         return inmixcls(data, *args, **kwargs)
        ...     def astype(self, dtype):
        ...         return type(self)([dtype(x) for x in self])

    Now ``SupportsAsType`` can be used to generate new classes or constuct
    instances of those classes.

    The most straightforward example is to subclass``SupportsAsType`` and
    inherit from the desired type.

        >>> class SupportsAsTypeList(list, SupportsAsType, data_cls=list):
        ...     pass

    ``SupportsAsTypeList`` can make list-like objects...

        >>> inst = SupportsAsTypeList([1, 2, 3])
        >>> inst
        [1, 2, 3]

    with the added method ``astype``.

        >>> inst.astype(float)
        [1.0, 2.0, 3.0]

    ``SupportsAsType`` can also be used as a class factory:

        >>> cls = SupportsAsType(list)
        >>> cls
        <class 'astropy.utils.metaclasses.SupportsAsTypeList'>

    When defining ``SupportsAsTypeList``, the class keyword-argument
    ``data_cls`` registered `list`, so ``SupportsAsTypeList`` will be stored for
    re-use.

        >>> cls is SupportsAsTypeList
        True

    This is another example of using ``SupportsAsType`` as a class factory:

        >>> cls = SupportsAsType(tuple)
        >>> cls
        <class 'abc.tupleSupportsAsType'>

        >>> inst = cls([1, 2, 3])
        >>> inst
        (1, 2, 3)
        >>> inst.astype(float)
        (1.0, 2.0, 3.0)

    ``SupportsAsType`` can also be used as an instance factory, utilizing the
    class-factory functionality in the background.

        >>> inst = SupportsAsType({1, })
        >>> type(inst)
        <class 'abc.setSupportsAsType'>
        >>> inst
        setSupportsAsType({1})
        >>> inst.astype(float)
        setSupportsAsType({1.0})
    """

    # ---------------------------------------------------------------
    # Construct class from metaclass

    __inmixbase__: type

    def __new__(
        cls: type[Self],
        name: str,
        bases: tuple[type, ...],
        namespace: dict[str, Any],
        base_cls: type | None=None,
        data_cls: type | None=None,
        **_: Any
    ) -> Self:
        # Make a new class. `__new__` is required to prevent class-level kwargs
        # from going to the superclass (ABCMeta), which does not take any
        # kwargs. Also, `__init_subclass__`, which may be defined on the classes
        # `FactoryMeta` is trying to make, will not have access to certain
        # class-level attributes that are made during the class creation
        # process.

        # inject data_cls into namespace before `__new__` so
        # `__init_subclass__` can use it.
        namespace["__intomixclass__"] = data_cls

        self = super().__new__(cls, name, bases, namespace)

        # check that this `cls` is the top-most in the MRO to be a
        # `FactoryMeta` type.
        # TODO! allow for subclasses of a
        # FactoryMeta class to be considered a baseclass.
        inmixbase = not any([isinstance(base, FactoryMeta) for base in bases])
        if inmixbase:
            self.__inmixbase__ = self
            """Base class for generating new intermixed classes."""

            self.__inmixed_base_classes = weakref.WeakKeyDictionary()
            """See ``_inmixed_base_classes``."""

            self.__inmixed_generated_classes = weakref.WeakKeyDictionary()
            """See ``_inmixed_generated_classes``."""

        elif not inspect.isclass(getattr(self, "__inmixbase__", None)):
            raise TypeError(f"type {self} does not have a defined '__inmixbase__'")

        return self

    def __init__(
        self: Self,
        name: str,
        bases: tuple[type, ...],
        namespace: dict[str, Any],
        base_cls: type | None=None,
        data_cls: type | None=None,
        **_: Any
    ) -> None:
        # Optionally register base class and data class information. This can be
        # done when creating the base inmix class, or subclasses thereof.
        if base_cls is not None:
            self.__inmixbase__.__inmixed_base_classes[base_cls] = weakref.ref(self)

        self.__intomixclass__: type | None = None
        if data_cls is not None:
            self.__intomixclass__ = data_cls  # (also injected into namespace in __new__)
            self.__inmixbase__.__inmixed_generated_classes[data_cls] = weakref.ref(self)

            if self.__doc__ is None:
                self.__doc__ = self._inmix_make__doc__(data_cls)

        super().__init__(name, bases, namespace)

    def _inmix_make__doc__(self, data_cls: type, /) -> str | None:
        """
        Make the docstring for the intermixed class.
        The default is to use the docstring of the data class.
        """
        return data_cls.__doc__

    # ---------------------------------------------------------------
    # Construct instance from class
    # This is when an intermixed class can be generated.

    def __call__(self: Self, *args: Any, **kwargs: Any) -> Self | Any:
        # Make an instance of the class. Before ``cls.__new__``` is called this
        # decides whether 1 of 2 things happens:
        # 1) Normal class instantiation proceeds
        # 2) "Factory" mode: defining either a new wrapped class or a making
        #    an instance of a wrapped class.
        if self is self.__inmixbase__:  # "Factory" mode
            # Defining a new class -- a sub-class of the argument, with an
            # injected inheritance.
            if not kwargs and len(args) == 1 and isinstance(args[0], type):
                return self._inmix_get_subclass(args[0])
            # Making an instance of a dependency-injected class.
            else:
                return self._inmix_make_instance(*args, **kwargs)

        # Normal class instantiation
        return super().__call__(*args, **kwargs)

    def _inmix_get_subclass(self: Self, data_cls: type, /) -> Self:
        """Get the inmix class for a given data class.

        Parameters
        ----------
        data_cls : type
            The given data class.

        Returns
        -------
        type
            The generated class.
        """
        # If the class is already correct, return unchanged.
        if issubclass(data_cls, self.__inmixbase__):
            return data_cls

        # Start by looking through the generated classes.
        inmixed = self._inmixed_generated_classes_get(data_cls)

        # If no generated class is good, a new class will have to be generated
        # (or the default class used).
        if inmixed is None:
            # Walk through MRO and find closest base data class. its own entry.
            for mro_item in data_cls.__mro__:
                base_cls = self._inmixed_base_classes_get(mro_item)
                if base_cls is not None:
                    break
            else:  # Hopefully this class is a good base class.
                if self.inmixed_default_cls is not None:
                    return self.inmixed_default_cls
                else:
                    base_cls = self

            # Create (and therefore register) new inmix subclass for the given
            # data_cls. Generally in the MRO the inmix class is inserted below
            # the data class, but the specifics are controlled by
            # `_inmix_prepare_type`.
            name, bases, namespace = self._inmix_prepare_type(data_cls, base_cls)
            inmixed = type(name, bases, namespace, data_cls=data_cls)

        return inmixed

    def _inmix_prepare_type(self, data_cls: type, base_cls: type, /) -> TypeArgs:
        """Prepare a subclass of the inmix class.

        Generally in the MRO the inmix class is after the data class.

        Parameters
        ----------
        data_cls : type
            The class of the data, e.g. `~numpy.ndarray`.
        base_cls : type
            Wrapper class, a subclass of the base inmix class.

        Returns
        -------
        name : str
            The name of the class that to be generated.
        bases : tuple[type, ...]
            The base class for the class to be generated.
        namespace : dict[str, Any]
            The namespace for `type`.

        See Also
        --------
        type
            For details of the arguments
        astropy.utils.metaclasses.FactoryMeta._inmix_get_subclass
            For where this method is called.
        """
        name: str = data_cls.__name__ + self.__name__
        bases: tuple[type, ...] = (data_cls, base_cls)
        return name, bases, {}

    def _inmix_make_instance(self: Self, data: Any, /,*args: Any, **kwargs: Any) -> Any:
        """Make an instance of the inmix subclass for the data type.

        Parameters
        ----------
        data : Any
            The data for which to make an instance of the relevant inmix
            subclass.
        *args, **kwargs : Any
            Arguments, positional and keyword, creating an instance of the
            inmix subclass for the data type.

        Returns
        -------
        object
        """
        inmixedcls = self._inmix_get_subclass(type(data))
        return inmixedcls(data, *args, **kwargs)

    # ---------------------------------------------------------------
    # Intermixed base classes

    @property
    def _inmixed_base_classes(self) -> weakref.WeakKeyDictionary:
        """Intermixed classes keyed by their unmixed counterparts.

        Base classes are instances of this metaclass that can serve as a base
        class when defining an intermixed class for a given data class.

        Returns
        -------
        `weakref.WeakKeyDictionary`
            If an unmixed class is deleted and garbage collected then the
            defined subclass will be removed as well. The values are
            `weakref.ReferenceType`, so are similarly weakly connected to the
            value, but keys are not removed if the value is missing. For proper
            weak key and value usage use ``_inmixed_base_classes_get``.

        See Also
        --------
        astropy.utils.metaclasses.FactoryMeta._inmixed_generated_classes
            Intermixed classes for a specific data class; they are not used as a
            base class for further generated subclasses.

        Notes
        -----
        For subclasses of these unmixed classes, intermixed counterparts can be
        generated and will be stored in ``_inmixed_generated_classes``.
        """
        return self.__inmixbase__.__inmixed_base_classes

    def _inmixed_base_classes_get(self, base_cls: type) -> type | None:
        """Get an intermixed class given its unmixed base class.

        Parameters
        ----------
        base_cls : type
            The unmixed class.

        Returns
        -------
        type or None
            The value corresponding to ``base_cls``. If `None`, ``base_cls`` is
            popped from the `weakref.WeakKeyDictionary`.
        """
        # Get weak reference
        inmixed_ref = self._inmixed_base_classes.get(base_cls)
        # Resolve reference
        inmixed = inmixed_ref() if isinstance(inmixed_ref, weakref.ReferenceType) else None
        # remove dead references
        if inmixed is None:
            self._inmixed_base_classes.pop(base_cls, None)

        return inmixed

    @property
    def inmixed_default_cls(self: Self) -> Self | None:
        """The default inmix base class.

        Returns
        -------
        type

        Notes
        -----
        In the reference implementation (`astropy.utils.masked.Masked`) this is
        the subclass of `astropy.utils.masked.Masked` for masking
        `numpy.ndarray` and array-like objects.
        """
        return None

    # ---------------------------------------------------------------
    # Intermixed generated classes

    @property
    def _inmixed_generated_classes(self) -> weakref.WeakKeyDictionary:
        """Intermixed classes keyed by their unmixed counterparts.

        Returns
        -------
        `~weakref.WeakKeyDictionary`
            If an unmixed class is deleted and garbage collected then the
            defined subclass will be removed as well. The values are
            `weakref.ReferenceType`, so are similarly weakly connected to the
            value, but keys are not removed if the value is missing. For proper
            weak key and value usage use ``_inmixed_generated_classes_get``.

        See Also
        --------
        astropy.utils.metaclasses.FactoryMeta._inmixed_base_classes
            Intermixed classes that can serve as a base class for another
            intermixed class.
        """
        return self.__inmixbase__.__inmixed_generated_classes

    def _inmixed_generated_classes_get(self, data_cls: type) -> type | None:
        """Get an intermixed class given its unmixed data class.

        Parameters
        ----------
        data_cls : type
            The unmixed class.

        Returns
        -------
        type or None
            The value corresponding to ``data_cls``. If `None`, ``data_cls`` is
            popped from the `weakref.WeakKeyDictionary`.
        """
        # Get weak reference
        inmixed_ref = self._inmixed_generated_classes.get(data_cls)
        # Resolve reference
        inmixed = inmixed_ref() if isinstance(inmixed_ref, weakref.ReferenceType) else None
        # remove dead references
        if inmixed is None:
            self._inmixed_generated_classes.pop(data_cls, None)

        return inmixed
