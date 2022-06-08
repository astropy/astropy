# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy's ABCs."""


import weakref
from abc import ABCMeta, abstractmethod
from typing import Tuple, Type, Dict, Any, Optional, TypeVar, Union

try:
    from typing import Self
except ImportError:
    try:
        from typing_extensions import Self
    except ImportError:
        Self = TypeVar("Self", bound="InheritanceInMixMeta")


__all__ = []  # InheritanceInMixMeta is in development


class InheritanceInMixMeta(ABCMeta):
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
    base_cls : type or None, optional keyword-only
        The un-mixed type when defining an intermixed class that should be
        used as a base class for future generated intermixed classes.
        This is useful when a type is base class for many derived classes,
        e.g. `numpy.ndarray` for the entire `astropy.units.Quantity` family.
    data_cls : type or None, optional keyword-only
        The un-mixed type when defining an intermixed class.
    **kwargs : Any
        A valid option is ``default_inmixed_class``. This should only be given
        to the base wrapper class. In Examples this is ``SupportsAsType``.

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
    `InheritanceInMixMeta` as a metaclass. The ``_inmix_make_instance`` defines
    how new instances are constructed. In this case it's just calling the same
    args and kwargs. We also add a method ``astype`` for casting each element of
    self to a float.

        >>> class SupportsAsType(metaclass=InheritanceInMixMeta):
        ...     @classmethod
        ...     def _inmix_make_instance(cls, data, *args, **kwargs):
        ...         inmixcls = cls._inmix_make_class(type(data))
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

    def __new__(cls: Type[Self], name: str, bases: Tuple[type, ...],
                namespace: Dict[str, Any], base_cls: Optional[type]=None,
                data_cls: Optional[type]=None, **kwargs: Any) -> Self:
        # Make a new class. `__new__` is required to prevent class-level
        # kwargs from going to the superclass (ABCMeta), which does not take
        # any kwargs. Also, `__init_subclass__`, which may be defined on the
        # classes `InheritanceInMixMeta` is trying to make, will not have access
        # to certain class-level attributes that are made during the class
        # creation process.

        # ABCMeta takes no kwargs. This is still passed to InheritanceInMixMeta.__init__.
        kwargs.pop("default_inmixed_class", None)

        # inject data_cls into namespace before `__new__` so
        # `__init_subclass__` can use it.
        if data_cls is not None:
            namespace["__intomixclass__"] = data_cls

        return super().__new__(cls, name, bases, namespace, **kwargs)

    def __init__(self: Self, name: str, bases: Tuple[type, ...],
                 namespace: Dict[str, Any], base_cls: Optional[type]=None,
                 data_cls: Optional[type]=None, **kwargs: Any) -> None:
        # check that this `cls` is the top-most in the MRO to be a
        # `InheritanceInMixMeta` type TODO! allow for subclasses of a
        # InheritanceInMixMeta class to be considered a baseclass.
        isbaseclass = all([not issubclass(base.__class__, InheritanceInMixMeta) for base in bases])
        if isbaseclass:
            self.__inmixbase__ = self
            """Base class for generating new intermixed classes."""

            self.__inmixed_default_defined_class: Union[str, Self, None]
            self.__inmixed_default_defined_class = kwargs.pop("default_inmixed_class", None)
            """
            Fallback wrapper class when there is not a good base class.
            See ``_inmixed_predefined_classes``.
            """

            self.__inmixed_predefined_classes = weakref.WeakKeyDictionary()
            """Explicitly defined intermixed classes keyed by their unmixed counterparts.

            For subclasses of these unmixed classes, intermixed counterparts can
            be generated and will be stored in ``_inmixed_generated_classes``.
            """

            self.__inmixed_generated_classes = weakref.WeakKeyDictionary()
            """Inter-mixed classes keyed by their unmixed counterparts."""

        # Optionally register base class and data class information.
        # This can be done when creating the base inmix class, or subclasses
        # thereof.
        if base_cls is not None:
            self.__inmixbase__._inmixed_predefined_classes[base_cls] = weakref.ref(self)

        self.__intomixclass__: Optional[type] = None
        if data_cls is not None:
            self.__intomixclass__ = data_cls  # (injected into namespace in __new__)
            self.__inmixbase__._inmixed_generated_classes[data_cls] = weakref.ref(self)

            if self.__doc__ is None:
                self.__doc__ = self._inmix_make__doc__(data_cls)

        super().__init__(name, bases, namespace, **kwargs)

    def __call__(self: Self, *args: Any, **kwargs: Any) -> Union[Self, Type[Self]]:
        # Make an instance of the class. Before `cls.__new__` is called this
        # decides whether 1 of 2 things happens:
        # 1) Normal class instantiation proceeds
        # 2) "Factory" mode: defining either a new wrapped class or a making
        #    an instance of a wrapped class.

        if self is self.__inmixbase__:  # "Factory" mode
            # Defining a new class -- a sub-class of the argument, with an
            # injected inheritance.
            if not kwargs and len(args) == 1 and isinstance(args[0], type):
                return self._inmix_make_class(args[0])
            # Making an instance of a dependency-injected class.
            else:
                return self._inmix_make_instance(*args, **kwargs)

        # Normal class instantiation
        return super().__call__(*args, **kwargs)

    @property
    def _inmixed_predefined_classes(self) -> weakref.WeakKeyDictionary:
        return self.__inmixbase__.__inmixed_predefined_classes

    @property
    def _inmixed_generated_classes(self) -> weakref.WeakKeyDictionary:
        return self.__inmixbase__.__inmixed_generated_classes

    @property
    def inmix_default_cls(self: Self) -> Optional[Self]:
        """The default wrapper base class.

        Returns
        -------
        type
            The contents of ``__inmixed_default_defined_class``, possibly replacing
            a string placeholder with the actual class.

        Notes
        -----
        In the reference implementation (`astropy.utils.masked.Masked`) this is
        the subclass of `astropy.utils.masked.Masked` for masking
        `numpy.ndarray` and array-like objects.
        """
        wcls = self.__inmixed_default_defined_class  # Stored on each cls, for class overrides

        # Delayed evaluation. Subclasses are defined after the main class, so
        # `__inmixed_default_defined_class` can't be a class kwarg on the main
        # class. If the string name was given, we find the correct base class
        # and store it for reuse.
        if isinstance(wcls, str):
            for kls_ref in self.__inmixbase__._inmixed_predefined_classes.values():
                kls = kls_ref()
                if kls is not None and kls.__qualname__ == wcls:
                    wcls = kls
                    break
            else:
                msg = (f"{wcls} is not the ``__qualname__`` of any "
                       f"inmix base classes known to {self.__inmixbase__}")
                raise TypeError(msg)
            self.__inmixed_default_defined_class = wcls

        return wcls

    def _inmix_prepare_type(self, data_cls: type, base_cls: type) -> Tuple[str, Tuple[type, ...], Dict[str, Any]]:
        """Prepare a subclass of the wrapper class.

        Generally in the MRO the wrapper class is after the data class.

        Parameters
        ----------
        data_cls : type
            The class of the data, e.g. `ndarray`.
        base_cls : type
            Wrapper class, a subclass of the base wrapper class.

        Returns
        -------
        name : str
        bases : tuple[type, ...]
        """
        name: str = data_cls.__name__ + self.__name__
        bases: Tuple[type, ...] = (data_cls, base_cls)
        return name, bases, {}

    def _inmix_make__doc__(self, data_cls: type) -> Optional[str]:
        """
        Make the docstring for the wrapped class.
        The default is to use the docstring of the data class.
        """
        return data_cls.__doc__

    def _inmix_make_class(self: Self, data_cls: type) -> Self:
        """Get the wrapper for a given data class.

        Parameters
        ----------
        data_cls : type
            The given data class.

        Returns
        -------
        type
            The wrapper class.
        """
        # If the class is already correct, return unchanged.
        if issubclass(data_cls, self.__inmixbase__):
            return data_cls

        # Start by looking through the generated classes
        inmixed_ref = self.__inmixbase__._inmixed_generated_classes.get(data_cls)
        inmixed = inmixed_ref() if isinstance(inmixed_ref, weakref.ReferenceType) else None

        # If no generated class is good, a new class will have to be generated
        # (or the default class used).
        if inmixed is None:
            # Walk through MRO and find closest base data class.
            # Note: right now, will basically always be ndarray, but
            # one could imagine needing some special care for one subclass,
            # which would then get its own entry.
            for mro_item in data_cls.__mro__:
                bcls_ref = self.__inmixbase__._inmixed_predefined_classes.get(mro_item)
                base_cls = bcls_ref() if isinstance(bcls_ref, weakref.ReferenceType) else None
                if base_cls is not None:
                    break
            else:
                # If there's a default base class, hope that can handle the data,
                # otherwise hopefully this class is a good base class.
                if self.inmix_default_cls is not None:
                    return self.inmix_default_cls
                else:
                    base_cls = self

            # Create (and therefore register) new wrapper subclass for the
            # given data_cls. Generally in the MRO the wrapper class is
            # inserted below the data class, but the specifics are controlled
            # by `_inmix_prepare_type`.
            name, bases, namespace = self._inmix_prepare_type(data_cls, base_cls)
            inmixed = type(name, bases, namespace, data_cls=data_cls)

        return inmixed

    @abstractmethod
    def _inmix_make_instance(self: Self, data: Any, *args: Any, **kwargs: Any) -> Self:
        """Make an instance of the wrapper subclass for the data type.

        Parameters
        ----------
        data : Any
            The data for which to make an instance of the relevant wrapper
            subclass.
        *args, **kwargs : Any
            Arguments, positional and keyword, creating an instance of the
            wrapper subclass for the data type.

        Returns
        -------
        object

        Examples
        --------
        The simplest example is

        .. code-block::

            inmixcls = cls._inmix_make_class(type(data))
            return inmixcls(data, *args, **kwargs)
        """
        raise NotImplementedError(f"{self.__inmixbase__} must define a "
                                  "`_inmix_make_instance` classmethod.")
