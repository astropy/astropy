# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Astropy's ABCs."""


from abc import ABCMeta, abstractmethod
from typing import Tuple, Type, Dict, Any, Optional, TypeVar, Union


__all__ = []  # InheritanceInMixMeta is in development


CWT = TypeVar("CWT", bound="InheritanceInMixMeta")


class InheritanceInMixMeta(ABCMeta):
    """Meta-class for generating classes of a data class.

    Parameters
    ----------
    base_cls : type or None, optional

    data_cls : type or None, optional
    **kwargs : Any
        A valid option is ``default_wrapped_class``. This should only be given
        to the base wrapper class. In Examples this is ``ExampleWrapper``.

    Warnings
    --------
    This is an experimental feature and subject to change. In particular,
    the private methods are not stable.

    Examples
    --------
    For real examples implemented in the code, see
    `astropy.utils.masked.Masked` or `astropy.uncertainty.Distribution`.
    Here we will work through a toy problem.

    Let's say we like how `numpy.ndarray` can cast elements of an array to
    different datatypes and we want to add this functionality to some builtin
    types.
    This starts by defining the inter-mix base class, using `InheritanceInMixMeta`
    as a metaclass. The ``_get_wrapper_subclass_instance`` defines how new
    instances are constructed. In this case it's just calling the same args
    and kwargs. We also add a method ``astype`` for casting each element
    of self to a float.

        >>> class SupportsAsType(metaclass=InheritanceInMixMeta):
        ...     @classmethod
        ...     def _get_wrapper_subclass_instance(cls, data, *args, **kwargs):
        ...         wrapper_cls = cls._get_wrapped_subclass(data.__class__)
        ...         return wrapper_cls(data, *args, **kwargs)
        ...     def astype(self, type):
        ...         return self.__class__([type(x) for x in self])

    Now implementations must be made by subclassing ``SupportsAsType``
    and also inheriting from the desired type.

        >>> class SupportsAsTypeList(SupportsAsType, list, data_cls=list):
        ...     pass

    The kwarg ``data_cls`` registered `list`, so the following is possible.

        >>> cls = SupportsAsType(list)
        >>> cls
        <class 'astropy.utils.metaclasses.SupportsAsTypeList'>

    ``SupportsAsTypeList`` can make list-like objects...

        >>> x = cls([1, 2, 3])
        >>> x
        [1, 2, 3]

    with the added method ``astype``.

        >>> x.astype(float)
        [1.0, 2.0, 3.0]
    """

    def __new__(metacls: Type[CWT], name: str, bases: Tuple[type, ...],
                namespace: Dict[str, Any], base_cls: Optional[type]=None,
                data_cls: Optional[type]=None, **kwargs: Any) -> CWT:
        # Make a new class. `__new__` is required to prevent class-level
        # kwargs from going to the superclass (ABCMeta), which does not take
        # any kwargs. Also, `__init_subclass__`, which may be defined on the
        # classes `InheritanceInMixMeta` is trying to make, will not have access
        # to certain class-level attributes that are made during the class
        # creation process.

        # ABCMeta takes no kwargs. This is still passed to InheritanceInMixMeta.__init__.
        kwargs.pop("default_wrapped_class", None)

        # inject data_cls into namespace before `__new__` so
        # `__init_subclass__` can use it.
        if data_cls is not None:
            namespace["_wrapped_data_cls"] = data_cls

        return super().__new__(metacls, name, bases, namespace, **kwargs)

    def __init__(cls: CWT, name: str, bases: Tuple[type, ...],
                 namespace: Dict[str, Any], base_cls: Optional[type]=None,
                 data_cls: Optional[type]=None, **kwargs: Any) -> None:
        # check that this `cls` is the top-most in the MRO to be a
        # `InheritanceInMixMeta` type TODO! allow for subclasses of a
        # InheritanceInMixMeta class to be considered a baseclass.
        isbaseclass = all([not issubclass(base.__class__, InheritanceInMixMeta) for base in bases])
        if isbaseclass:
            cls._wrapper_class_: CWT = cls
            """Base class."""

            cls._wrapper_default_base_class_: Union[str, CWT, None]
            cls._wrapper_default_base_class_ = kwargs.pop("default_wrapped_class", None)
            """
            Fallback wrapper class when there is not a good wrapper base class.
            See ``_wrapper_base_classes``.
            """

            cls._wrapper_base_classes: Dict[type, CWT] = {}
            """Explicitly defined wrapper classes keyed by their unwrapped counterparts.

            For subclasses of these unwrapped classes, wrapped counterparts can
            be generated.
            """

            cls._wrapper_generated_subclasses: Dict[type, CWT] = {}
            """Wrapped classes keyed by their unwrapped data counterparts."""

        # Optionally register base class and data class information.
        # This can be done when creating the base wrapper class, or subclasses
        # thereof.
        if base_cls is not None:
            cls._wrapper_class_._wrapper_base_classes[base_cls] = cls

        if data_cls is not None:
            cls._wrapped_data_cls = data_cls  # (injected into namespace in __new__)
            cls._wrapper_class_._wrapper_generated_subclasses[data_cls] = cls

            if cls.__doc__ is None:
                cls.__doc__ = cls._make_wrapped__doc__(data_cls)

        super().__init__(name, bases, namespace, **kwargs)

    def __call__(cls: CWT, *args: Any, **kwargs: Any) -> Union[CWT, Type[CWT]]:
        # Make an instance of the class. Before `cls.__new__` is called this
        # decides whether 1 of 2 things happens:
        # 1) Normal class instantiation proceeds
        # 2) "Factory" mode: defining either a new wrapped class or a making
        #    an instance of a wrapped class.

        if cls is cls._wrapper_class_:  # "Factory" mode
            # Defining a new class -- a sub-class of the argument, with an
            # injected inheritance.
            if not kwargs and len(args) == 1 and isinstance(args[0], type):
                return cls._get_wrapped_subclass(args[0])
            # Making an instance of a dependency-injected class.
            else:
                return cls._get_wrapper_subclass_instance(*args, **kwargs)

        # Normal class instantiation
        return super().__call__(*args, **kwargs)

    @abstractmethod
    def _get_wrapper_subclass_instance(cls: CWT, data: Any, *args: Any, **kwargs: Any) -> CWT:
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
        """
        raise NotImplementedError(f"{cls._wrapper_class_} must define a "
                                  "`_get_wrapper_subclass_instance` classmethod.")

    def _prepare_wrapper_subclass(cls: CWT, data_cls: type, base_cls: type
                                 ) -> Tuple[str, Tuple[type, ...]]:
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
        name: str = data_cls.__name__ + cls.__name__
        bases: Tuple[type, ...] = (data_cls, base_cls)
        return name, bases

    @property
    def wrapper_default_base_class(cls: CWT) -> Optional[CWT]:
        """The default wrapper base class.

        Returns
        -------
        type
            The contents of ``_wrapper_default_base_class_``, possibly replacing
            a string placeholder with the actual class.

        Notes
        -----
        In the reference implementation (`astropy.utils.masked.Masked`) this is
        the subclass of `astropy.utils.masked.Masked` for masking
        `numpy.ndarray` and array-like objects.
        """
        wcls = cls._wrapper_default_base_class_  # Stored on each cls, for class overrides

        # Delayed evaluation. Subclasses are defined after the main class, so
        # `_wrapper_default_base_class_` can't be a class kwarg on the main
        # class. If the string name was given, we find the correct base class
        # and store it for reuse.
        if isinstance(wcls, str):
            bclss = {kls.__qualname__: kls for kls
                     in cls._wrapper_class_._wrapper_base_classes.values()}
            wcls = bclss[wcls]
            cls._wrapper_default_base_class_ = wcls

        return wcls

    def _make_wrapped__doc__(cls: CWT, data_cls: type) -> str:
        """
        Make the docstring for the wrapped class.
        The default is to use the docstring of the data class.
        """
        return data_cls.__doc__

    def _get_wrapped_subclass(cls: CWT, data_cls: type) -> CWT:
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
        if issubclass(data_cls, cls._wrapper_class_):
            return data_cls

        wrapper_cls = cls._wrapper_class_._wrapper_generated_subclasses.get(data_cls)
        if wrapper_cls is None:
            # Walk through MRO and find closest base data class.
            # Note: right now, will basically always be ndarray, but
            # one could imagine needing some special care for one subclass,
            # which would then get its own entry.
            for mro_item in data_cls.__mro__:
                base_cls = cls._wrapper_class_._wrapper_base_classes.get(mro_item)
                if base_cls is not None:
                    break
            else:
                # Hope the default base class can handle this
                return cls.wrapper_default_base_class

            # Create (and therefore register) new wrapper subclass for the
            # given data_cls. Generally in the MRO the wrapper class is
            # inserted below the data class, but the specifics are controlled
            # by `_prepare_wrapper_subclass`.
            name, bases = cls._prepare_wrapper_subclass(data_cls, base_cls)
            wrapper_cls = type(name, bases, {}, data_cls=data_cls)

        return wrapper_cls
