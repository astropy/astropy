# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module containing specialized collection classes.
"""

from abc import ABCMeta, abstractmethod
from typing import Tuple


class HomogeneousList(list):
    """
    A subclass of list that contains only elements of a given type or
    types.  If an item that is not of the specified type is added to
    the list, a `TypeError` is raised.
    """
    def __init__(self, types, values=[]):
        """
        Parameters
        ----------
        types : sequence of types
            The types to accept.

        values : sequence, optional
            An initial set of values.
        """
        self._types = types
        super().__init__()
        self.extend(values)

    def _assert(self, x):
        if not isinstance(x, self._types):
            raise TypeError(
                f"homogeneous list must contain only objects of type '{self._types}'")

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __setitem__(self, idx, value):
        if isinstance(idx, slice):
            value = list(value)
            for item in value:
                self._assert(item)
        else:
            self._assert(value)
        return super().__setitem__(idx, value)

    def append(self, x):
        self._assert(x)
        return super().append(x)

    def insert(self, i, x):
        self._assert(x)
        return super().insert(i, x)

    def extend(self, x):
        for item in x:
            self._assert(item)
            super().append(item)


class ClassWrapperMeta(ABCMeta):
    """Meta-class for creation of wrapper classes around different data classes.

    Warnings
    --------
    This is an experimental feature and subject to change. In particular,
    the private methods are not stable.

    Examples
    --------
    We want to make a class that wraps stuff.

        >>> class ExampleWrapper(metaclass=ClassWrapperMeta):
        ...     @classmethod
        ...     def _get_wrapper_subclass_instance(cls, data, *args, **kwargs):
        ...         wrapper_cls = cls._get_wrapped_subclass(data.__class__)
        ...         return wrapper_cls(data, *args, **kwargs)
    """

    def __new__(metacls, name, bases, namespace, base_cls=None, data_cls=None,
                **kwargs):
        # Make a new class. `__new__` is required to prevent class-level
        # kwargs from going to the superclass (ABCMeta), which does not take
        # any kwargs. Also, `__init_subclass__`, which may be defined on the
        # classes `ClassWrapperMeta` is trying to make, will not have access
        # to certain class-level attributes that are made during the class
        # creation process.

        # ABCMeta takes no kwargs. This is still passed to ClassWrapperMeta.__init__.
        kwargs.pop("default_wrapped_class", None)

        # inject data_cls into namespace before `__new__` so
        # `__init_subclass__` can use it.
        if data_cls is not None:
            namespace["_wrapped_data_cls"] = data_cls

        return super().__new__(metacls, name, bases, namespace, **kwargs)

    def __init__(cls, name, bases, namespace, base_cls=None, data_cls=None,
                 **kwargs):

        # check that this `cls` is the top-most in the MRO to be a `ClassWrapperMeta` type
        # TODO! allow for subclasses of a ClassWrapperMeta class to be considered
        # a baseclass.
        isbaseclass = all([not issubclass(base.__class__, ClassWrapperMeta) for base in bases])
        if isbaseclass:
            cls._wrapper_class_ = cls
            """Base class."""

            cls._wrapper_default_base_class_ = kwargs.pop("default_wrapped_class", None)
            """
            Fallback wrapper class when there is not a good wrapper base class.
            See ``_wrapper_base_classes``.
            """

            cls._wrapper_base_classes = {}
            """Explicitly defined wrapper classes keyed by their unwrapped counterparts.

            For subclasses of these unwrapped classes, wrapped counterparts can
            be generated.
            """

            cls._wrapper_generated_subclasses = {}
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

    def __call__(cls, *args, **kwargs):
        # Make an instance of the class. Before `cls.__new__` is called this
        # decides whether 1 of 2 things happens:
        # 1) Normal class instantiation proceeds
        # 2) "Factory" mode: defining either a new wrapped class or a making
        #    an instance of a wrapped class.

        if cls is cls._wrapper_class_:  # "Factory" mode
            # Defining a new wrapped class.
            if not kwargs and len(args) == 1 and isinstance(args[0], type):
                return cls._get_wrapped_subclass(args[0])
            # Making an instance of a wrapped class
            else:
                return cls._get_wrapper_subclass_instance(*args, **kwargs)

        # Normal class instantiation
        return super().__call__(*args, **kwargs)

    @abstractmethod
    def _get_wrapper_subclass_instance(cls, data, *args, **kwargs) -> object:
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

    def _prepare_wrapper_subclass(cls, data_cls: type, base_cls: type) -> Tuple[str, Tuple[type, ...]]:
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
        name = data_cls.__name__ + cls.__name__
        bases = (data_cls, base_cls)
        return name, bases

    @property
    def wrapper_default_base_class(cls) -> type:
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

    def _make_wrapped__doc__(cls, data_cls: type) -> str:
        """
        Make the docstring for the wrapped class.
        The default is to use the docstring of the data class.
        """
        return data_cls.__doc__

    def _get_wrapped_subclass(cls, data_cls: type) -> type:
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
