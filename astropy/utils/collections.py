# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module containing specialized collection classes.
"""

import abc


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


class ClassWrapperMeta(abc.ABCMeta):

    def __new__(metacls, name, bases, namespace, base_cls=None, data_cls=None,
                **kwargs):
        # inject data_cls into namespace before __new__ so __init_subclass__
        # can use it.
        if data_cls is not None:
            namespace["_wrapped_data_cls"] = data_cls

        # ABCMeta takes no kwargs. This is still passed to __init__
        kwargs.pop("default_wrapped_class", None)

        return super().__new__(metacls, name, bases, namespace, **kwargs)

    def __init__(cls, name, bases, namespace, base_cls=None, data_cls=None,
                 **kwargs):

        # check that this `cls` is the top-most in the MRO to be a `ClassWrapperMeta` type
        isbaseclass = all([not issubclass(base.__class__, ClassWrapperMeta) for base in bases])
        if isbaseclass:
            cls._wrapper_class_ = cls
            """Base class."""

            cls._wrapper_default_base_class_ = kwargs.pop("default_wrapped_class", None)
            """Backup explicitly defined wrapper class."""

            cls._wrapper_base_classes = {}
            """Explicitly defined wrapper classes keyed by their unwrapped counterparts.

            For subclasses of these unwrapped classes, wrapped counterparts can be generated.
            """

            cls._wrapper_generated_subclasses = {}
            """Wrapped classes keyed by their unwrapped data counterparts.
            """

        if base_cls is not None:
            cls._wrapper_class_._wrapper_base_classes[base_cls] = cls

        if data_cls is not None:
            cls._wrapped_data_cls = data_cls  # injected into namespace in __new__
            cls._wrapper_class_._wrapper_generated_subclasses[data_cls] = cls

            if cls.__doc__ is None:
                cls.__doc__ = cls._make_wrapped__doc__(data_cls)

        super().__init__(name, bases, namespace, **kwargs)

    def __call__(cls, *args, **kwargs):
        if cls is cls._wrapper_class_:
            # Initializing with Masked itself means we're in "factory mode".
            if not kwargs and len(args) == 1 and isinstance(args[0], type):
                # Create a new masked class.
                return cls._get_wrapped_subclass(args[0])
            else:
                return cls._get_wrapper_subclass_instance(*args, **kwargs)

        # Otherwise we're a subclass and should just pass information on.
        return super().__call__(*args, **kwargs)

    @abc.abstractmethod
    def _get_wrapper_subclass_instance(cls, data):
        pass

    @abc.abstractmethod
    def _make_wrapper_subclass(cls, data_cls, base_cls):
        pass

    @property
    def wrapper_default_base_class(cls):
        wcls = cls._wrapper_default_base_class_  # Stored on each cls, for class overrides

        # Delayed evaluation. Subclasses are defined after the main class, so
        # can't be a class kwarg on the main class. If the string name was
        # given, we find the correct base class and store it for reuse.
        if isinstance(wcls, str):
            bclss = {kls.__qualname__: kls for kls
                     in cls._wrapper_class_._wrapper_base_classes.values()}
            wcls = bclss[wcls]
            cls._wrapper_default_base_class_ = wcls

        return wcls

    def _make_wrapped__doc__(cls, data_cls):
        return data_cls.__doc__

    def _get_wrapped_subclass(cls, data_cls):
        """Get the wrapper for a given data class."""
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
            # given data_cls. In the MRO the wrapper class is inserted below
            # the data class.
            wrapper_cls = cls._make_wrapper_subclass(data_cls, base_cls)

        return wrapper_cls
