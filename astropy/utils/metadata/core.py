# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Classes for handling metadata."""

__all__ = ["MetaAttribute", "MetaData"]

import inspect
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import is_dataclass


class MetaData:
    """
    A descriptor for classes that have a ``meta`` property.

    This can be set to any valid :class:`~collections.abc.Mapping`.

    Parameters
    ----------
    doc : `str`, optional
        Documentation for the attribute of the class.
        Default is ``""``.

        .. versionadded:: 1.2

    copy : `bool`, optional
        If ``True`` the value is deepcopied before setting, otherwise it
        is saved as reference.
        Default is ``True``.

        .. versionadded:: 1.2

    default_factory : Callable[[], Mapping], optional keyword-only
        The factory to use to create the default value of the ``meta``
        attribute.  This must be a callable that returns a `Mapping` object.
        Default is `OrderedDict`, creating an empty `OrderedDict`.

        .. versionadded:: 6.0

    Examples
    --------
    ``MetaData`` can be used as a descriptor to define a ``meta`` attribute`.

        >>> class Foo:
        ...     meta = MetaData()
        ...     def __init__(self, meta=None):
        ...         self.meta = meta

    ``Foo`` can be instantiated with a ``meta`` argument.

        >>> foo = Foo(meta={'a': 1, 'b': 2})
        >>> foo.meta
        {'a': 1, 'b': 2}

    The default value of ``meta`` is an empty :class:`~collections.OrderedDict`.
    This can be set by passing ``None`` to the ``meta`` argument.

        >>> foo = Foo()
        >>> foo.meta
        OrderedDict()

    If an :class:`~collections.OrderedDict` is not a good default metadata type then
    the ``default_factory`` keyword can be used to set the default to a different
    `Mapping` type, when the class is defined.'

        >>> class Bar:
        ...     meta = MetaData(default_factory=dict)
        ...     def __init__(self, meta=None):
        ...         self.meta = meta

        >>> Bar().meta
        {}

    When accessed from the class ``.meta`` returns `None` since metadata is
    on the class' instances, not the class itself.

        >>> print(Foo.meta)
        None
    """

    def __init__(self, doc="", copy=True, *, default_factory=OrderedDict):
        self.__doc__ = doc
        self.copy = copy
        self._default_factory = default_factory

    @property
    def default_factory(self):
        return self._default_factory

    def __get__(self, instance, owner):
        # class attribute access. Often, descriptors just return `self`, but if the
        # owning class is a `dataclass`, the expectation is that the default is
        # returned. In our case, this is None, triggering the creation of a dict-like in
        # `__set__`.
        if instance is None:
            return None
        # instance attribute access
        if not hasattr(instance, "_meta"):
            self.__set__(instance, None)
        return instance._meta

    def __set__(self, instance, value):
        # The 'default' value is `None`, but we want to set it to an empty `Mapping`
        # if it is `None` so that we can always assume it is a `Mapping` and not have
        # to check for `None` everywhere.
        if value is None:
            value = self.default_factory()
        # We don't want to allow setting the meta attribute to a non-dict-like object.
        # NOTE: with mypyc compilation this can be removed.
        elif not isinstance(value, Mapping):
            raise TypeError("meta attribute must be dict-like")
        # This is called when the dataclass is instantiated with a `meta` argument.
        else:
            value = deepcopy(value) if self.copy else value

        if is_dataclass(instance) and instance.__dataclass_params__.frozen:
            object.__setattr__(instance, "_meta", value)
        else:
            instance._meta = value


class MetaAttribute:
    """
    Descriptor to define custom attribute which gets stored in the object
    ``meta`` dict and can have a defined default.

    This descriptor is intended to provide a convenient way to add attributes
    to a subclass of a complex class such as ``Table`` or ``NDData``.

    This requires that the object has an attribute ``meta`` which is a
    dict-like object.  The value of the MetaAttribute will be stored in a
    new dict meta['__attributes__'] that is created when required.

    Classes that define MetaAttributes are encouraged to support initializing
    the attributes via the class ``__init__``.  For example::

        for attr in list(kwargs):
            descr = getattr(self.__class__, attr, None)
            if isinstance(descr, MetaAttribute):
                setattr(self, attr, kwargs.pop(attr))

    The name of a ``MetaAttribute`` cannot be the same as any of the following:

    - Keyword argument in the owner class ``__init__``
    - Method or attribute of the "parent class", where the parent class is
      taken to be ``owner.__mro__[1]``.

    Parameters
    ----------
    default : Any, optional
        Default value for the attribute, by default `None`.
    """

    def __init__(self, default=None):
        self.default = default

    def __get__(self, instance, owner):
        # When called without an instance, return self to allow access
        # to descriptor attributes.
        if instance is None:
            return self

        # If default is None and value has not been set already then return None
        # without doing touching meta['__attributes__'] at all. This helps e.g.
        # with the Table._hidden_columns attribute so it doesn't auto-create
        # meta['__attributes__'] always.
        if self.default is None and self.name not in instance.meta.get(
            "__attributes__", {}
        ):
            return None

        # Get the __attributes__ dict and create if not there already.
        attributes = instance.meta.setdefault("__attributes__", {})
        try:
            value = attributes[self.name]
        except KeyError:
            if self.default is not None:
                attributes[self.name] = deepcopy(self.default)
            # Return either specified default or None
            value = attributes.get(self.name)
        return value

    def __set__(self, instance, value):
        # Get the __attributes__ dict and create if not there already.
        attributes = instance.meta.setdefault("__attributes__", {})
        attributes[self.name] = value

    def __delete__(self, instance):
        # Remove this attribute from meta['__attributes__'] if it exists.
        if "__attributes__" in instance.meta:
            attrs = instance.meta["__attributes__"]
            if self.name in attrs:
                del attrs[self.name]
            # If this was the last attribute then remove the meta key as well
            if not attrs:
                del instance.meta["__attributes__"]

    def __set_name__(self, owner, name):
        params = [
            param.name
            for param in inspect.signature(owner).parameters.values()
            if param.kind
            not in (inspect.Parameter.VAR_KEYWORD, inspect.Parameter.VAR_POSITIONAL)
        ]

        # Reject names from existing params or best guess at parent class
        if name in params or hasattr(owner.__mro__[1], name):
            raise ValueError(f"{name} not allowed as {self.__class__.__name__}")

        self.name = name

    def __repr__(self):
        return f"<{self.__class__.__name__} name={self.name} default={self.default}>"
