# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""A module containing utilities for working with dataclasses."""

__all__ = ["hasdocstringfield"]

from dataclasses import dataclass

from astropy.utils.compat.misc import PYTHON_LT_3_10

if PYTHON_LT_3_10:
    _dataclass_kwargs = {}
else:
    _dataclass_kwargs = {"slots": True}


@dataclass(frozen=True, **_dataclass_kwargs)
class DocstringConnectedProperty:
    """Descriptor returning different docstrings based on access from class or instance.

    Users are unlikely to need to use this class directly, but it is used by
    :func:`~astropy.utils.dataclasses.hasdocstringfield` to allow for a slotted
    dataclass to have a docstring that is different when accessed from the class or an
    instance.
    """

    name: str
    cls__doc__: str

    def __get__(self, instance, owner):
        return getattr(instance, self.name) if instance is not None else self.cls__doc__


def hasdocstringfield(name: str):
    """Decorator to add a docstring to a field in a slotted dataclass.

    This decorator is intended to be used with slotted dataclasses, which do not
    allow for ``__doc__`` to be a field.  This decorator allows for a different
    field, specified by ``name`` (e.g. "doc"), to also set the docstring of the
    slotted class.

    Parameters
    ----------
    name : str
        The name of the field that should be used to set the docstring of the
        slotted class.

    Examples
    --------
    >>> from dataclasses import dataclass
    >>> from astropy.utils.compat.misc import PYTHON_LT_3_10
    >>> from astropy.utils.dataclasses import hasdocstringfield

    >>> @dataclass(**({"slots": True} if not PYTHON_LT_3_10 else {}))
    ... @hasdocstringfield('doc')
    ... class Foo:
    ...     '''This is the docstring of the class Foo.'''
    ...     doc: str

    >>> Foo.__doc__
    'This is the docstring of the class Foo.'

    >>> foo = Foo('This is the docstring of an instance of Foo.')
    >>> foo.doc
    'This is the docstring of an instance of Foo.'

    >>> foo.__doc__
    'This is the docstring of an instance of Foo.'

    >>> foo.__doc__ is foo.doc
    True

    Notes
    -----
    This decorator must be used *before* the ``dataclass`` decorator, as in the
    example above.  This is because the ``dataclass`` decorator will make the class
    slotted, which will prevent the ``__doc__`` attribute from being set on the
    class and its instances.
    """

    def decorator(cls):
        cls.__doc__ = DocstringConnectedProperty(name, cls.__doc__)
        return cls

    return decorator
