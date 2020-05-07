# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from asdf.types import CustomType, ExtensionTypeMeta

__all__ = ['AstropyType', 'AstropyAsdfType']

# Names of AstropyType or AstropyAsdfType subclasses that are base classes
# and aren't used directly for serialization.
_TYPE_BASE_CLASS_NAMES = {'PolynomialTypeBase'}

_astropy_types = set()
_astropy_asdf_types = set()


class AstropyTypeMeta(ExtensionTypeMeta):
    """
    Keeps track of `AstropyType` subclasses that are created so that they can
    be stored automatically by astropy extensions for ASDF.
    """
    def __new__(mcls, name, bases, attrs):
        cls = super().__new__(mcls, name, bases, attrs)
        # Classes using this metaclass are automatically added to the list of
        # astropy extensions
        if cls.__name__ not in _TYPE_BASE_CLASS_NAMES:
            if cls.organization == 'astropy.org' and cls.standard == 'astropy':
                _astropy_types.add(cls)
            elif cls.organization == 'stsci.edu' and cls.standard == 'asdf':
                _astropy_asdf_types.add(cls)

        return cls


class AstropyType(CustomType, metaclass=AstropyTypeMeta):
    """
    This class represents types that have schemas and tags that are defined by
    Astropy.

    IMPORTANT: This parent class should **not** be used for types that have
    schemas that are defined by the ASDF standard.
    """
    organization = 'astropy.org'
    standard = 'astropy'


class AstropyAsdfType(CustomType, metaclass=AstropyTypeMeta):
    """
    This class represents types that have schemas that are defined in the ASDF
    standard, but have tags that are implemented within astropy.

    IMPORTANT: This parent class should **not** be used for types that also
    have schemas that are defined by astropy.
    """
    organization = 'stsci.edu'
    standard = 'asdf'
