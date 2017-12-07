# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import six

from asdf.asdftypes import CustomType, ExtensionTypeMeta


__all__ = ['AstropyType', 'AstropyAsdfType']


_astropy_types = set()
_astropy_asdf_types = set()


class AstropyTypeMeta(ExtensionTypeMeta):
    """
    Keeps track of `AstropyType` subclasses that are created so that they can
    be stored automatically by astropy extensions for ASDF.
    """
    def __new__(mcls, name, bases, attrs):
        cls = super(AstropyTypeMeta, mcls).__new__(mcls, name, bases, attrs)
        # Classes using this metaclass are automatically added to the list of
        # astropy extensions
        if cls.organization == 'astropy.org' and cls.standard == 'astropy':
            _astropy_types.add(cls)
        elif cls.organization == 'stsci.edu' and cls.standard == 'asdf':
            _astropy_asdf_types.add(cls)

        return cls


@six.add_metaclass(AstropyTypeMeta)
class AstropyType(CustomType):
    """
    This class represents types that have schemas and tags that are defined by
    Astropy.

    IMPORTANT: This parent class should **not** be used for types that have
    schemas that are defined by the ASDF standard.
    """
    organization = 'astropy.org'
    standard = 'astropy'


@six.add_metaclass(AstropyTypeMeta)
class AstropyAsdfType(CustomType):
    """
    This class represents types that have schemas that are defined in the ASDF
    standard, but have tags that are implemented within astropy.

    IMPORTANT: This parent class should **not** be used for types that also
    have schemas that are defined by astropy.
    """
    organization = 'stsci.edu'
    standard = 'asdf'
