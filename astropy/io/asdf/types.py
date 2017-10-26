# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import six

from asdf.asdftypes import CustomType, ExtensionTypeMeta


__all__ = ['AstropyAsdfType']


_astropy_asdf_types = set()


class AstropyTypeMeta(ExtensionTypeMeta):
    """
    Keeps track of `AstropyType` subclasses that are created so that they can
    be stored automatically by AstropyExtension
    """
    def __new__(mcls, name, bases, attrs):
        cls = super(AstropyTypeMeta, mcls).__new__(mcls, name, bases, attrs)
        # Classes using this metaclass are automatically added to the list of
        # astropy extensions
        if cls.organization == 'stsci.edu' and cls.standard == 'asdf':
            _astropy_asdf_types.add(cls)

        return cls


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
