# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDDataBase class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from abc import ABCMeta, abstractproperty, abstractmethod

from ..extern.six import add_metaclass

__all__ = ['NDDataBase']


@add_metaclass(ABCMeta)
class NDDataBase(object):
    """
    Base metaclass that defines the interface for `NDData`-like classes.

    Classes that wish to use this interface without inheriting the restrictions
    and internal representations of `~astropy.nddata.NDData` should subclass
    ``NDDataBase`` instead.

    TODO: Copy parameter docstrings from NDData here(?!)
    """

    @abstractmethod
    def __init__(self):
        pass

    @abstractproperty
    def data(self):
        """
        The data.
        """
        pass

    @abstractproperty
    def mask(self):
        """
        Mask for the data, if any.
        """
        return None

    @abstractproperty
    def unit(self):
        """
        Unit for the data, if any.
        """
        return None

    @abstractproperty
    def wcs(self):
        """
        WCS for the data, if any.
        """
        return None

    @abstractproperty
    def meta(self):
        """
        Meta information about the data, if any.
        """
        return None

    @abstractproperty
    def uncertainty(self):
        """
        Uncertainty in the data, if any.
        """
        return None
