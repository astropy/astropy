# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDDataBase class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from abc import ABCMeta, abstractproperty, abstractmethod

from ..extern import six

__all__ = ['NDDataBase']


@six.add_metaclass(ABCMeta)
class NDDataBase(object):
    """
    Base metaclass that defines the interface for NDData

    Classes that wish to use this interface without inheriting from
    `~astropy.nddata.NDData` should subclass ``NDDataBase`` instead.

    All properties and methods except uncertainty must be override by derived
    classes.
    """

    @abstractmethod
    def __init__(self):
        self._uncertainty = None

    @abstractproperty
    def data(self):
        """
        The data; should be capable of behaving like a numpy array, though it
        need not actually be a numpy array.
        """
        pass

    @abstractproperty
    def mask(self):
        """
        Mask for the data, following the numpy convention that ``True`` means
        the data should not be used.
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
        Metadata, if any, must be dict-like.
        """
        return None

    @abstractproperty
    def uncertainty(self):
        """
        Uncertainty in the data.
        """
        return None
