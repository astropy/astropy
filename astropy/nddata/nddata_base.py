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
    Base metaclass that defines the interface for N-dimensional datasets with
    associated meta informations used in `astropy`.

    NDDataBase makes no attempt of restricting any attribute or defines how
    these are stored. Therefore all properties and methods have to be
    overridden in subclasses. See `~astropy.nddata.NDData` for a subclass that
    defines this interface on `~numpy.ndarray`-like based ``data``.
    """

    @abstractmethod
    def __init__(self):
        pass

    @abstractproperty
    def data(self):
        """
        The stored dataset.
        """
        pass

    @abstractproperty
    def mask(self):
        """
        Mask for the dataset.

        Masks should follow the `numpy` convention that valid data points are
        marked by ``False`` and invalid ones with ``True``.
        """
        return None

    @abstractproperty
    def unit(self):
        """
        Unit for the dataset.
        """
        return None

    @abstractproperty
    def wcs(self):
        """
        A world coordinate system (WCS) for the dataset.
        """
        return None

    @abstractproperty
    def meta(self):
        """
        Meta information about the dataset.

        Should be `dict`-like.
        """
        return None

    @abstractproperty
    def uncertainty(self):
        """
        Uncertainty in the dataset.

        Should have an attribute ``uncertainty_type`` that defines what kind of
        uncertainty is stored, such as ``'std'`` for standard deviation or
        ``'var'`` for variance.
        """
        return None
