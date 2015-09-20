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

    # uncertainty and its setter are implemented as concrete to enforce the
    # logic in the uncertainty setter. For a long discussion of the problems
    # with trying to implement them as abstract (particularly the setter but
    # not the getter), see http://bugs.python.org/issue11610
    #
    # In python >= 3.3 it would be easy to decorate one of these (setter or
    # getter) as abstract but not the other.
    @property
    def uncertainty(self):
        """
        Uncertainty in the data.

        Uncertainty must have an attribute ``uncertainty_type`` that is
        a string.
        """
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            if (not hasattr(value, 'uncertainty_type') or
                    not isinstance(value.uncertainty_type, six.string_types)):

                raise TypeError('Uncertainty must have attribute '
                                'uncertainty_type whose type is string.')
        self._uncertainty = value
        #Try setting the parent for the uncertainty. This is currently an
        #enforced requirement for `StdDevUncertainty` for error propagation.
        try:
            self.uncertainty.parent_nddata = self
        except:
            pass
