# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDDataBase class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from abc import ABCMeta, abstractproperty

from ..extern import six


@six.add_metaclass(ABCMeta)
class NDDataBase(object):

    @abstractproperty
    def data(self):
        raise NotImplementedError("Concrete classes must implement data")

    @abstractproperty
    def mask(self):
        return None

    @abstractproperty
    def unit(self):
        return None

    @abstractproperty
    def wcs(self):
        return None

    @abstractproperty
    def meta(self):
        return None

    # uncertainty and its setter are implemented as concrete to enforce the
    # logic in the uncertainty setter. For a long discussion of the problems
    # with trying to implement them as abstract (particularly the setter), see
    # http://bugs.python.org/issue11610
    #
    # In python >= 3.3 it would be easy to decorate one of these (setter or
    # getter) as abstract but not the other.
    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            try:
                not_good_uncertainty = not isinstance(value.uncertainty_type,
                                                      six.string_types)
            except AttributeError:
                not_good_uncertainty = True
            finally:
                if not_good_uncertainty:
                    raise TypeError('Uncertainty must have attribute '
                                    'uncertainty_type whose type is string.')
        self._uncertainty = value
