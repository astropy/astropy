# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, print_function, unicode_literals

import abc

import numpy as np

from ..utils.compat import ignored
from ..units import Quantity
from ..extern import six

__all__ = ['MissingDataAssociationException',
           'IncompatibleUncertaintiesException', 'NDUncertainty',
           'StdDevUncertainty']


class IncompatibleUncertaintiesException(Exception):
    """
    This exception should be used to indicate cases in which uncertainties
    with two different classes can not be propagated.
    """


class MissingDataAssociationException(Exception):
    """
    This exception should be used to indicate that an uncertainty instance has
    not been associated with a parent `~astropy.nddata.NDData` object.
    """


class NDUncertainty(object):
    '''
    This is the base class for uncertainty classes used with NDData. It is
    implemented as an abstract class and should never be directly
    instantiated.

    Classes inheriting from NDData should overload the ``propagate_*``
    methods, keeping the call signature the same. The propagate methods can
    assume that a `parent_nddata` attribute is present which links to the parent_nddata
    dataset, and take an `~astropy.nddata.NDData` instance as the positional
    argument, *not* an `~astropy.nddata.NDUncertainty` instance, because the
    `~astropy.nddata.NDData` instance can be used to access both the data and
    the uncertainties (some propagations require the data values).
    '''

    __metaclass__ = abc.ABCMeta

    # Indicates whether the class supports the propagation of correlated
    # uncertainties
    supports_correlated = False

    @property
    def parent_nddata(self):
        if self._parent_nddata is None:
            raise MissingDataAssociationException("uncertainty is not associated with an NDData object")
        else:
            return self._parent_nddata

    @parent_nddata.setter
    def parent_nddata(self, value):
        self._parent_nddata = value

    @property
    def uncertainty_type(self):
        return self._uncertainty_type

    @uncertainty_type.setter
    def uncertainty_type(self, value):
        if not isinstance(value, six.string_types):
            raise ValueError('uncertainty_type must be a string.')
        self._uncertainty_type = value

    @abc.abstractmethod
    def propagate_add(self, other_nddata, result_data):
        '''
        Propagate uncertainties for addition.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to add the uncertainties
        '''

    @abc.abstractmethod
    def propagate_subtract(self, other_nddata, result_data):
        '''
        Propagate uncertainties for subtraction.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to add the uncertainties
        '''

    @abc.abstractmethod
    def propagate_multiply(self, other_nddata, result_data):
        '''
        Propagate uncertainties for multiplication.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty
        '''

    @abc.abstractmethod
    def propagate_divide(self, other_nddata, result_data):
        '''
        Propagate uncertainties for division.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty
        '''


class StdDevUncertainty(NDUncertainty):
    '''
    A class for standard deviation uncertainties
    '''

    support_correlated = False

    def __init__(self, array=None, copy=True):
        self._unit = None
        self.uncertainty_type = 'std'
        if array is None:
            self.array = None
        elif isinstance(array, StdDevUncertainty):
            self.array = np.array(array.array, copy=copy, subok=True)
        elif isinstance(array, Quantity):
            self.array = np.array(array.value, copy=copy, subok=True)
            self._unit = array.unit
        else:
            self.array = np.array(array, copy=copy, subok=True)

    @property
    def parent_nddata(self):
        try:
            if self._parent_nddata is None:
                raise MissingDataAssociationException("Uncertainty is not associated with an NDData object")
            else:
                return self._parent_nddata
        except AttributeError:
            raise MissingDataAssociationException("Uncertainty is not associated with an NDData object")

    @parent_nddata.setter
    def parent_nddata(self, value):
        if self.array is None or value is None:
            self._parent_nddata = value
        else:
            if value.shape != self.array.shape:
                raise ValueError("parent shape does not match array data shape")
        self._parent_nddata = value

    @property
    def array(self):
        return self._array

    @array.setter
    def array(self, value):
        if value is not None:
            with ignored(MissingDataAssociationException):
                if value.shape != self.parent_nddata.shape:
                    raise ValueError(
                            "array shape does not match parent data shape")
        self._array = value

    def propagate_add(self, other_nddata, result_data):
        '''
        Propagate uncertainties for addition.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the uncertainties
        '''

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = np.sqrt(self.array ** 2 + other_nddata.uncertainty.array ** 2)

        return result_uncertainty

    def __getitem__(self, item):
        '''
            Slice the standard deviation array using standard numpy slicing
        '''

        new_array = self.array[item]
        return self.__class__(new_array, copy=False)

    def propagate_subtract(self, other_nddata, result_data):
        '''
        Propagate uncertainties for subtraction.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the uncertainties
        '''

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = np.sqrt(self.array ** 2 + other_nddata.uncertainty.array ** 2)

        return result_uncertainty

    def propagate_multiply(self, other_nddata, result_data):
        '''
        Propagate uncertainties for mutliplication.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the uncertainties
        '''

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = (np.sqrt((self.array / self.parent_nddata.data) ** 2
                                            + (other_nddata.uncertainty.array / other_nddata.data) ** 2)
                                    * result_data)

        return result_uncertainty

    def propagate_divide(self, other_nddata, result_data):
        '''
        Propagate uncertainties for division.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the uncertainties
        '''

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = (np.sqrt((self.array / self.parent_nddata.data) ** 2
                                            + (other_nddata.uncertainty.array / other_nddata.data) ** 2)
                                    * result_data)

        return result_uncertainty
