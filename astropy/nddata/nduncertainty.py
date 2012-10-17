import abc

import numpy as np

__all__ = ['MissingDataAssociationException', 'IncompatibleErrorsException', 'NDError', 'StandardDeviationError']


class IncompatibleErrorsException(Exception):
    """
    This exception should be used to indicate cases in which errors with two different classes can not be propagated.
    """
    pass


class MissingDataAssociationException(Exception):
    """
    This exception should be used to indicate that an error instance has not been associated with a parent `~astropy.nddata.nddata.NDData` object.
    """
    pass


class NDError(object):
    '''
    This is the base class for error classes used with NDData. It is
    implemented as an abstract class and should never be directly
    instantiated.

    Classes inheriting from NDData should overload the ``propagate_*``
    methods, keeping the call signature the same. The propagate methods can
    assume that a `parent_nddata` attribute is present which links to the parent_nddata
    dataset, and take an `~astropy.nddata.NDData` instance as the positional
    argument, *not* an `~astropy.nddata.NDError` instance, because the
    `~astropy.nddata.NDData` instance can be used to access both the data and
    the errors (some propagations require the data values).
    '''

    __metaclass__ = abc.ABCMeta

    @property
    def parent_nddata(self):
        if self._parent_nddata is None:
            raise MissingDataAssociationException("Error is not associated with an NDData object")
        else:
            return self._parent_nddata

    @parent_nddata.setter
    def parent_nddata(self, value):
        self._parent_nddata = value

    @abc.abstractmethod
    def propagate_add(self, other_nddata, result_data):
        '''
        Propagate errors for addition.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrorsException
            Raised if the method does not know how to add the errors
        '''
        pass

    @abc.abstractmethod
    def propagate_subtract(self, other_nddata, result_data):
        '''
        Propagate errors for subtraction.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrorsException
            Raised if the method does not know how to add the errors
        '''
        pass

    @abc.abstractmethod
    def propagate_multiply(self, other_nddata, result_data):
        '''
        Propagate errors for mutliplication.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error
        '''
        pass

    @abc.abstractmethod
    def propagate_divide(self, other_nddata, result_data):
        '''
        Propagate errors for division.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error
        '''
        pass


class StandardDeviationError(NDError):
    '''
    A class for standard deviation errors
    '''

    def __init__(self, array=None, copy=True):
        if array is None:
            self.array = None
        else:
            self.array = np.array(array, copy=copy, subok=True)

    @property
    def parent_nddata(self):
        try:
            if self._parent_nddata is None:
                raise MissingDataAssociationException("Error is not associated with an NDData object")
            else:
                return self._parent_nddata
        except AttributeError:
            raise MissingDataAssociationException("Error is not associated with an NDData object")


    @parent_nddata.setter
    def parent_nddata(self, value):
        if self.array is None or value is None:
            self._parent_nddata = value
        else:
            if value.shape != self.array.shape:
                raise ValueError("parent shape does not match array data shape")

    @property
    def array(self):
        return self._array

    @array.setter
    def array(self, value):
        if value is not None:
            try:
                if value.shape != self.parent_nddata.shape:
                    raise ValueError("array shape does not match parent data shape")
            except MissingDataAssociationException:
                pass
        self._array = value

    def propagate_add(self, other_nddata, result_data):
        '''
        Propagate errors for addition.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrorsException
            Raised if the method does not know how to propagate the errors
        '''

        if not isinstance(other_nddata.error, StandardDeviationError):
            raise IncompatibleErrorsException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.error.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_error = StandardDeviationError()
        result_error.array = np.sqrt(self.array ** 2 + other_nddata.error.array ** 2)

        return result_error

    def __getitem__(self, item):
        '''
            Slice the standard deviation array using standard numpy slicing
        '''

        new_array = self.array[item]
        return self.__class__(new_array, copy=False)

    def propagate_subtract(self, other_nddata, result_data):
        '''
        Propagate errors for subtraction.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrorsException
            Raised if the method does not know how to propagate the errors
        '''

        if not isinstance(other_nddata.error, StandardDeviationError):
            raise IncompatibleErrorsException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.error.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_error = StandardDeviationError()
        result_error.array = np.sqrt(self.array ** 2 + other_nddata.error.array ** 2)

        return result_error

    def propagate_multiply(self, other_nddata, result_data):
        '''
        Propagate errors for mutliplication.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrorsException
            Raised if the method does not know how to propagate the errors
        '''

        if not isinstance(other_nddata.error, StandardDeviationError):
            raise IncompatibleErrorsException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.error.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_error = StandardDeviationError()
        result_error.array = np.sqrt((self.array / self.data) ** 2
                               + (other_nddata.error.array / other_nddata.data) ** 2) \
                               * result_data

        return result_error

    def propagate_divide(self, other_nddata, result_data):
        '''
        Propagate errors for division.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrorsException
            Raised if the method does not know how to propagate the errors
        '''

        if not isinstance(other_nddata.error, StandardDeviationError):
            raise IncompatibleErrorsException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.error.array is None:
            raise ValueError("standard deviation values are not set in other_nddata")

        result_error = StandardDeviationError()
        result_error.array = np.sqrt((self.array / self.data) ** 2
                               + (other_nddata.error.array / other_nddata.data) ** 2) \
                               * result_data

        return result_error
