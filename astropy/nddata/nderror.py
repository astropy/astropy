import abc

import numpy as np

__all__ = ['MissingDataAssociation', 'IncompatibleErrors', 'NDError', 'StandardDeviationError']


class IncompatibleErrors(Exception):
    pass


class MissingDataAssociation(Exception):
    pass


class NDError(object):
    '''
    This is the base class for error classes used with NDData. It is
    implemented as an abstract class and should never be directly
    instantiated.

    Classes inheriting from NDData should overload the ``propagate_*``
    methods, keeping the call signature the same. The propagate methods can
    assume that a `parent` attribute is present which links to the parent
    dataset, and take an `~astropy.nddata.NDData` instance as the positional
    argument, *not* an `~astropy.nddata.NDError` instance, because the
    `~astropy.nddata.NDData` instance can be used to access both the data and
    the errors (some propagations require the data values).
    '''

    __metaclass__ = abc.ABCMeta

    @property
    def parent(self):
        if self._parent is None:
            raise MissingDataAssociation("Error is not associated with an NDData object")
        else:
            return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = value

    @abc.abstractmethod
    def propagate_add(self, operand, result):
        '''
        Propagate errors for addition.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrors
            Raised if the method does not know how to add the errors
        '''
        pass

    @abc.abstractmethod
    def propagate_subtract(self, operand, result):
        '''
        Propagate errors for subtraction.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrors
            Raised if the method does not know how to add the errors
        '''
        pass

    @abc.abstractmethod
    def propagate_multiply(self, operand, result):
        '''
        Propagate errors for mutliplication.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error
        '''
        pass

    @abc.abstractmethod
    def propagate_divide(self, operand, result):
        '''
        Propagate errors for division.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

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

    def __init__(self, array=None, parent=None):

        # First initialize internal values to None
        self._parent = None
        self._array = None

        self.array = array
        self.parent = parent

    @property
    def parent(self):
        if self._parent is None:
            raise MissingDataAssociation("Error is not associated with an NDData object")
        else:
            return self._parent

    @parent.setter
    def parent(self, value):
        if self.array is None or value is None:
            self._parent = value
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
                if value.shape != self.parent.shape:
                    raise ValueError("array shape does not match parent data shape")
            except MissingDataAssociation:
                pass
        self._array = value

    def propagate_add(self, operand, result):
        '''
        Propagate errors for addition.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrors
            Raised if the method does not know how to add the errors
        '''

        if not isinstance(operand.error, StandardDeviationError):
            raise IncompatibleErrors

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if operand.error.array is None:
            raise ValueError("standard deviation values are not set in operand")

        result_error = StandardDeviationError(parent=result)
        result_error.array = np.sqrt(self.array ** 2 + operand.error.array ** 2)

        return result_error

    def propagate_subtract(self, operand, result):
        '''
        Propagate errors for subtraction.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrors
            Raised if the method does not know how to add the errors
        '''

        if not isinstance(operand.error, StandardDeviationError):
            raise IncompatibleErrors

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if operand.error.array is None:
            raise ValueError("standard deviation values are not set in operand")

        result_error = StandardDeviationError(parent=result)
        result_error.array = np.sqrt(self.array ** 2 + operand.error.array ** 2)

        return result_error

    def propagate_multiply(self, operand, result):
        '''
        Propagate errors for mutliplication.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error
        '''

        if not isinstance(operand.error, StandardDeviationError):
            raise IncompatibleErrors

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if operand.error.array is None:
            raise ValueError("standard deviation values are not set in operand")

        result_error = StandardDeviationError(parent=result)
        result_error.array = np.sqrt((self.array / self.data) ** 2
                               + (operand.error.array / operand.data) ** 2) \
                               * result.data

        return result_error

    def propagate_divide(self, operand, result):
        '''
        Propagate errors for division.

        Parameters
        ----------
        operand : NDData instance
            The data for the second operand in a + b
        result : NDData instance
            The data object that is the result of the addition

        Returns
        -------
        result_error : NDError instance
            The resulting error
        '''
        if not isinstance(operand.error, StandardDeviationError):
            raise IncompatibleErrors

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if operand.error.array is None:
            raise ValueError("standard deviation values are not set in operand")

        result_error = StandardDeviationError(parent=result)
        result_error.array = np.sqrt((self.array / self.data) ** 2
                               + (operand.error.array / operand.data) ** 2) \
                               * result.data

        return result_error
