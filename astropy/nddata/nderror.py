import abc

__all__ = ['IncompatibleErrors', 'NDError']


class IncompatibleErrors(Exception):
    pass


class NDError(object):
    '''
    This is the base class for error classes used with NDData. It is
    implemented as an abstract class and should never be directly
    instantiated.

    Classes inheriting from NDData should implement an __init__ method that
    takes the parent NDData object as the first argument, and should overload
    the ``propagate_*`` methods, keeping the call signature the same.
    '''

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, parent):
        self._parent = parent

    @abc.abstractmethod
    def propagate_add(self, data, error):
        '''
        Propagate errors for addition.

        Parameters
        ----------
        data : NDData instance
            The data for the second operand in a + b
        error : NDError instance
            The error on the second operand in a + b

        Returns
        -------
        result : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrors
            Raised if the method does not know how to add the errors
        '''
        pass

    @abc.abstractmethod
    def propagate_subtract(self, data, error):
        '''
        Propagate errors for subtraction.

        Parameters
        ----------
        data : NDData instance
            The data for the second operand in a + b
        error : NDError instance
            The error on the second operand in a - b

        Returns
        -------
        result : NDError instance
            The resulting error

        Raises
        ------
        IncompatibleErrors
            Raised if the method does not know how to add the errors
        '''
        pass

    @abc.abstractmethod
    def propagate_multiply(self, data, error):
        '''
        Propagate errors for mutliplication.

        Parameters
        ----------
        data : NDData instance
            The data for the second operand in a + b
        error : NDError instance
            The error on the second operand in a * b

        Returns
        -------
        result : NDError instance
            The resulting error
        '''
        pass

    @abc.abstractmethod
    def propagate_divide(self, data, error):
        '''
        Propagate errors for division.

        Parameters
        ----------
        data : NDData instance
            The data for the second operand in a + b
        error : NDError instance
            The error on the second operand in a / b

        Returns
        -------
        result : NDError instance
            The resulting error
        '''
        pass
