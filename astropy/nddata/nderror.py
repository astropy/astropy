import abc

__all__ = ['NDError']


class NDError(object):
    '''
    This is the base class for error classes used with NDData. It is
    implemented as an abstract class and should never be directly
    instantiated.
    '''

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def add(self, error):
        '''
        Add current instance errors with `error`.

        Parameters
        ----------
        error : NDError instance
            The error to add to the current instance

        Returns
        -------
        result : NDError instance
            The resulting error
        '''
        pass

    @abc.abstractmethod
    def subtract(self, error):
        '''
        Subtract `error` from current instance errors.

        Parameters
        ----------
        error : NDError instance
            The error to subtract from the current instance

        Returns
        -------
        result : NDError instance
            The resulting error
        '''
        pass
