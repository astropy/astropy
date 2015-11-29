# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Asynchronous VO service requests."""
from __future__ import absolute_import, division, print_function, unicode_literals

# LOCAL
from ...utils.compat.futures import ThreadPoolExecutor


__all__ = ['AsyncBase']


class AsyncBase(object):
    """Base class for asynchronous VO service requests
    using :py:class:`concurrent.futures.ThreadPoolExecutor`.

    Service request will be forced to run in silent
    mode by setting ``verbose=False``. Warnings are controlled
    by :py:mod:`warnings` module.

    .. note::

        Methods of the attributes can be accessed directly,
        with priority given to ``executor``.

    Parameters
    ----------
    func : function
        The function to run.

    args, kwargs
        Arguments and keywords accepted by the service request
        function to be called asynchronously.

    Attributes
    ----------
    executor : :py:class:`concurrent.futures.ThreadPoolExecutor`
        Executor running the function on single thread.

    future : :py:class:`concurrent.futures.Future`
        Asynchronous execution created by ``executor``.

    """
    def __init__(self, func, *args, **kwargs):
        kwargs['verbose'] = False
        self.executor = ThreadPoolExecutor(1)
        self.future = self.executor.submit(func, *args, **kwargs)

    def __getattr__(self, what):
        """Expose ``executor`` and ``future`` methods."""
        try:
            return getattr(self.executor, what)
        except AttributeError:
            return getattr(self.future, what)

    def get(self, timeout=None):
        """Get result, if available, then shut down thread.

        Parameters
        ----------
        timeout : int or float
            Wait the given amount of time in seconds before
            obtaining result. If not given, wait indefinitely
            until function is done.

        Returns
        -------
        result
            Result returned by the function.

        Raises
        ------
        Exception
            Errors raised by :py:class:`concurrent.futures.Future`.

        """
        try:
            result = self.future.result(timeout=timeout)
        except Exception as e:  # pragma: no cover
            result = None
            raise e
        finally:
            self.executor.shutdown(wait=False)
            return result
