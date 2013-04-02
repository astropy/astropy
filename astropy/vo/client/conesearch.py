# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Support basic VO Cone Search capabilities."""
from __future__ import print_function, division

# THIRD-PARTY
import numpy as np

# LOCAL
from . import vos_catalog
from ...config.configuration import ConfigurationItem
from ...logger import log
from ...utils.data import REMOTE_TIMEOUT
from ...utils.compat.futures import ThreadPoolExecutor
from ...utils.timer import timefunc, RunTimePredictor


__all__ = ['AsyncConeSearch', 'conesearch', 'list_catalogs', 'predict_search',
           'conesearch_timer']

CONESEARCH_DBNAME = ConfigurationItem('conesearch_dbname', 'conesearch_good',
                                      'Conesearch database name.')


class ConeSearchError(Exception):  # pragma: no cover
    pass


class AsyncConeSearch(object):
    """Perform a Cone Search asynchronously using
    :py:class:`concurrent.futures.ThreadPoolExecutor`.

    Cone Search will be forced to run in silent
    mode. Warnings are controled by :py:mod:`warnings`
    module.

    .. note::

        Methods of the attributes can be accessed directly,
        with priority given to ``executor``.

    Parameters
    ----------
    args, kwargs : see :func:`conesearch`

    Attributes
    ----------
    executor : :py:class:`concurrent.futures.ThreadPoolExecutor`
        Executor running :func:`conesearch` on single thread.

    future : :py:class:`concurrent.futures.Future`
        Asynchronous execution created by ``executor``.

    Examples
    --------
    >>> async_search = conesearch.AsyncConeSearch(6.088, -72.086, 0.5)

    Check search status:

    >>> async_search.running()
    True
    >>> async_search.done()
    False


    Get search results after a 30-second wait (not to be
    confused with ``astropy.utils.data.REMOTE_TIMEOUT`` that
    governs individual Cone Search queries). If search is still not
    done after 30 seconds, ``TimeoutError`` is raised. Otherwise,
    Cone Search result is returned and can be manipulated as in
    :ref:`Simple Cone Search Examples <vo-sec-scs-examples>`.
    If no ``timeout`` keyword given, it waits until completion:

    >>> async_result = async_search.get(timeout=30)
    >>> cone_arr = async_result.array.data
    >>> cone_arr.size
    36386

    """
    def __init__(self, *args, **kwargs):
        kwargs['verbose'] = False
        self.executor = ThreadPoolExecutor(1)
        self.future = self.executor.submit(conesearch, *args, **kwargs)

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
        result : `astropy.io.votable.tree.Table` object
            First table from first successful VO service request.

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


def conesearch(ra, dec, sr, verb=1, **kwargs):
    """Do a Cone Search on the given catalog.

    Parameters
    ----------
    ra, dec : float
        Right-ascension and declination in the ICRS coordinate
        system for the position of the center of the cone to
        search, given in decimal degrees.

    sr : float
        Radius of the cone to search, given in decimal degrees.

    verb : {1, 2, 3}
        Verbosity indicating how many columns are to be returned
        in the resulting table. Support for this parameter by
        a Cone Search service implementation is optional. If the
        service supports the parameter:

            1. Return the bare minimum number of columns that
               the provider considers useful in describing the
               returned objects.
            2. Return a medium number of columns between the
               minimum and maximum (inclusive) that are
               considered by the provider to most typically
               useful to the user.
            3. Return all of the columns that are available for
               describing the objects.

        If not supported, the service should ignore the parameter
        and always return the same columns for every request.

    catalog_db
        May be one of the following, in order from easiest to
        use to most control:

            - `None`: A database of
              ``astropy.vo.client.conesearch.CONESEARCH_DBNAME``
              catalogs is downloaded from
              ``astropy.vo.client.vos_catalog.BASEURL``.  The first
              catalog in the database to successfully return a result is used.

            - *catalog name*: A name in the database of
              ``astropy.vo.client.conesearch.CONESEARCH_DBNAME``
              catalogs at ``astropy.vo.client.vos_catalog.BASEURL`` is used.
              For a list of acceptable names, see :func:`list_catalogs`.

            - *url*: The prefix of a URL to a IVOA Service for
              ``astropy.vo.client.conesearch.CONESEARCH_DBNAME``.
              Must end in either '?' or '&'.

            - `VOSCatalog` object: A specific catalog manually downloaded and
              selected from the database (see :ref:`vo-sec-client-vos`).

            - Any of the above 3 options combined in a list, in which case
              they are tried in order.

    pedantic : bool or `None`
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        :py:mod:`warnings` module.
        When not provided, uses the configuration setting
        ``astropy.io.votable.table.PEDANTIC``, which defaults to `False`.

    verbose : bool
        Verbose output.

    cache : bool
        Use caching for VO Service database. Access to actual VO
        websites referenced by the database still needs internet
        connection.

    Returns
    -------
    obj : `astropy.io.votable.tree.Table` object
        First table from first successful VO service request.

    Raises
    ------
    ConeSearchError
        When invalid inputs are passed into Cone Search.

    VOSError
        If VO service request fails.

    """
    # Validate arguments
    ra = _local_conversion(float, ra)
    dec = _local_conversion(float, dec)
    sr = _local_conversion(float, sr)
    verb = _local_conversion(int, verb)
    if verb not in (1, 2, 3):  # pragma: no cover
        raise ConeSearchError('Verbosity must be 1, 2, or 3')

    args = {'RA': ra, 'DEC': dec, 'SR': sr, 'VERB': verb}

    return vos_catalog.call_vo_service(CONESEARCH_DBNAME(),
                                       kwargs=args, **kwargs)


def list_catalogs(**kwargs):
    """Return the available Cone Search catalogs as a list of strings.
    These can be used for the ``catalog_db`` argument to
    :func:`conesearch`.

    Parameters
    ----------
    cache : bool
        Use caching for VO Service database. Access to actual VO
        websites referenced by the database still needs internet
        connection.

    verbose : bool
        Show download progress bars.

    pattern : str or `None`
        If given string is anywhere in a catalog name, it is
        considered a matching catalog. It accepts patterns as
        in :py:mod:`fnmatch` and is case-insensitive.
        By default, all catalogs are returned.

    sort : bool
        Sort output in alphabetical order. If not sorted, the
        order depends on dictionary hashing. Default is `True`.

    Returns
    -------
    arr : list of str
        List of catalog names.

    """
    return vos_catalog.list_catalogs(CONESEARCH_DBNAME(), **kwargs)


def predict_search(url, *args, **kwargs):
    """Predict the run time needed and the number of objects
    for a Cone Search for the given access URL, position, and
    radius.

    Run time prediction uses `astropy.utils.timer.RunTimePredictor`.
    Baseline searches are done with starting and ending radii at
    0.05 and 0.5 of the given radius, respectively.

    Extrapolation on good data uses least-square straight line fitting,
    assuming linear increase of search time and number of objects
    with radius, which might not be accurate for some cases. If
    there are less than 3 data points in the fit, it fails.

    Warnings (controlled by :py:mod:`warnings`) are given when:

        #. Fitted slope is negative.
        #. Any of the estimated results is negative.
        #. Estimated run time exceeds ``astropy.utils.data.REMOTE_TIMEOUT``.

    .. note::

        If ``verbose=True``, extra log info will be provided.
        But unlike :func:`conesearch_timer`, timer info is suppressed.

        If ``plot=True``, plot will be displayed.
        Plotting uses :mod:`matplotlib`.

        The predicted results are just *rough* estimates.

        Prediction is done using :func:`conesearch`. Prediction for
        `AsyncConeSearch` is not supported.

    Parameters
    ----------
    url : string
        Cone Search access URL to use.

    args, kwargs : see :func:`conesearch`
        Extra keyword ``plot`` is allowed and only used by this
        function and not :func:`conesearch`.

    Returns
    -------
    t_est : float
        Estimated time in seconds needed for the search.

    n_est : int
        Estimated number of objects the search will yield.

    Raises
    ------
    AssertionError
        If prediction fails.

    ConeSearchError
        If input paramters are invalid.

    VOSError
        If VO service request fails.

    """
    if len(args) != 3 or args[2] <= 0:  # pragma: no cover
        raise ConeSearchError(
            'conesearch must have exactly 3 arguments and search radius '
            'has to be > 0.')

    plot = kwargs.get('plot', False)
    if 'plot' in kwargs:  # pragma: no cover
        del kwargs['plot']

    ra, dec, sr = args

    kwargs['catalog_db'] = url
    cs_pred = RunTimePredictor(conesearch, ra, dec, **kwargs)

    # Search properties for timer extrapolation
    num_datapoints = 10  # Number of desired data points for extrapolation
    sr_min = 0.05 * sr  # Min radius to start the timer
    sr_max = 0.5 * sr   # Max radius to stop the timer
    sr_step = (1.0 / num_datapoints) * (sr_max - sr_min)  # Radius step

    # Slowly increase radius to get data points for extrapolation
    sr_arr = np.arange(sr_min, sr_max + sr_step, sr_step)
    cs_pred.time_func(sr_arr)

    # Predict run time
    t_coeffs = cs_pred.do_fit()
    t_est = cs_pred.predict_time(sr)

    if t_est < 0 or t_coeffs[0] < 0:  # pragma: no cover
        log.warn('Estimated runtime ({0} s) is non-physical with slope of '
                 '{1}'.format(t_est, t_coeffs[0]))
    elif t_est > REMOTE_TIMEOUT():  # pragma: no cover
        log.warn('Estimated runtime is longer than timeout of '
                 '{0} s'.format(REMOTE_TIMEOUT()))

    # Predict number of objects
    sr_arr = sorted(cs_pred.results)  # Orig with floating point error
    n_arr = [cs_pred.results[key].array.size for key in sr_arr]
    n_coeffs = np.polyfit(sr_arr, n_arr, 1)
    n_fitfunc = np.poly1d(n_coeffs)
    n_est = int(round(n_fitfunc(sr)))

    if n_est < 0 or n_coeffs[0] < 0:  # pragma: no cover
        log.warn('Estimated #objects ({0}) is non-physical with slope of '
                 '{1}'.format(n_est, n_coeffs[0]))

    if plot:  # pragma: no cover
        import matplotlib.pyplot as plt

        xlabeltext = 'radius (deg)'
        sr_fit = np.append(sr_arr, sr)
        n_fit = n_fitfunc(sr_fit)

        cs_pred.plot(xlabeltext=xlabeltext)

        fig, ax = plt.subplots()
        ax.plot(sr_arr, n_arr, 'kx-', label='Actual')
        ax.plot(sr_fit, n_fit, 'b--', label='Fit')
        ax.scatter([sr], [n_est], marker='o', c='r', label='Predicted')
        ax.set_xlabel(xlabeltext)
        ax.set_ylabel('#objects')
        ax.legend(loc='best', numpoints=1)
        plt.draw()

    return t_est, n_est


@timefunc(1)
def conesearch_timer(*args, **kwargs):
    """Time a single Cone Search using `astropy.utils.timer.timefunc`
    with single try and verbose timer.

    Parameters
    ----------
    args, kwargs : see :func:`conesearch`

    Returns
    -------
    t : float
        Run time in seconds.

    obj : `astropy.io.votable.tree.Table` object
        First table from first successful VO service request.

    """
    return conesearch(*args, **kwargs)


def _local_conversion(func, x):
    """Try ``func(x)`` and replace ``ValueError`` with ``ConeSearchError``."""
    try:
        y = func(x)
    except ValueError as e:  # pragma: no cover
        raise ConeSearchError(str(e))
    else:
        return y
