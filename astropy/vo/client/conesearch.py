# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Support basic VO conesearch capabilities.

Based on the `Simple Cone Search Version 1.03 Recommendation
<http://www.ivoa.net/Documents/REC/DAL/ConeSearch-20080222.html>`_.

Available databases are generated on the server side by
`astropy.vo.server.validate.check_conesearch_sites`.
Default database is 'conesearch.json', which can be changed
locally in a session via AstroPy configuration system.

*CONFIGURABLE PROPERTIES*

These properties are set via Astropy configuration system:

    * `astropy.vo.client.conesearch_dbname`
    * Also depends on properties set by `astropy.vo.client.vos_catalog`

Examples
--------
>>> from astropy.vo.client import conesearch

Perform cone search for 0.5 degree radius around 47 Tuc
(RA 6.088 deg, DEC -72.086 deg) with minimum verbosity,
if supported. The first catalog in the database to
successfully return a result that conforms to IVOA standards
is used. If running this for the first time, a copy of the
catalogs database will be downloaded to local cache.
To run this again without using cached data, set `cache`
keyword to `False` (local cache will not be updated):

>>> result = conesearch.conesearch(6.088, -72.086, 0.5, pedantic=True)

Extract Numpy array containing the matched objects. See
`numpy.ndarray` for available operations:

>>> cone_arr = result.array.data
>>> col_names = cone_arr.dtype.names
>>> n_rec = cone_arr.size
>>> ra_list = cone_arr[ col_names[1] ] # This depends on the catalog
>>> first_row = cone_arr[0]
>>> last_row = cone_arr[-1]
>>> first_ten_rows = cone_arr[:10]

Perform the same cone search as above but asynchronously:

>>> async_search = conesearch.AsyncConeSearch(
        6.088, -72.086, 0.5, pedantic=True)

Check search status:

>>> async_search.is_done()

Get search results with 30-sec time-out. If still not
done after 30 seconds, `None` is returned. Otherwise,
conesearch result is returned and can be manipulated as
above. If no `timeout` given, it waits till completion:

>>> async_result = async_search.get(timeout=30)
>>> cone_arr = async_result.array.data

If search is taking too long and going nowhere,
it can be forced to terminate (not recommended):

>>> async_search.terminate()

Estimate the execution time and the number of results for
the cone search above. The function naively assumes a
linear model, which might not be accurate for some cases.
It also uses the normal cone search function, not the
asynchronous version:

>>> t_est, n_est = conesearch.predict_search(
        result.url, 6.088, -72.086, 0.5, pedantic=True)

For debugging purpose, one can obtain the actual values
and compare with the prediction above. Keep in mind that
running this for every prediction would defeat the purpose
of the prediction itself:

>>> t_real, n_real = conesearch.conesearch_timer(
        6.088, -72.086, 0.5, catalog_db=result.url, pedantic=True)

If one is unable to obtain any results using the default
cone search database that only contains sites that cleanly
passed validation, one can use the AstroPy configuration
system to use another database containing sites with
validation warnings. One should use these sites with
caution:

>>> conesearch.CONESEARCH_DBNAME.set('conesearch_warn')

Find catalog names containing 'STSCI' and sort them
alphabetically:

>>> stsci_cats = conesearch.list_catalogs(match_string='stsci', sort=True)

Perform the same cone search as above using only STSCI
catalogs in this new database. Get results even if VO table
does not comform to IVOA standards, i.e., ignore validation
warnings. As the first catalog in the database to
successfully return a result is used, the order of catalog
names given to `catalog_db` is important:

>>> result = conesearch.conesearch(
        6.088, -72.086, 0.5, catalog_db=stsci_cats, pedantic=False)
>>> cone_arr = result.array.data

To see the catalog information for the access URL used
above (also see `astropy.vo.client.vos_catalog`):

>>> from astropy.vo.client import vos_catalog
>>> my_db = vos_catalog.get_remote_catalog_db(conesearch.CONESEARCH_DBNAME())
>>> my_cat = my_db.get_catalog_by_url(result.url)
>>> print(my_cat)

"""
from __future__ import print_function, division

# STDLIB
import time

# THIRD PARTY
import numpy as np

# LOCAL
from . import vos_catalog
from ...logger import log
from ...utils.misc import Future

# LOCAL CONFIG
from ...config.configuration import ConfigurationItem

__all__ = ['AsyncConeSearch', 'conesearch', 'list_catalogs', 'predict_search',
           'conesearch_timer']

CONESEARCH_DBNAME = ConfigurationItem('conesearch_dbname', 'conesearch',
                                      'Conesearch database name.')


class ConeSearchError(Exception):  # pragma: no cover
    pass


class AsyncConeSearch(Future):
    """
    Perform a cone search asynchronously using
    :py:class:`threading.Thread`.

    Cone search will be forced to run in silent
    mode. Warnings are controled by :mod:`warnings`
    module.

    .. seealso::

        `astropy.utils.misc.Future`

    Parameters
    ----------
    args, kwargs : see `conesearch`

    Returns
    -------
    value : see `conesearch`

    """
    def __init__(self, *args, **kwargs):
        kwargs['verbose'] = False
        Future.__init__(self, conesearch, *args, **kwargs)


def conesearch(ra, dec, sr, verb=1, **kwargs):
    """
    Do a cone search on the given catalog.

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

    kwargs : keywords for `astropy.vo.client.vos_catalog.call_vo_service`

    Returns
    -------
    value : see `astropy.vo.client.vos_catalog.call_vo_service`

    Raises
    ------
    ConeSearchError
        When invalid inputs are passed into cone search.

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
    """
    Return the available conesearch catalogs as a list of strings.
    These can be used for the *catalog_db* argument to
    :func:`conesearch`.

    Parameters
    ----------
    kwargs : keywords for `astropy.vo.client.vos_catalog.list_catalogs`

    """
    return vos_catalog.list_catalogs(CONESEARCH_DBNAME(), **kwargs)


def predict_search(url, *args, **kwargs):
    """
    Predict the execution time needed and the number of results
    for a cone search for the given access URL, position, and
    radius.

    A search is first done with zero radius to establish network
    latency time, which does not account for variability in
    latency. Since the first call can take longer than usual
    due to non-network activities, latency time is defined as
    half of this search time.

    Subsequent searches are done with starting and ending radii
    at 0.05 and 0.5 of the given radius, respectively. Searches
    that take longer than network latency time will be used to
    extrapolate the search time and number of results. Searches
    will stop when search time reaches `astropy.vo.client.vos_timeout`
    or radius reaches half of the given value.

    Extrapolation uses least-square straight line fitting,
    assuming linear increase of search time and number of results
    with radius, which might not be accurate for some cases. If
    there are less than three data points in the fit, it fails.

    Warnings (see :mod:`warnings`) are given upon:

        #. Fitted slope is negative.
        #. Estimated time exceeds `astropy.vo.client.vos_timeout`.
        #. Estimated number of results is 0 or exceeds 10,000.

    .. note::

        If `verbose` keyword is set to `True`, extra log INFO
        and plot will be displayed. The plot uses :mod:`matplotlib`.

        The predicted results are just *rough* estimates.

        Prediction is done using `conesearch`. Prediction for
        `AsyncConeSearch` is not supported.

    Parameters
    ----------
    url : string
        Cone Search access URL to use.

    args, kwargs : see `conesearch`

    Returns
    -------
    t_est : float
        Estimated time in seconds needed for the search.

    n_est : int
        Estimated number of results the search will yield.

    Raises
    ------
    ConeSearchError
        If input args are invalid or prediction fails.

    """
    if len(args) != 3 or args[2] <= 0:  # pragma: no cover
        raise ConeSearchError(
            'conesearch must have exactly 3 arguments and search radius '
            'has to be > 0.')

    ra, dec, sr = args
    verbose = kwargs.get('verbose', True)

    # First search with radius=0 to establish network latency time
    kwargs['catalog_db'] = url
    t_0, n_cur = conesearch_timer(ra, dec, 0, **kwargs)

    # Search properties for timer extrapolation
    min_datapoints = 3  # Minimum successful searches needed for extrapolation
    num_datapoints = 10  # Number of desired data points for extrapolation
    t_min = 0.5 * t_0  # Min time to be considered not due to network latency
    t_max = vos_catalog.TIMEOUT()  # Max time to be considered too long
    n_min = 1      # Min number of results to be considered valid (inclusive)
    n_max = 10000  # Max number of results to be considered valid (inclusive)
    sr_min = 0.05 * sr  # Min radius to start the timer
    sr_max = 0.5 * sr   # Max radius to stop the timer
    sr_step = (1.0 / num_datapoints) * (sr_max - sr_min)  # Radius step

    if verbose:  # pragma: no cover
        log.info('predict_search latency time = {} s'.format(t_min))
        log.info('predict_search radius step = {} deg'.format(sr_step))

    # Slowly increase radius to get data points for extrapolation
    t_cur = t_min
    sr_cur = sr_min
    sr_arr, t_arr, n_arr = [], [], []
    while t_cur < t_max and sr_cur < sr_max:
        t_cur, n_cur = conesearch_timer(ra, dec, sr_cur, **kwargs)

        if verbose:  # pragma: no cover
            log.info('predict_search took {} s with {} results '
                     'at {} deg radius'.format(t_cur, n_cur, sr_cur))

        if t_cur >= t_min:
            sr_arr.append(sr_cur)
            t_arr.append(t_cur)
            n_arr.append(n_cur)

        sr_cur += sr_step

    n_datapoints = len(sr_arr)
    if n_datapoints < min_datapoints:  # pragma: no cover
        raise ConeSearchError('predict_search only has {} data points; '
                              'unable to continue.'.format(n_datapoints))

    sr_arr = np.array(sr_arr)
    t_arr = np.array(t_arr)
    n_arr = np.array(n_arr)

    # Predict execution time
    t_est, t_fit = _extrapolate(sr_arr, t_arr, sr, ymax=t_max,
                                name='execution time', unit='s')

    # Predict number of results
    n_est, n_fit = _extrapolate(sr_arr, n_arr, sr, ymin=n_min, ymax=n_max,
                                name='number of results')

    if verbose:  # pragma: no cover
        from matplotlib import pyplot as plt
        fig = plt.figure()

        ax1 = fig.add_subplot(211)
        _plot_predictions(ax1, sr_arr, t_arr, t_fit, sr, t_est, 't (s)')

        ax2 = fig.add_subplot(212)
        _plot_predictions(ax2, sr_arr, n_arr, n_fit, sr, n_est, 'N')
        ax2.set_xlabel('radius (deg)')

        plt.draw()

    return t_est, int(n_est)


def conesearch_timer(*args, **kwargs):
    """
    Time a single conesearch. For use by `predict_search`.

    Parameters
    ----------
    args, kwargs : see `conesearch`

    Returns
    -------
    t : float
        Execution time in seconds.

    n : int
        Number of results.

    """
    t_beg = time.time()
    try:
        out_votable = conesearch(*args, **kwargs)
    except vos_catalog.VOSError:
        n = 0
    else:
        n = out_votable.array.size
    t_end = time.time()
    return t_end - t_beg, n


def _local_conversion(func, x):
    """Try `func(x)` and replace `ValueError` with `ConeSearchError`."""
    try:
        y = func(x)
    except ValueError as e:  # pragma: no cover
        raise ConeSearchError(e.message)
    return y


def _extrapolate(x_arr, y_arr, x, ymin=None, ymax=None, name='data', unit=''):
    """For use by `predict_search`."""
    a = np.polyfit(x_arr, y_arr, 1)
    p = np.poly1d(a)
    y = p(x)
    y_fit = p(x_arr)

    if a[0] < 0:
        log.warn('Fitted slope for {} is negative ({})'.format(name, a[0]))

    if ymin is not None and y < ymin:  # pragma: no cover
        log.warn('Predicted {} is less than {} {}'.format(name, ymin, unit))
    elif ymax is not None and y > ymax:  # pragma: no cover
        log.warn('Predicted {} is more than {} {}'.format(name, ymax, unit))

    return y, y_fit


def _plot_predictions(ax, x_arr, y_arr, y_fit, x, y,
                      ylabel):  # pragma: no cover
    """For use by `predict_search`."""
    ax.plot(x_arr, y_arr, 'kx-')
    ax.plot(x_arr, y_fit, 'b--')
    ax.scatter([x], [y], marker='o', c='r')
    ax.set_ylabel(ylabel)
