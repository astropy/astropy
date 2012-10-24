# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Support basic VO conesearch capabilities.

Based on the `Simple Cone Search Version 1.03 Recommendation
<http://www.ivoa.net/Documents/REC/DAL/ConeSearch-20080222.html>`_.

.. note::

    See `astropy.vo.client.vos_catalog` for configurable
    properties.

Examples
--------
>>> from astropy.vo.client import conesearch

Get sorted catalog names containing 'SDSS' (not case-sensitive).
If running this for the first time, a copy of the catalogs
database will be downloaded to local cache. See `astropy.config`
for caching behavior:

>>> all_sdss_cat = conesearch.list_catalogs(match_string='sdss', sort=True)

Perform cone search for 0.1 arcsec radius around Alpha Centauri
(ICRS RA=14:39:36.204 DEC=-60:50:08.23), which are already converted
to decimal degrees in the example below. Only use catalogs containing
'SDSS' obtained from the example above, with minimum verbosity.
Get results even if VO table does not comform to IVOA standards.
The first catalog in the database to successfully return a result is
used, i.e., order of catalog names is important:

>>> result = conesearch.conesearch(219.900850, -60.835619, 2.78e-05,
                                   catalog_db=all_sdss_cat, pedantic=False)

Extract Numpy recarray containing the matched objects. See
`numpy.recarray` for available operations:

>>> cone_arr = result.array
>>> col_names = cone_arr.dtype.names
>>> n_rec = cone_arr.size
>>> ra_list = cone_arr[ col_names[1] ] # This depends on the catalog
>>> first_row = cone_arr[0]
>>> last_row = cone_arr[-1]
>>> first_ten_rows = cone_arr[:10]

Reload catalog names containing 'SDSS' from VO Service
Database online instead of the cached version (cache will not be
updated):

>>> all_sdss_cat = conesearch.list_catalogs(match_string='sdss', sort=True,
                                            cache=False)

Perform same cone search as above but asynchronously and
only use catalogs that conform to IVOA standards:

>>> async_search = conesearch.AsyncConeSearch(
        219.900850, -60.835619, 2.78e-05, catalog_db=all_sdss_cat,
        pedantic=True)

Check search status:

>>> async_search.is_done()

Wait for 30 seconds and get search results. If not done
after 30 seconds, `None` is returned. Otherwise, conesearch result
is returned and can be manipulated as above.

>>> async_result = async_search.get(timeout=30)
>>> cone_arr = async_result.array

If search is taking too long and going nowhere,
it can be forced to terminate:

>>> async_search.terminate()

Estimate the execution time and the number of results for the cone
search above. This uses the normal cone search function, not the
asynchronous version:

>>> t_est, n_est = conesearch.predict_search(
        219.900850, -60.835619, 2.78e-05, catalog_db=all_sdss_cat,
        pedantic=False)

Estimate the execution time and the number of results for a cone
search at a specific URL with different RA, DEC, and search radius:

>>> url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=NOMAD&'
>>> t_est, n_est = conesearch.predict_search(
        120, 20, 0.5, catalog_db=url, pedantic=False)

For debugging purpose, one can obtain the actual values and compare
with the prediction above. Keep in mind that running this for every
prediction would defeat the purpose of the prediction itself:

>>> t_real, n_real, url_used = conesearch.conesearch_timer(
        120, 20, 0.5, catalog_db=url, pedantic=False)

"""
from __future__ import print_function, division

# STDLIB
import time

# THIRD PARTY
import numpy

# LOCAL
from . import vos_catalog
from ...logger import log
from ...utils.misc import Future

__all__ = ['AsyncConeSearch', 'conesearch', 'list_catalogs', 'predict_search',
           'conesearch_timer']

_SERVICE_TYPE = 'conesearch'

class ConeSearchError(Exception):
    pass

class AsyncConeSearch(Future):
    def __init__(self, *args, **kwargs):
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
    ra  = _local_conversion(float, ra)
    dec = _local_conversion(float, dec)
    sr  = _local_conversion(float, sr)
    verb = _local_conversion(int, verb)
    if verb not in (1, 2, 3):
        raise ConeSearchError('Verbosity must be 1, 2, or 3')

    args = {'RA': ra, 'DEC': dec, 'SR': sr, 'VERB': verb}

    return vos_catalog.call_vo_service(_SERVICE_TYPE, kwargs=args, **kwargs)

def list_catalogs(**kwargs):
    """
    Return the available conesearch catalogs as a list of strings.
    These can be used for the *catalog_db* argument to
    :func:`conesearch`.

    Parameters
    ----------
    kwargs : keywords for `astropy.vo.client.vos_catalog.list_catalogs`

    """
    return vos_catalog.list_catalogs(_SERVICE_TYPE, **kwargs)

def predict_search(*args, **kwargs):
    """
    Predict the execution time needed and the number of results
    for a cone search for the given position and radius.

    A search is first done with zero radius to establish network
    latency time, which does not account for variability in latency.
    Since the first call can take longer than usual due to
    non-network activities, latency time is defined as half of
    this search time.

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
    if len(args) != 3 or args[2] <= 0:
        raise ConeSearchError(
            'conesearch must have exactly 3 arguments and search radius '
            'has to be > 0.')

    ra, dec, sr = args
    verbose = kwargs.get('verbose', True)

    # First search with radius=0 to establish network latency time
    # and obtain access URL.
    t_0, n_cur, url = conesearch_timer(ra, dec, 0, **kwargs)
    kwargs['catalog_db'] = url

    # Search properties for timer extrapolation
    min_datapoints = 3  # Minimum successful searches needed for extrapolation
    num_datapoints = 10 # Number of desired data points for extrapolation
    t_min = 0.5 * t_0  # Min time to be considered not due to network latency
    t_max = vos_catalog.TIMEOUT()  # Max time to be considered too long
    n_min = 1      # Min number of results to be considered valid (inclusive)
    n_max = 10000  # Max number of results to be considered valid (inclusive)
    sr_min = 0.05 * sr  # Min radius to start the timer
    sr_max = 0.5 * sr   # Max radius to stop the timer
    sr_step = (1.0/num_datapoints) * (sr_max - sr_min) # Radius step

    if verbose:
        log.info('predict_search latency time = {} s'.format(t_min))
        log.info('predict_search radius step = {} deg'.format(sr_step))

    # Slowly increase radius to get data points for extrapolation
    t_cur = t_min
    sr_cur = sr_min
    sr_arr, t_arr, n_arr = [], [], []
    while t_cur < t_max and sr_cur < sr_max:
        t_cur, n_cur, url = conesearch_timer(ra, dec, sr_cur, **kwargs)

        if verbose:
            log.info('predict_search took {} s with {} results '
                     'at {} deg radius'.format(t_cur, n_cur, sr_cur))

        if t_cur >= t_min:
            sr_arr.append(sr_cur)
            t_arr.append(t_cur)
            n_arr.append(n_cur)

        sr_cur += sr_step

    n_datapoints = len(sr_arr)
    if n_datapoints < min_datapoints:
        raise ConeSearchError('predict_search only has {} data points; '
                              'unable to continue.'.format(n_datapoints))

    sr_arr = numpy.array(sr_arr)
    t_arr = numpy.array(t_arr)
    n_arr = numpy.array(n_arr)

    # Predict execution time
    t_est, t_fit = _extrapolate(sr_arr, t_arr, sr, ymax=t_max,
                                name='execution time', unit='s')

    # Predict number of results
    n_est, n_fit = _extrapolate(sr_arr, n_arr, sr, ymin=n_min, ymax=n_max,
                                name='number of results')

    if verbose:
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
    For use by `predict_search`.

    Parameters
    ----------
    args, kwargs : see `conesearch`

    Returns
    -------
    t : float
        Execution time in seconds.

    n : int
        Number of results.

    url : string
        Access URL used for search. This is useful to
        identify which cone search service to use if
        `catalog_db` selection yields multiple services.

    """
    t_beg = time.time()
    out_votable = conesearch(*args, **kwargs)
    t_end = time.time()
    return t_end - t_beg, out_votable.array.size, out_votable.url

def _local_conversion(func, x):
    """Try `func(x)` and replace `ValueError` with `ConeSearchError`."""
    try:
        y = func(x)
    except ValueError as e:
        raise ConeSearchError(e.message)
    return y

def _extrapolate(x_arr, y_arr, x, ymin=None, ymax=None, name='data', unit=''):
    """For use by `predict_search`."""
    m, c = numpy.polyfit(x_arr, y_arr, 1)
    y = m * x + c
    y_fit = m * x_arr + c

    if m < 0:
        log.warn('Fitted slope for {} is negative ({})'.format(name, m))

    if ymin is not None and y < ymin:
        log.warn('Predicted {} is less than {} {}'.format(name, ymin, unit))
    elif ymax is not None and y > ymax:
        log.warn('Predicted {} is more than {} {}'.format(name, ymax, unit))

    return y, y_fit

def _plot_predictions(ax, x_arr, y_arr, y_fit, x, y, ylabel):
    """For use by `predict_search`."""
    ax.plot(x_arr, y_arr, 'kx-')
    ax.plot(x_arr, y_fit, 'b--')
    ax.scatter([x], [y], marker='o', c='r')
    ax.set_ylabel(ylabel)
