# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Support VO Simple Cone Search capabilities."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import warnings

# THIRD-PARTY
import numpy as np

# LOCAL
from . import vos_catalog
from .async import AsyncBase
from .exceptions import ConeSearchError, VOSError
from ... import units as u
from ...coordinates import (ICRS, BaseCoordinateFrame, Longitude, Latitude,
                            SkyCoord)
from ...units import Quantity
from ...utils.timer import timefunc, RunTimePredictor
from ...utils.exceptions import AstropyUserWarning
from ...utils import data

__all__ = ['AsyncConeSearch', 'conesearch', 'AsyncSearchAll', 'search_all',
           'list_catalogs', 'predict_search', 'conesearch_timer']

# Skip these doctests for now;
# TODO: Add the ability to add py.test markers (such as remote_data) to
# doctests
__doctest_skip__ = ['AsyncConeSearch', 'AsyncSearchAll']


class AsyncConeSearch(AsyncBase):
    """Perform a Cone Search asynchronously and returns the result of the
    first successful query.

    .. note::

        See `~astropy.vo.client.async.AsyncBase` for more details.

    Parameters
    ----------
    args, kwargs : see :func:`conesearch`

    Examples
    --------
    >>> from astropy import coordinates as coord
    >>> from astropy import units as u
    >>> c = coord.ICRS(6.0223 * u.degree, -72.0814 * u.degree)
    >>> async_search = conesearch.AsyncConeSearch(
    ...     c, 0.5 * u.degree,
    ...     catalog_db='The PMM USNO-A1.0 Catalogue (Monet 1997) 1')

    Check search status:

    >>> async_search.running()
    True
    >>> async_search.done()
    False

    Get search results after a 30-second wait (not to be
    confused with `astropy.utils.data.Conf.remote_timeout` that
    governs individual Cone Search queries). If search is still not
    done after 30 seconds, `TimeoutError` is raised. Otherwise,
    Cone Search result is returned and can be manipulated as in
    :ref:`Simple Cone Search Examples <vo-sec-scs-examples>`.
    If no ``timeout`` keyword given, it waits until completion:

    >>> async_result = async_search.get(timeout=30)
    >>> cone_arr = async_result.array.data
    >>> cone_arr.size
    36184

    """
    def __init__(self, *args, **kwargs):
        super(AsyncConeSearch, self).__init__(conesearch, *args, **kwargs)


def conesearch(center, radius, verb=1, **kwargs):
    """Perform Cone Search and returns the result of the
    first successful query.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`, `~astropy.coordinates.BaseCoordinateFrame`, or sequence of length 2
        Position of the center of the cone to search. It may be specified as an
        object from the :ref:`astropy-coordinates` package, or as a length 2
        sequence. If a sequence, it is assumed to be ``(RA, DEC)`` in the
        ICRS coordinate frame, given in decimal degrees.

    radius : float or `~astropy.units.quantity.Quantity`
        Radius of the cone to search:

            - If float is given, it is assumed to be in decimal degrees.
            - If astropy quantity is given, it is internally converted
              to degrees.

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
              `astropy.vo.Conf.conesearch_dbname` catalogs is
              downloaded from `astropy.vo.Conf.vos_baseurl`.  The
              first catalog in the database to successfully return a
              result is used.

            - *catalog name*: A name in the database of
              `astropy.vo.Conf.conesearch_dbname` catalogs at
              `astropy.vo.Conf.vos_baseurl` is used.  For a list of
              acceptable names, use :func:`list_catalogs`.

            - *url*: The prefix of a URL to a IVOA Service for
              `astropy.vo.Conf.conesearch_dbname`.  Must
              end in either '?' or '&'.

            - `~astropy.vo.client.vos_catalog.VOSCatalog` object: A specific
              catalog manually downloaded and selected from the database
              (see :ref:`vo-sec-client-vos`).

            - Any of the above 3 options combined in a list, in which case
              they are tried in order.

    pedantic : bool or `None`
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning. Warnings may be controlled using
        :py:mod:`warnings` module. When not provided, uses the
        configuration setting `astropy.io.votable.Conf.pedantic`,
        which defaults to `False`.

    verbose : bool
        Verbose output.

    cache : bool
        Use caching for VO Service database. Access to actual VO
        websites referenced by the database still needs internet
        connection.

    Returns
    -------
    obj : `astropy.io.votable.tree.Table`
        First table from first successful VO service request.

    Raises
    ------
    ConeSearchError
        When invalid inputs are passed into Cone Search.

    VOSError
        If VO service request fails.

    """
    from .. import conf

    # Validate RA and DEC
    ra, dec = _validate_coord(center)

    # Validate search radius
    sr = _validate_sr(radius)

    # Validate verbosity
    verb = _local_conversion(int, verb)
    if verb not in (1, 2, 3):  # pragma: no cover
        raise ConeSearchError('Verbosity must be 1, 2, or 3')

    args = {'RA': ra, 'DEC': dec, 'SR': sr, 'VERB': verb}

    return vos_catalog.call_vo_service(conf.conesearch_dbname,
                                       kwargs=args, **kwargs)


class AsyncSearchAll(AsyncBase):
    """Perform a Cone Search asynchronously, storing all results
    instead of just the result from first successful query.

    .. note::

        See `~astropy.vo.client.async.AsyncBase` for more details.

    Parameters
    ----------
    args, kwargs : see :func:`search_all`

    Examples
    --------
    >>> from astropy import coordinates as coord
    >>> from astropy import units as u
    >>> c = coord.ICRS(6.0223 * u.degree, -72.0814 * u.degree)
    >>> async_search = conesearch.AsyncSearchAll(c, 0.5 * u.degree)

    Check search status:

    >>> async_search.running()
    True
    >>> async_search.done()
    False

    Get a dictionary of all search results after a 30-second wait (not
    to be confused with `astropy.utils.data.Conf.remote_timeout` that
    governs individual Cone Search queries). If search is still not
    done after 30 seconds, `TimeoutError` is raised. Otherwise, a
    dictionary is returned and can be manipulated as in :ref:`Simple
    Cone Search Examples <vo-sec-scs-examples>`.  If no ``timeout``
    keyword given, it waits until completion:

    >>> async_allresults = async_search.get(timeout=30)
    >>> all_catalogs = list(async_allresults)
    >>> first_cone_arr = async_allresults[all_catalogs[0]].array.data
    >>> first_cone_arr.size
    36184

    """
    def __init__(self, *args, **kwargs):
        AsyncBase.__init__(self, search_all, *args, **kwargs)


def search_all(*args, **kwargs):
    """Perform Cone Search and returns the results of
    all successful queries.

    .. warning::

        Could potentially take up significant run time and
        computing resources.

    Parameters
    ----------
    args, kwargs :
        Arguments and keywords accepted by :func:`conesearch`.

    Returns
    -------
    all_results : dict of `astropy.io.votable.tree.Table` objects
        A dictionary of tables from successful VO service requests,
        with keys being the access URLs. If none is successful,
        an empty dictionary is returned.

    Raises
    ------
    ConeSearchError
        When invalid inputs are passed into Cone Search.

    """
    from .. import conf

    all_results = {}

    catalog_db = kwargs.get('catalog_db', None)
    if 'catalog_db' in kwargs:
        kwargs.pop('catalog_db')

    cache = kwargs.get('cache', True)
    verbose = kwargs.get('verbose', True)

    catalogs = vos_catalog._get_catalogs(conf.conesearch_dbname, catalog_db,
                                         cache=cache, verbose=verbose)

    for name, catalog in catalogs:
        try:
            result = conesearch(catalog_db=catalog, *args, **kwargs)
        except VOSError:
            pass
        else:
            all_results[result.url] = result

    return all_results


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
    from .. import conf

    return vos_catalog.list_catalogs(conf.conesearch_dbname, **kwargs)


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
        #. Estimated run time exceeds
           `astropy.utils.data.Conf.remote_timeout`.

    .. note::

        If ``verbose=True``, extra log info will be provided.
        But unlike :func:`conesearch_timer`, timer info is suppressed.

        If ``plot=True``, plot will be displayed.
        Plotting uses `matplotlib <http://matplotlib.sourceforge.net/>`_.

        The predicted results are just *rough* estimates.

        Prediction is done using :func:`conesearch`. Prediction for
        `AsyncConeSearch` is not supported.

    Parameters
    ----------
    url : str
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
        If input parameters are invalid.

    VOSError
        If VO service request fails.

    """
    if len(args) != 2:  # pragma: no cover
        raise ConeSearchError('conesearch must have exactly 2 arguments')

    plot = kwargs.get('plot', False)
    if 'plot' in kwargs:  # pragma: no cover
        del kwargs['plot']

    center, radius = args
    sr = _validate_sr(radius)
    if sr <= 0:
        raise ConeSearchError('Search radius must be > 0 degrees')

    kwargs['catalog_db'] = url
    cs_pred = RunTimePredictor(conesearch, center, **kwargs)

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

    if t_est < 0 or t_coeffs[1] < 0:  # pragma: no cover
        warnings.warn(
            'Estimated runtime ({0} s) is non-physical with slope of '
            '{1}'.format(t_est, t_coeffs[1]), AstropyUserWarning)
    elif t_est > data.conf.remote_timeout:  # pragma: no cover
        warnings.warn(
            'Estimated runtime is longer than timeout of '
            '{0} s'.format(data.conf.remote_timeout), AstropyUserWarning)

    # Predict number of objects
    sr_arr = sorted(cs_pred.results)  # Orig with floating point error
    n_arr = [cs_pred.results[key].array.size for key in sr_arr]
    n_coeffs = np.polyfit(sr_arr, n_arr, 1)
    n_fitfunc = np.poly1d(n_coeffs)
    n_est = int(round(n_fitfunc(sr)))

    if n_est < 0 or n_coeffs[0] < 0:  # pragma: no cover
        warnings.warn('Estimated #objects ({0}) is non-physical with slope of '
                      '{1}'.format(n_est, n_coeffs[0]), AstropyUserWarning)

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
    with a single try and a verbose timer.

    Parameters
    ----------
    args, kwargs : see :func:`conesearch`

    Returns
    -------
    t : float
        Run time in seconds.

    obj : `astropy.io.votable.tree.Table`
        First table from first successful VO service request.

    """
    return conesearch(*args, **kwargs)


def _local_conversion(func, x):
    """Try ``func(x)`` and replace `ValueError` with
    ``ConeSearchError``."""
    try:
        y = func(x)
    except ValueError as e:  # pragma: no cover
        raise ConeSearchError(str(e))
    else:
        return y


def _validate_coord(center):
    """Validate coordinates."""
    if isinstance(center, SkyCoord):
        icrscoord = center.transform_to(ICRS).frame
    elif isinstance(center, BaseCoordinateFrame):
        icrscoord = center.transform_to(ICRS)
    else:
        icrscoord = ICRS(Longitude(center[0], unit=u.degree),
                         Latitude(center[1], unit=u.degree))

    return icrscoord.ra.degree, icrscoord.dec.degree


def _validate_sr(radius):
    """Validate search radius."""
    if isinstance(radius, Quantity):
        sr_angle = radius.to(u.degree)
    else:
        sr_angle = radius * u.degree

    return sr_angle.value
