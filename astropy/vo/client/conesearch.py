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

Get sorted catalog names containing 'SDSS' (not case-sensitive):

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

Several ways to check search status. Any error raised
does not terminate the process:

>>> async_search.is_alive()
>>> async_search.ready()
>>> async_search.successful()

Get search results. If no results are returned after
30 seconds, an error is raised but this does not
terminate the process. Otherwise, conesearch result
is returned and can be manipulated as above.

>>> async_result = async_search.get(timeout=30)
>>> cone_arr = async_result.array

If search is taking too long and going nowhere,
it can be forced to terminate:

>>> async_search.terminate()

"""
from __future__ import print_function, division

# STDLIB
import multiprocessing
import operator

# LOCAL
from . import vos_catalog

__all__ = ['AsyncConeSearch', 'conesearch', 'list_catalogs']

_SERVICE_TYPE = 'conesearch'

class ConeSearchError(Exception):
    pass

class AsyncConeSearch(object):
    def __init__(self, *args, **kwargs):
        """
        Perform a cone search asynchronously using
        :py:class:`multiprocessing.pool.AsyncResult`.

        Cone search will be forced to run in silent
        mode. Warnings are controled by :mod:`warnings`
        module.

        Parameters
        ----------
        args, kwargs : see `conesearch`

        """
        kwargs['verbose'] = False
        self.proc = multiprocessing.Pool()
        self.result = self.proc.apply_async(conesearch, args, kwargs)

    def __getattr__(self, what):
        """Expose :py:class:`multiprocessing.pool.AsyncResult` attributes."""
        return getattr(self.result, what)

    def is_alive(self):
        """Check if any of pool processes is running."""
        return reduce(operator.or_, [p.is_alive() for p in self.proc._pool])

    def terminate(self):
        """Terminate the process."""
        self.proc.terminate()

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

def _local_conversion(func, x):
    """Try `func(x)` and replace `ValueError` with `ConeSearchError`."""
    try:
        y = func(x)
    except ValueError as e:
        raise ConeSearchError(e.message)
    return y
