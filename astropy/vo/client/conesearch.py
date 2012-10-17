# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Support basic VO conesearch capabilities.

Based on the `Simple Cone Search Version 1.03 Recommendation
<http://www.ivoa.net/Documents/REC/DAL/ConeSearch-20080222.html>`_.

Examples
--------
>>> UNTIL HERE

"""
from __future__ import print_function, division

# LOCAL
from . import vos_catalog
from ...utils.misc import dict_soft_update

_SERVICE_TYPE = 'conesearch'
#_SERVICE_TYPE = 'conesearch_simple' # Initial version by mdroe

class ConeSearchError(Exception):
    pass

def conesearch(ra, dec, sr, catalog_db=None, pedantic=None, verb=1, **kwargs):
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

    catalog_db : str, optional
        See `astropy.vo.client.vos_catalog.call_vo_service`.

    pedantic : bool or `None`
        See `astropy.io.votable.table.parse`.

    verb : {1, 2, 3}
        Verbosity indicating how many columns are to be returned
        in the resulting table.  Support for this parameter by
        a Cone Search service implementation is optional. If the
        service supports the parameter, then when the value is 1,
        the response should include the bare minimum of columns
        that the provider considers useful in describing the
        returned objects. When the value is 3, the service should
        return all of the columns that are available for describing
        the objects. A value of 2 is intended for requesting a
        medium number of columns between the minimum and maximum
        (inclusive) that are considered by the provider to most
        typically useful to the user. When not provided, the server
        should respond as if `verb` is 2. If not supported, the
        service should ignore the parameter and should always
        return the same columns for every request.

    **kwargs : dictionary
        See `astropy.vo.client.vos_catalog.call_vo_service`.

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

    args = {}
    dict_soft_update(args, kwargs)
    dict_soft_update(args, {
            'RA': ra,
            'DEC': dec,
            'SR': sr,
            'VERB': verb})

    return vos_catalog.call_vo_service(
        _SERVICE_TYPE, catalog_db=catalog_db,
        pedantic=pedantic, kwargs=args)

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
