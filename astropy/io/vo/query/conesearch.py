"""
Support basic VO conesearch capabilities.

Based on the `Simple Cone Search Version 1.03 Recommendation
<http://www.ivoa.net/Documents/REC/DAL/ConeSearch-20080222.html>`_.

"""
from __future__ import division, absolute_import

# ASTROPY
from ...utils.misc import dict_soft_update

# LOCAL
from . import vos_catalog

def conesearch(ra, dec, sr, catalog_db=None, pedantic=None, verb=1, **kwargs):
    """
    Do a conesearch on the given catalog.

    Parameters
    ----------
    ra, dec : float
        Right-ascension and declination in the ICRS coordinate
        system for the position of the center of the cone to
        search, given in decimal degrees.

    sr : float
        Radius of the cone to search, given in decimal degrees.

    catalog_db : str, optional
        See `astropy.io.vo.query.vos_catalog.catalog_db`.

    pedantic : bool, optional
        See `astropy.io.vo.table.parse.pedantic`.

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

    **kwargs : keywords
        Additional keywords may be provided to pass along to the
        server. These arguments are specific to the particular
        catalog(s) being queried.

    """
    # Validate arguments
    ra = float(ra)
    dec = float(dec)
    sr = float(sr)
    verb = int(verb)
    assert verb in (1, 2, 3)

    args = {}
    dict_soft_update(args, kwargs)
    dict_soft_update(args, {
            'RA': ra,
            'DEC': dec,
            'SR': sr,
            'VERB': verb})

    return vos_catalog.call_vo_service(
        'conesearch', catalog_db=catalog_db,
        pedantic=pedantic, kwargs=args)


def list_catalogs():
    """
    Return the available conesearch catalogs as a list of strings.
    These can be used for the *catalog_db* argument to
    :func:`conesearch`.
    """
    return vos_catalog.list_catalogs('conesearch')
