# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Common utilities for accessing VO simple services.

When run for the first time, catalogs database will be
downloaded to a local cache. Otherwise, cached data
are used. To by-pass cached data, some functions have
a `cache` keyword that can be set to `False`.

*CONFIGURABLE PROPERTIES*

These properties are set via Astropy configuration system:

    * `astropy.io.votable.pedantic`
    * `astropy.vo.client.vos_baseurl`
    * `astropy.vo.client.vos_timeout`

Examples
--------
Get all catalogs from a database named 'conesearch_good'.
For more cone search examples, see `astropy.vo.client.conesearch`:

>>> from astropy.vo.client import vos_catalog
>>> my_db = vos_catalog.get_remote_catalog_db('conesearch_good')

Find catalog names containing 'usno*a2':

>>> cat_names = my_db.list_catalogs(pattern='usno*a2', sort=True)

Get information for first catalog from above. Catalog
fields may vary:

>>> my_cat = my_db.get_catalog(cat_names[0])
>>> print(my_cat)
{
    \"validate_network_error\": null,
    \"capabilityClass\": \"ConeSearch\",
    \"updated\": \"2011-09-14T20:20:21\",
    \"capabilityValidationLevel\": \"\",
    \"maxRecords\": 9999,
    # ...
    \"validate_xmllint\": true
}
>>> cat_keys = my_cat.keys()
>>> cat_url = my_cat['url']
>>> max_rad = my_cat['maxRadius']

One can also get information for a catalog using its URL.
If URL yields multiple catalogs, only first match is given:

>>> my_cat2 = my_db.get_catalog_by_url(cat_url)

To get all the matching catalogs by URL:

>>> cat_list = [cat for key, cat in my_db.get_catalogs_by_url(cat_url)]

To get all catalogs in the database:

>>> all_cat_list = [cat for key, cat in my_db.get_catalogs()]

Other functions listed in this module but not in the examples
above are used indirectly by other modules in `astropy.vo`
package.

"""
from __future__ import print_function, division

# STDLIB
import json

# LOCAL
from ...config.configuration import ConfigurationItem
from ...io.votable import table
from ...io.votable.exceptions import vo_raise, vo_warn, E19, W24, W25
from ...utils import webquery
from ...utils.console import color_print
from ...utils.data import get_readable_fileobj


__all__ = ['VOSCatalog', 'VOSDatabase', 'get_remote_catalog_db',
           'call_vo_service', 'list_catalogs']

__dbversion__ = 1

VO_PEDANTIC = table.PEDANTIC()

BASEURL = ConfigurationItem('vos_baseurl',
                            'http://stsdas.stsci.edu/astrolib/vo_databases/',
                            'URL where VO Service database file is stored.')

TIMEOUT = ConfigurationItem('vos_timeout', 30.0,
                            'Timeout in seconds for VO Service query')


class VOSError(Exception):  # pragma: no cover
    pass


class VOSCatalog(object):
    """
    A class to represent VO Service Catalog.

    Parameters
    ----------
    tree : JSON tree

    """
    def __init__(self, tree):
        self._tree = tree

    def __repr__(self):
        """Pretty print."""
        return json.dumps(self._tree, indent=4)  # pragma: no cover

    def __getattr__(self, what):
        """Expose dictionary attributes."""
        return getattr(self._tree, what)

    def __getitem__(self, what):
        """Expose dictionary key look-up."""
        return self._tree[what]


class VOSDatabase(VOSCatalog):
    """
    A class to represent a collection of `VOSCatalog`.

    Parameters
    ----------
    tree : JSON tree

    Raises
    ------
    VOSError
        If given `tree` does not have 'catalogs' key.

    """
    def __init__(self, tree):
        self._tree = tree

        if tree['__version__'] > __dbversion__:  # pragma: no cover
            vo_warn(W24)

        if not 'catalogs' in tree:  # pragma: no cover
            raise VOSError("Invalid VO service catalog database")

        self._catalogs = tree['catalogs']

    def get_catalogs(self):
        """Iterator to get all catalogs."""
        for key, val in self._catalogs.items():
            yield key, VOSCatalog(val)

    def get_catalogs_by_url(self, url):
        """Like `get_catalogs` but using access URL look-up."""
        for key, cat in self.get_catalogs():
            if cat['url'] == url:
                yield key, cat

    def get_catalog(self, name):
        """
        Get one catalog of given name.

        Parameters
        ----------
        name : str
            Primary key identifying the catalog.

        Returns
        -------
        obj : `VOSCatalog` object

        Raises
        ------
        VOSError
            If catalog is not found.

        """
        if not name in self._catalogs:
            raise VOSError("No catalog '{0}' found.".format(name))
        return VOSCatalog(self._catalogs[name])

    def get_catalog_by_url(self, url):
        """
        Like `get_catalog` but using access URL look-up.
        On multiple matches, only first match is returned.

        """
        out_cat = None
        for key, cat in self.get_catalogs_by_url(url):
            out_cat = cat
            break
        if out_cat is None:
            raise VOSError("No catalog with URL '{0}' found.".format(url))
        return out_cat

    def list_catalogs(self, pattern=None, sort=False):
        """
        List of catalog names.

        Parameters
        ----------
        pattern : str or `None`
            If given string is anywhere in a catalog name, it is
            considered a matching catalog. It accepts patterns as
            in :mod:`fnmatch` and is case-insensitive.
            By default, all catalogs are returned.

        sort : bool
            Sort output in alphabetical order. If not sorted, the
            order depends on dictionary hashing.

        """
        all_catalogs = self._catalogs.keys()

        if pattern is None or len(all_catalogs) == 0:
            out_arr = all_catalogs
        else:
            import fnmatch
            import re
            pattern = re.compile(fnmatch.translate('*' + pattern + '*'),
                                 re.IGNORECASE)
            out_arr = [s for s in all_catalogs if pattern.match(s)]

        if sort:
            out_arr.sort()

        return out_arr


def get_remote_catalog_db(dbname, cache=True):
    """
    Get a database of VO services (which is a JSON file) from a remote
    location.

    Parameters
    ----------
    dbname : str
        Prefix of JSON file to download from `astropy.vo.client.vos_baseurl`.

    cache : bool
        Use caching for VO Service database. Access to actual VO
        websites referenced by the database still needs internet
        connection.

    Returns
    -------
    obj : `VOSDatabase` object

    """
    with get_readable_fileobj(BASEURL() + dbname + '.json',
                              encoding='utf8', cache=cache) as fd:
        tree = json.load(fd)

    return VOSDatabase(tree)


def _vo_service_request(url, pedantic, kwargs):
    req = webquery.webget_open(url, timeout=TIMEOUT(), **kwargs)
    try:
        tab = table.parse(req, filename=req.geturl(), pedantic=pedantic)
    finally:
        req.close()

    return vo_tab_parse(tab, url, kwargs)


def vo_tab_parse(tab, url, kwargs):
    """
    In case of errors from the server, a complete and correct
    'stub' VOTable file may still be returned.  This is to
    detect that case.

    Parameters
    ----------
    tab : `astropy.io.votable.tree.VOTableFile` object

    url : string
        URL used to obtain `tab`.

    kwargs : dict
        Keywords used to obtain `tab`, if any.

    Returns
    -------
    out_tab : `astropy.io.votable.tree.Table` object

    Raises
    ------
    IndexError
        Table iterator fails.

    VOSError
        Server returns error message or invalid table.

    """
    for param in tab.iter_fields_and_params():
        if param.ID.lower() == 'error':
            raise VOSError("Catalog server '{0}' returned error '{1}'".format(
                url, param.value))
        break

    if tab.resources == []:
        vo_raise(E19)

    for info in tab.resources[0].infos:
        if info.name == 'QUERY_STATUS' and info.value != 'OK':
            if info.content is not None:
                long_descr = ':\n{0}'.format(info.content)
            else:
                long_descr = ''
            raise VOSError("Catalog server '{0}' returned status "
                           "'{1}'{2}".format(url, info.value, long_descr))
        break

    out_tab = tab.get_first_table()

    kw_sr = [k for k in kwargs if 'sr' == k.lower()]
    if len(kw_sr) == 0:
        sr = 0
    else:
        sr = kwargs.get(kw_sr[0])

    if sr != 0 and out_tab.array.size <= 0:
        raise VOSError("Catalog server '{0}' returned {1} result".format(
            url, out_tab.array.size))

    out_tab.url = url  # Track the URL
    return out_tab


def call_vo_service(service_type, catalog_db=None, pedantic=None,
                    verbose=True, cache=True, kwargs={}):
    """
    Makes a generic VO service call.

    Parameters
    ----------
    service_type : str
        Name of the type of service, e.g., 'conesearch'.
        Used in error messages and to select a catalog database
        if `catalog_db` is not provided.

    catalog_db
        May be one of the following, in order from easiest to
        use to most control:

        - `None`
              A database of `service_type` catalogs is downloaded from
              `astropy.vo.client.vos_baseurl`.  The first catalog
              in the database to successfully return a result is used.

        - *catalog name*
              A name in the database of `service_type` catalogs
              at `astropy.vo.client.vos_baseurl` is used.  For a list
              of acceptable names, see :func:`list_catalogs`.

        - *url*
              The prefix of a *url* to a IVOA Service for `service_type`.
              Must end in either '?' or '&'.

        - `VOSCatalog` object
              A specific catalog manually downloaded and selected
              from the database using the APIs in
              `~astropy.vo.client.vos_catalog`.

        - Any of the above 3 options combined in a list, in which case
          they are tried in order.

    pedantic : bool or `None`
        See  `astropy.io.votable.table.parse`.

    verbose : bool
        Verbose output.

    cache : bool
        See `get_remote_catalog_db`.

    kwargs : dictionary
        Keyword arguments to pass to the catalog service.
        No checking is done that the arguments are accepted by
        the service, etc.

    Returns
    -------
    obj : `astropy.io.votable.tree.Table` object
        First table from first successful VO service request.

    Raises
    ------
    VOSError
        If VO service request fails.

    """
    if catalog_db is None:
        catalog_db = get_remote_catalog_db(service_type, cache=cache)
        catalogs = catalog_db.get_catalogs()
    elif isinstance(catalog_db, VOSDatabase):
        catalogs = catalog_db.get_catalogs()
    elif isinstance(catalog_db, (VOSCatalog, basestring)):
        catalogs = [(None, catalog_db)]
    elif isinstance(catalog_db, list):
        for x in catalog_db:
            assert (isinstance(x, (VOSCatalog, basestring)) and
                    not isinstance(x, VOSDatabase))
        catalogs = [(None, x) for x in catalog_db]
    else:  # pragma: no cover
        raise VOSError('catalog_db must be a catalog database, '
                       'a list of catalogs, or a catalog')

    if pedantic is None:  # pragma: no cover
        pedantic = VO_PEDANTIC

    for name, catalog in catalogs:
        if isinstance(catalog, basestring):
            if catalog.startswith("http"):
                url = catalog
            else:
                remote_db = get_remote_catalog_db(service_type, cache=cache)
                catalog = remote_db.get_catalog(catalog)
                url = catalog['url']
        else:
            url = catalog['url']

        if verbose:  # pragma: no cover
            color_print('Trying {0}'.format(url), 'green')

        try:
            return _vo_service_request(url, pedantic, kwargs)
        except Exception as e:
            vo_warn(W25, (url, str(e)))

    raise VOSError('None of the available catalogs returned valid results.')


def list_catalogs(service_type, cache=True, **kwargs):
    """
    List the catalogs available for the given service type.

    Parameters
    ----------
    service_type : {'conesearch_good', 'conesearch_warn'}
        Only Cone Search is supported for now.

    cache : bool
        See `get_remote_catalog_db`.

    kwargs : keywords for `VOSDatabase.list_catalogs`

    Returns
    -------
    arr : list of str

    """
    return get_remote_catalog_db(service_type,
                                 cache=cache).list_catalogs(**kwargs)
