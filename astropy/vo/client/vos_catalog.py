# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Common utilities for accessing VO simple services.

When run for the first time, catalogs database will be
downloaded to a local cache. See `astropy.config` for
caching behavior.

*CONFIGURABLE PROPERTIES*

These properties are set via Astropy configuration system:

    * `astropy.io.votable.pedantic`
    * `astropy.vo.client.vos_baseurl`
    * `astropy.vo.client.vos_timeout`

Examples
--------
>>> from astropy.vo.client import vos_catalog

Get all catalogs from a database named 'conesearch'.
For more cone search examples, see `astropy.vo.client.conesearch`:

>>> my_db = vos_catalog.get_remote_catalog_db('conesearch')

Find catalog names containing 'noao':

>>> cat_names = my_db.list_catalogs(match_string='noao', sort=True)

Get information for first catalog from above. Catalog
fields may vary:

>>> my_cat = my_db.get_catalog(cat_names[0])
>>> print my_cat
>>> cat_keys = my_cat.keys()
>>> cat_url = my_cat['url']
>>> max_rec = my_cat['maxRecords']

"""
from __future__ import print_function, division

# STDLIB
import io
import json

# THIRD PARTY
import numpy as np

# LOCAL
from ...io.votable import table
from ...io.votable.exceptions import vo_warn, W24, W25
from ...io.votable.util import IS_PY3K
from ...utils import webquery
from ...utils.console import color_print

# LOCAL CONFIG
from ...config.configuration import ConfigurationItem, get_config_items
from ...config.data import get_data_fileobj
from ...logger import log

__all__ = ['VOSCatalog', 'VOSDatabase',
           'get_remote_catalog_db', 'call_vo_service', 'list_catalogs']

__dbversion__ = 1

BASEURL = ConfigurationItem('vos_baseurl',
                            'http://stsdas.stsci.edu/astrolib/vo_databases/',
                            'URL where VO Service database file is stored.')

TIMEOUT = ConfigurationItem('vos_timeout', 30.0,
                            'Timeout in seconds for VO Service query')

# LOCAL CONFIG FROM DIFFERENT MODULE
VO_CFG = get_config_items('astropy.io.votable')
try:
    VO_PEDANTIC = VO_CFG['PEDANTIC']
except KeyError as e:
    # Could use set, but do not want to mess with general vo config
    log.warn('astropy.io.votable.pedantic not found in config, '
             'defaulting to True')
    VO_PEDANTIC = True


class VOSError(Exception):
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
        return json.dumps(self._tree, indent=4)

    def __getattr__(self, what):
        """Expose dictionary attributes."""
        return getattr(self._tree, what)

    def __getitem__(self, what):
        """Expose dictionary key look-up."""
        return self._tree[what]


class VOSDatabase(object):
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

        if tree['__version__'] > __dbversion__:
            vo_warn(W24)

        if not 'catalogs' in tree:
            raise VOSError("Invalid VO service catalog database")

        self._catalogs = tree['catalogs']

    def get_catalogs(self):
        """Iterator to get all catalogs."""
        for key, val in self._catalogs.items():
            yield key, VOSCatalog(val)

    def get_catalog(self, name):
        """
        Get one catalog of given name.

        Parameters
        ----------
        name : str
            Primary key identifying the catalog.

        Returns
        -------
        value : `VOSCatalog` object

        Raises
        ------
        VOSError
            If catalog is not found.

        """
        if not name in self._catalogs:
            raise VOSError("No catalog '{}' found.".format(name))
        return VOSCatalog(self._catalogs[name])

    def list_catalogs(self, match_string=None, sort=False):
        """
        List of catalog names.

        Parameters
        ----------
        match_string : str or `None`
            If given string is anywhere in a catalog name, it is
            considered a matching catalog. It is not case-sensitive.
            By default, all catalogs are returned.

        sort : bool
            Sort output in alphabetical order. If not sorted, the
            order depends on dictionary hashing.

        """
        all_catalogs = self._catalogs.keys()

        if match_string is None:
            out_arr = all_catalogs
        else:
            all_cat_arr = np.array(all_catalogs)
            all_cat_ucase = np.char.upper(all_cat_arr)
            i = np.char.count(all_cat_ucase,
                              match_string.upper()).astype('bool')
            out_arr = list(all_cat_arr[i])

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
    value : `VOSDatabase` object

    """
    fd = get_data_fileobj(BASEURL() + dbname + '.json', cache=cache)
    if IS_PY3K:
        wrapped_fd = io.TextIOWrapper(fd, 'utf8')
    else:
        wrapped_fd = fd
    try:
        tree = json.load(wrapped_fd)
    finally:
        fd.close()

    return VOSDatabase(tree)


def _vo_service_request(url, pedantic, kwargs):
    req = webquery.webget_open(url, timeout=TIMEOUT(), **kwargs)
    try:
        tab = table.parse(req, filename=req.geturl(), pedantic=pedantic)
    finally:
        req.close()

    # In case of errors from the server, a complete and correct
    # "stub" VOTable file may still be returned.  This is to
    # detect that case.
    for param in tab.iter_fields_and_params():
        if param.ID.lower() == 'error':
            raise VOSError("Catalog server '{}' returned error '{}'".format(
                url, param.value))
        break

    for info in tab.resources[0].infos:
        if info.name == 'QUERY_STATUS' and info.value != 'OK':
            if info.content is not None:
                long_descr = ':\n{}'.format(info.content)
            else:
                long_descr = ''
            raise VOSError("Catalog server '{}' returned status '{}'{}".format(
                url, info.value, long_descr))
        break

    out_tab = tab.get_first_table()
    if out_tab.array.size <= 0:
        raise VOSError("Catalog server '{}' returned {} result".format(
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
        if one is not provided.

    catalog_db : may be one of the following, in order from easiest to
    use to most control:

        - `None`: A database of conesearch catalogs is downloaded from
           STScI.  The first catalog in the database to successfully
           return a result is used.

        - *catalog name*: A name in the database of conesearch catalogs
          at STScI is used.  For a list of acceptable names, see
          :func:`list_catalogs`.

        - *url*: The prefix of a *url* to a IVOA Cone Search Service.
          Must end in either `?` or `&`.

        - A `~vos_catalog.VOSCatalog` instance: A specific catalog
          manually downloaded and selected from the database using the
          APIs in the :mod:`vos_catalog` module.

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
        the service etc.

    Returns
    -------
    value : `astropy.io.votable.tree.Table` object
        First table from first successful VO service request.

    Raises
    ------
    VOSError
        If VO service request fails.

    """
    if catalog_db is None:
        catalog_db = get_remote_catalog_db(service_type, cache=cache)
        catalogs = catalog_db.get_catalogs()
    elif isinstance(catalog_db, (VOSCatalog, basestring)):
        catalogs = [(None, catalog_db)]
    elif isinstance(catalog_db, list):
        for x in catalog_db:
            assert isinstance(
                x, (VOSCatalog, basestring))
        catalogs = [(None, x) for x in catalog_db]
    elif isinstance(catalog_db, VOSDatabase):
        catalogs = catalog_db.get_catalogs()
    else:
        raise VOSError('catalog_db must be a catalog database, '
                       'a list of catalogs, or a catalog')

    if pedantic is None:
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

        if verbose:
            color_print('Trying {}'.format(url), 'green')

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
    service_type : {'basic', 'conesearch', 'image', 'ssa'}
        Only 'conesearch' is supported for now.
        Others are for testing only.

    cache : bool
        See `get_remote_catalog_db`.

    kwargs : keywords for `VOSDatabase.list_catalogs`

    Returns
    -------
    value : list of str

    """
    return get_remote_catalog_db(service_type,
                                 cache=cache).list_catalogs(**kwargs)
