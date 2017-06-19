# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Common utilities for accessing VO simple services."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six
from ...extern.six.moves import urllib

# STDLIB
import fnmatch
import json
import os
import re
import socket
import warnings

from collections import defaultdict
from copy import deepcopy

# LOCAL
from .exceptions import (VOSError, MissingCatalog, DuplicateCatalogName,
                         DuplicateCatalogURL, InvalidAccessURL)
from .. import conf as vo_conf
from ...io.votable import parse_single_table, table, tree
from ...io.votable import conf as votable_conf
from ...io.votable.exceptions import vo_raise, vo_warn, E19, W24, W25
from ...utils.console import color_print
from ...utils.data import get_readable_fileobj
from ...utils.data import conf as data_conf
from ...utils.exceptions import AstropyUserWarning
from ...utils.misc import JsonCustomEncoder
from ...utils.xml.unescaper import unescape_all
from ...utils.decorators import deprecated, deprecated_renamed_argument


__all__ = ['VOSBase', 'VOSCatalog', 'VOSDatabase', 'get_remote_catalog_db',
           'call_vo_service', 'list_catalogs']

__dbversion__ = 1


@deprecated(
    '2.0', alternative='astroquery.vo_conesearch.vos_catalog.VOSBase')
class VOSBase(object):
    """Base class for `VOSCatalog` and `VOSDatabase`.

    Parameters
    ----------
    tree : JSON tree

    """
    def __init__(self, tree):
        self._tree = tree

    def __getattr__(self, what):
        """Expose dictionary attributes."""
        return getattr(self._tree, what)

    def __getitem__(self, what):
        """Expose dictionary key look-up."""
        return self._tree[what]

    def __setitem__(self, what, value):
        """Expose dictionary key assignment."""
        self._tree[what] = value

    def __iter__(self):
        """Expose dictionary iteration."""
        return iter(self._tree)

    def dumps(self):
        """Dump the contents into a string.

        Returns
        -------
        s : str
            Contents as JSON string dump.

        """
        return json.dumps(self._tree, cls=JsonCustomEncoder, sort_keys=True,
                          indent=4)


@deprecated(
    '2.0', alternative='astroquery.vo_conesearch.vos_catalog.VOSCatalog')
class VOSCatalog(VOSBase):
    """A class to represent VO Service Catalog.

    Parameters
    ----------
    tree : JSON tree

    Raises
    ------
    VOSError
        Missing necessary key(s).

    """
    _compulsory_keys = ('title', 'url')

    def __init__(self, tree):
        super(VOSCatalog, self).__init__(tree)

        for key in self._compulsory_keys:
            if key not in self._tree:
                raise VOSError('Catalog must have "{0}" key.'.format(key))

    def __str__(self):  # pragma: no cover
        """Show the most important and unique things about a catalog."""
        out_str = '\n'.join(['{0}: {1}'.format(key, self._tree[key])
                             for key in self._compulsory_keys
                             if key in self._tree])
        return out_str

    def delete_attribute(self, key):
        """Delete given metadata key and its value from the catalog.

        Parameters
        ----------
        key : str
            Metadata key to delete.

        Raises
        ------
        KeyError
            Key not found.

        VOSError
            Key must exist in catalog, therefore cannot be deleted.

        """
        if key in self._compulsory_keys:
            raise VOSError('{0} must exist in catalog, therefore cannot be '
                           'deleted.'.format(key))
        del self._tree[key]

    @classmethod
    def create(cls, title, url, **kwargs):
        """Create a new VO Service Catalog with user parameters.

        Parameters
        ----------
        title : str
            Title of the catalog.

        url : str
            Access URL of the service. This is used to build queries.

        kwargs : dict
            Additional metadata as keyword-value pairs describing the catalog,
            except 'title' and 'url'.

        Returns
        -------
        cat : `VOSCatalog`
            VO Service Catalog.

        Raises
        ------
        TypeError
            Multiple values given for keyword argument.

        """
        tree = {'title': title, 'url': url}
        tree.update(kwargs)
        return cls(tree)


@deprecated(
    '2.0', alternative='astroquery.vo_conesearch.vos_catalog.VOSDatabase')
class VOSDatabase(VOSBase):
    """A class to represent a collection of `VOSCatalog`.

    Parameters
    ----------
    tree : JSON tree

    Raises
    ------
    VOSError
        If given ``tree`` does not have 'catalogs' key
        or catalog is invalid.

    """
    def __init__(self, tree):
        if 'catalogs' not in tree:
            raise VOSError("Invalid VO service catalog database")

        super(VOSDatabase, self).__init__(tree)
        self._catalogs = tree['catalogs']

        if self.version > __dbversion__:  # pragma: no cover
            vo_warn(W24)

        # Maps access URL to primary key(s).
        # URL is the real key, but we chose title because it is more readable
        # when written out to JSON.
        self._url_keys = defaultdict(list)
        for key, cat in self.get_catalogs():
            self._url_keys[cat['url']].append(key)

    def __str__(self):  # pragma: no cover
        """Show the most important and unique things about a database."""
        return '\n'.join(sorted(self._catalogs))

    def __len__(self):
        """Return the number of catalogs in database."""
        return len(self._catalogs)

    @property
    def version(self):
        """Database version number."""
        return self._tree['__version__']

    def get_catalogs(self):
        """Iterator to get all catalogs."""
        for key, val in self._catalogs.items():
            yield key, VOSCatalog(val)

    def get_catalogs_by_url(self, url):
        """Like :func:`get_catalogs` but using access URL look-up."""
        keys = self._url_keys[url]
        for key in keys:
            yield key, VOSCatalog(self._catalogs[key])

    def get_catalog(self, name):
        """Get one catalog of given name.

        Parameters
        ----------
        name : str
            Primary key identifying the catalog.

        Returns
        -------
        obj : `VOSCatalog`

        Raises
        ------
        MissingCatalog
            If catalog is not found.

        """
        if name not in self._catalogs:
            raise MissingCatalog("No catalog '{0}' found.".format(name))

        return VOSCatalog(self._catalogs[name])

    def get_catalog_by_url(self, url):
        """Like :func:`get_catalog` but using access URL look-up.
        On multiple matches, only first match is returned.

        """
        keys = self._url_keys[url]

        if len(keys) < 1:
            raise MissingCatalog("No catalog with URL '{0}' found.".format(url))

        return VOSCatalog(self._catalogs[keys[0]])

    @staticmethod
    def _match_pattern(all_keys, pattern, sort):
        """Used by :func:`list_catalogs` and :func:`list_catalogs_by_url`."""
        if pattern is None or len(all_keys) == 0:
            out_arr = all_keys
        else:
            pattern = re.compile(fnmatch.translate('*' + pattern + '*'),
                                 re.IGNORECASE)
            out_arr = [s for s in all_keys if pattern.match(s)]

        if sort:
            out_arr.sort()

        return out_arr

    def list_catalogs(self, pattern=None, sort=True):
        """List catalog names.

        Parameters
        ----------
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
        out_arr : list of str
            List of catalog names.

        """
        return self._match_pattern(list(self._catalogs), pattern, sort)

    def list_catalogs_by_url(self, pattern=None, sort=True):
        """Like :func:`list_catalogs` but using access URL."""
        out_arr = self._match_pattern(list(self._url_keys), pattern, sort)

        # Discard URL that maps to nothing
        return [a for a in out_arr if len(self._url_keys[a]) > 0]

    def add_catalog(self, name, cat, allow_duplicate_url=False):
        """Add a catalog to database.

        Parameters
        ----------
        name : str
            Primary key for the catalog.

        cat : `VOSCatalog`
            Catalog to add.

        allow_duplicate_url : bool
            Allow catalog with duplicate access URL?

        Raises
        ------
        VOSError
            Invalid catalog.

        DuplicateCatalogName
            Catalog with given name already exists.

        DuplicateCatalogURL
            Catalog with given access URL already exists.

        """
        if not isinstance(cat, VOSCatalog):
            raise VOSError('{0} is not a VO Service Catalog.'.format(cat))

        if name in self._catalogs:
            raise DuplicateCatalogName('{0} already exists.'.format(name))

        url = cat['url']
        names = self._url_keys[url]
        if len(names) > 0 and not allow_duplicate_url:
            raise DuplicateCatalogURL(
                '{0} already exists: {1}'.format(url, names))

        self._catalogs[name] = deepcopy(cat._tree)
        self._url_keys[url].append(name)

    def add_catalog_by_url(self, name, url, **kwargs):
        """Like :func:`add_catalog` but the catalog is created with
        only the given name and access URL.

        Parameters
        ----------
        name : str
            Primary key for the catalog.

        url : str
            Access URL of the service. This is used to build queries.

        kwargs : dict
            Keywords accepted by :func:`add_catalog`.

        """
        self.add_catalog(name, VOSCatalog.create(name, url), **kwargs)

    def delete_catalog(self, name):
        """Delete a catalog from database with given name.

        Parameters
        ----------
        name : str
            Primary key identifying the catalog.

        Raises
        ------
        MissingCatalog
            If catalog is not found.

        """
        if name not in self._catalogs:
            raise MissingCatalog('{0} not found.'.format(name))

        self._url_keys[self._catalogs[name]['url']].remove(name)
        del self._catalogs[name]

    def delete_catalog_by_url(self, url):
        """Like :func:`delete_catalog` but using access URL.
        On multiple matches, all matches are deleted.

        """
        keys = sorted(self._url_keys[url])  # Makes a copy of list

        if len(keys) < 1:
            raise MissingCatalog('{0} not found.'.format(url))

        for key in keys:
            self.delete_catalog(key)

    def merge(self, other, **kwargs):
        """Merge two database together.

        Parameters
        ----------
        other : `VOSDatabase`
            The other database to merge.

        kwargs : dict
            Keywords accepted by :func:`add_catalog`.

        Returns
        -------
        db : `VOSDatabase`
            Merged database.

        Raises
        ------
        VOSError
            Invalid database or incompatible version.

        """
        if not isinstance(other, VOSDatabase):
            raise VOSError('{0} is not a VO database.'.format(other))

        if other.version != self.version:
            raise VOSError('Incompatible database version: {0}, '
                           '{1}'.format(self.version, other.version))

        db = VOSDatabase.create_empty()
        for old_db in (self, other):
            for key, cat in old_db.get_catalogs():
                db.add_catalog(key, cat, **kwargs)

        return db

    @deprecated_renamed_argument('clobber', 'overwrite', '1.3')
    def to_json(self, filename, overwrite=False):
        """
        Write database content to a JSON file.

        Parameters
        ----------
        filename : str
            JSON file.

        overwrite : bool
            If ``True``, overwrite the output file if it exists.

            .. versionchanged:: 1.3
               ``overwrite`` replaces the deprecated ``clobber`` argument.

        Raises
        ------
        OSError
            If the file exists and ``overwrite`` is ``False``.
        """
        if os.path.exists(filename) and not overwrite:
            raise OSError('{0} exists.'.format(filename))

        with open(filename, 'w') as fd:
            fd.write(self.dumps())

    @classmethod
    def create_empty(cls):
        """Create an empty database of VO services.

        Empty database format::

            {
                "__version__": 1,
                "catalogs" : {
                }
            }

        Returns
        -------
        db : `VOSDatabase`
            Empty database.

        """
        return cls({'__version__': __dbversion__, 'catalogs': {}})

    @classmethod
    def from_json(cls, filename, **kwargs):
        """Create a database of VO services from a JSON file.

        Example JSON format for Cone Search::

            {
                "__version__": 1,
                "catalogs" : {
                    "My Cone Search": {
                        "capabilityClass": "ConeSearch",
                        "title": "My Cone Search",
                        "url": "http://foo/cgi-bin/search?CAT=bar&",
                        ...
                    },
                    "Another Cone Search": {
                        ...
                    }
                }
            }

        Parameters
        ----------
        filename : str
            JSON file.

        kwargs : dict
            Keywords accepted by
            :func:`~astropy.utils.data.get_readable_fileobj`.

        Returns
        -------
        db : `VOSDatabase`
            Database from given file.

        """
        with get_readable_fileobj(filename, **kwargs) as fd:
            tree = json.load(fd)

        return cls(tree)

    @classmethod
    def from_registry(cls, registry_url, timeout=60, **kwargs):
        """Create a database of VO services from VO registry URL.

        This is described in detail in :ref:`vo-sec-validator-build-db`,
        except for the ``validate_xxx`` keys that are added by the
        validator itself.

        Parameters
        ----------
        registry_url : str
            URL of VO registry that returns a VO Table.
            For example, see ``astropy.vo.validator.validate.CS_MSTR_LIST``.
            Pedantic is automatically set to `False` for parsing.

        timeout : number
            Temporarily set `astropy.utils.data.Conf.remote_timeout` to
            this value to avoid time out error while reading the
            entire registry.

        kwargs : dict
            Keywords accepted by
            :func:`~astropy.utils.data.get_readable_fileobj`.

        Returns
        -------
        db : `VOSDatabase`
            Database from given registry.

        Raises
        ------
        VOSError
            Invalid VO registry.

        """
        # Download registry as VO table
        with data_conf.set_temp('remote_timeout', timeout):
            with get_readable_fileobj(registry_url, **kwargs) as fd:
                tab_all = parse_single_table(fd, pedantic=False)

        # Registry must have these fields
        compulsory_fields = ['title', 'accessURL']
        cat_fields = tab_all.array.dtype.names
        for field in compulsory_fields:
            if field not in cat_fields:  # pragma: no cover
                raise VOSError('"{0}" is missing from registry.'.format(field))

        title_counter = defaultdict(int)
        title_fmt = '{0} {1}'
        db = cls.create_empty()

        # Each row in the table becomes a catalog
        for arr in tab_all.array.data:
            cur_cat = {}
            cur_key = ''

            # Process each field and build the catalog.
            # Catalog is completely built before being thrown out
            # because codes need less changes should we decide to
            # allow duplicate URLs in the future.
            for field in cat_fields:

                # For primary key, a number needs to be appended to the title
                # because registry can have multiple entries with the same
                # title but different URLs.
                if field == 'title':
                    cur_title = arr['title']
                    title_counter[cur_title] += 1  # Starts with 1

                    if isinstance(cur_title, bytes):  # pragma: py3
                        cur_key = title_fmt.format(cur_title.decode('utf-8'),
                                                   title_counter[cur_title])
                    else:  # pragma: py2
                        cur_key = title_fmt.format(cur_title,
                                                   title_counter[cur_title])

                # Special handling of access URL, otherwise no change.
                if field == 'accessURL':
                    cur_cat['url'] = unescape_all(arr['accessURL'])
                else:
                    cur_cat[field] = arr[field]

            # New field to track duplicate access URLs.
            cur_cat['duplicatesIgnored'] = 0

            # Add catalog to database, unless duplicate access URL exists.
            # In that case, the entry is thrown out and the associated
            # counter is updated.
            dup_keys = db._url_keys[cur_cat['url']]
            if len(dup_keys) < 1:
                db.add_catalog(
                    cur_key, VOSCatalog(cur_cat), allow_duplicate_url=False)
            else:
                db._catalogs[dup_keys[0]]['duplicatesIgnored'] += 1
                warnings.warn(
                    '{0} is thrown out because it has same access URL as '
                    '{1}.'.format(cur_key, dup_keys[0]), AstropyUserWarning)

        return db


@deprecated(
    '2.0',
    alternative='astroquery.vo_conesearch.vos_catalog.get_remote_catalog_db')
def get_remote_catalog_db(dbname, cache=True, verbose=True):
    """Get a database of VO services (which is a JSON file) from a remote
    location.

    Parameters
    ----------
    dbname : str
        Prefix of JSON file to download from
        `astropy.vo.Conf.vos_baseurl`.

    cache : bool
        Use caching for VO Service database. Access to actual VO
        websites referenced by the database still needs internet
        connection.

    verbose : bool
        Show download progress bars.

    Returns
    -------
    db : `VOSDatabase`
        A database of VO services.

    """
    return VOSDatabase.from_json(
        urllib.parse.urljoin(vo_conf.vos_baseurl, dbname + '.json'),
        encoding='utf8', cache=cache,
        show_progress=verbose)


def _get_catalogs(service_type, catalog_db, **kwargs):
    """Expand ``catalog_db`` to a list of catalogs.

    Parameters
    ----------
    service_type, catalog_db
        See :func:`call_vo_service`.

    kwargs : dict
        Keywords accepted by :func:`get_remote_catalog_db`.

    Returns
    -------
    catalogs : list of tuple
        List of catalogs in the form of ``(key, VOSCatalog)``.

    Raises
    ------
    VOSError
        Invalid ``catalog_db``.

    """
    if catalog_db is None:
        catalog_db = get_remote_catalog_db(service_type, **kwargs)
        catalogs = catalog_db.get_catalogs()
    elif isinstance(catalog_db, VOSDatabase):
        catalogs = catalog_db.get_catalogs()
    elif isinstance(catalog_db, (VOSCatalog, six.string_types)):
        catalogs = [(None, catalog_db)]
    elif isinstance(catalog_db, list):
        for x in catalog_db:
            if (not isinstance(x, (VOSCatalog, six.string_types)) or
                    isinstance(x, VOSDatabase)):
                raise VOSError('catalog_db must be a catalog database, '
                       'a list of catalogs, or a catalog')
        catalogs = [(None, x) for x in catalog_db]
    else:  # pragma: no cover
        raise VOSError('catalog_db must be a catalog database, '
                       'a list of catalogs, or a catalog')

    return catalogs


def _vo_service_request(url, pedantic, kwargs, cache=True, verbose=False):
    """This is called by :func:`call_vo_service`.

    Raises
    ------
    InvalidAccessURL
        Invalid access URL.

    """
    if len(kwargs) and not url.endswith(('?', '&')):
        raise InvalidAccessURL("url should already end with '?' or '&'")

    query = []
    for key, value in six.iteritems(kwargs):
        query.append('{0}={1}'.format(
            urllib.parse.quote(key), urllib.parse.quote_plus(str(value))))

    parsed_url = url + '&'.join(query)
    with get_readable_fileobj(parsed_url, encoding='binary', cache=cache,
                              show_progress=verbose) as req:
        tab = table.parse(req, filename=parsed_url, pedantic=pedantic)

    return vo_tab_parse(tab, url, kwargs)


@deprecated(
    '2.0', alternative='astroquery.vo_conesearch.vos_catalog.vo_tab_parse')
def vo_tab_parse(tab, url, kwargs):
    """In case of errors from the server, a complete and correct
    'stub' VOTable file may still be returned.  This is to
    detect that case.

    Parameters
    ----------
    tab : `astropy.io.votable.tree.VOTableFile`

    url : str
        URL used to obtain ``tab``.

    kwargs : dict
        Keywords used to obtain ``tab``, if any.

    Returns
    -------
    out_tab : `astropy.io.votable.tree.Table`

    Raises
    ------
    IndexError
        Table iterator fails.

    VOSError
        Server returns error message or invalid table.

    """
    for param in tab.iter_fields_and_params():
        if param.ID is not None and param.ID.lower() == 'error':
            if isinstance(param, tree.Param):
                e = param.value
            else:  # pragma: no cover
                e = ''
            raise VOSError("Catalog server '{0}' returned error '{1}'".format(
                url, e))

    for info in tab.infos:
        if info.name is not None and info.name.lower() == 'error':
            raise VOSError("Catalog server '{0}' returned error '{1}'".format(
                url, info.value))

    if tab.resources == []:  # pragma: no cover
        vo_raise(E19)

    for info in tab.resources[0].infos:
        if ((info.name == 'QUERY_STATUS' and info.value != 'OK') or
                (info.name is not None and info.name.lower() == 'error')):
            if info.content is not None:  # pragma: no cover
                long_descr = ':\n{0}'.format(info.content)
            else:
                long_descr = ''
            raise VOSError("Catalog server '{0}' returned status "
                           "'{1}'{2}".format(url, info.value, long_descr))

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


@deprecated(
    '2.0', alternative='astroquery.vo_conesearch.vos_catalog.call_vo_service')
def call_vo_service(service_type, catalog_db=None, pedantic=None,
                    verbose=True, cache=True, kwargs={}):
    """Makes a generic VO service call.

    Parameters
    ----------
    service_type : str
        Name of the type of service, e.g., 'conesearch_good'.
        Used in error messages and to select a catalog database
        if ``catalog_db`` is not provided.

    catalog_db
        May be one of the following, in order from easiest to
        use to most control:

            - `None`: A database of ``service_type`` catalogs is
              downloaded from `astropy.vo.Conf.vos_baseurl`.  The
              first catalog in the database to successfully return a
              result is used.

            - *catalog name*: A name in the database of
              ``service_type`` catalogs at
              `astropy.vo.Conf.vos_baseurl` is used.  For a list of
              acceptable names, use :func:`list_catalogs`.

            - *url*: The prefix of a URL to a IVOA Service for
              ``service_type``. Must end in either '?' or '&'.

            - `VOSCatalog` object: A specific catalog manually downloaded and
              selected from the database (see :ref:`vo-sec-client-vos`).

            - Any of the above 3 options combined in a list, in which case
              they are tried in order.

    pedantic : bool or `None`
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        :py:mod:`warnings` module.  When not provided, uses the
        configuration setting `astropy.io.votable.Conf.pedantic`, which
        defaults to `False`.

    verbose : bool
        Verbose output.

    cache : bool
        Use caching for VO Service database. Access to actual VO
        websites referenced by the database still needs internet
        connection.

    kwargs : dictionary
        Keyword arguments to pass to the catalog service.
        No checking is done that the arguments are accepted by
        the service, etc.

    Returns
    -------
    obj : `astropy.io.votable.tree.Table`
        First table from first successful VO service request.

    Raises
    ------
    VOSError
        If VO service request fails.

    """
    n_timed_out = 0
    catalogs = _get_catalogs(service_type, catalog_db, cache=cache,
                             verbose=verbose)

    if pedantic is None:  # pragma: no cover
        pedantic = votable_conf.pedantic

    for name, catalog in catalogs:
        if isinstance(catalog, six.string_types):
            if catalog.startswith("http"):
                url = catalog
            else:
                remote_db = get_remote_catalog_db(service_type, cache=cache,
                                                  verbose=verbose)
                catalog = remote_db.get_catalog(catalog)
                url = catalog['url']
        else:
            url = catalog['url']

        if verbose:  # pragma: no cover
            color_print('Trying {0}'.format(url), 'green')

        try:
            return _vo_service_request(url, pedantic, kwargs, cache=cache,
                                       verbose=verbose)
        except Exception as e:
            vo_warn(W25, (url, str(e)))
            if hasattr(e, 'reason') and isinstance(e.reason, socket.timeout):
                n_timed_out += 1

    err_msg = 'None of the available catalogs returned valid results.'
    if n_timed_out > 0:
        err_msg += ' ({0} URL(s) timed out.)'.format(n_timed_out)
    raise VOSError(err_msg)


@deprecated(
    '2.0', alternative='astroquery.vo_conesearch.vos_catalog.list_catalogs')
def list_catalogs(service_type, cache=True, verbose=True, **kwargs):
    """List the catalogs available for the given service type.

    Parameters
    ----------
    service_type : str
        Name of the type of service, e.g., 'conesearch_good'.

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
    return get_remote_catalog_db(service_type, cache=cache,
                                 verbose=verbose).list_catalogs(**kwargs)
