# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Inspect results from `astropy.vo.validator.validate`."""
from __future__ import absolute_import, division, print_function, unicode_literals

# STDLIB
import sys

# LOCAL
from ..client.vos_catalog import get_remote_catalog_db


__all__ = ['ConeSearchResults']


class ConeSearchResults(object):
    """A class to store Cone Search validation results.

    Attributes
    ----------
    dbtypes : list
        Cone Search database identifiers.

    dbs : dict
        Stores `~astropy.vo.client.vos_catalog.VOSDatabase`
        for each ``dbtypes``.

    catkeys : dict
        Stores sorted catalog keys for each ``dbtypes``.

    Parameters
    ----------
    cache : bool
       Read from cache, if available.
       Default is `False` to ensure the latest data are read.

    verbose : bool
        Show download progress bars.

    """
    def __init__(self, cache=False, verbose=True):
        self.dbtypes = ['good', 'warn', 'exception', 'error']
        self.dbs = {}
        self.catkeys = {}

        for typ in self.dbtypes:
            self.dbs[typ] = get_remote_catalog_db(
                'conesearch_' + typ, cache=cache, verbose=verbose)
            self.catkeys[typ] = self.dbs[typ].list_catalogs()

    def tally(self, fout=None):
        """
        Tally databases.

        Parameters
        ----------
        fout : output stream
            Default is screen output.

        """
        if fout is None:  # pragma: no cover
            fout = sys.stdout

        str_list = []
        n_tot = 0

        for typ in self.dbtypes:
            n_cur = len(self.catkeys[typ])
            n_tot += n_cur
            str_list.append('{0}: {1} catalog(s)'.format(typ, n_cur))

        if len(str_list) > 0:
            str_list.append('total: {0} catalog(s)\n'.format(n_tot))
            fout.write('\n'.join(str_list))

    def list_cats(self, typ, fout=None, ignore_noncrit=False):
        """
        List catalogs in given database.

        Listing contains:

            #. Catalog key
            #. Cone search access URL
            #. Warning codes
            #. Warning descriptions

        Parameters
        ----------
        typ : str
            Any value in ``self.dbtypes``.

        fout : output stream
            Default is screen output.

        ignore_noncrit : bool
            Exclude warnings in
            `astropy.vo.validator.Conf.noncritical_warnings`.  This is
            useful to see why a catalog failed validation.

        """
        if fout is None:  # pragma: no cover
            fout = sys.stdout

        assert typ in self.dbtypes
        str_list = []

        for cat in self.catkeys[typ]:
            cat_db = self.dbs[typ].get_catalog(cat)

            if ignore_noncrit:
                out_wt = _exclude_noncrit(cat_db['validate_warning_types'])
                out_ws = _exclude_noncrit(cat_db['validate_warnings'])
            else:
                out_wt = cat_db['validate_warning_types']
                out_ws = cat_db['validate_warnings']

            # Warning types contains None if some other Exception was thrown.
            # There should be only 1 occurence for each warning type.
            # But will put in a loop anyway, just in case.
            while None in out_wt:  # pragma: no cover
                out_wt[out_wt.index(None)] = 'None'

            str_list += [cat, cat_db['url']]
            if len(out_wt) > 0:
                str_list.append(','.join(out_wt))
            if len(out_ws) > 0:
                str_list.append('\n'.join(out_ws))
            str_list[-1] += '\n'

        if len(str_list) > 0:
            fout.write('\n'.join(str_list))

    def print_cat(self, key, fout=None):
        """
        Display a single catalog of given key.

        If not found, nothing is written out.

        Parameters
        -----------
        key : str
            Catalog key.

        fout : output stream
            Default is screen output.

        """
        if fout is None:  # pragma: no cover
            fout = sys.stdout

        str_list = []

        for typ in self.dbtypes:
            if key in self.catkeys[typ]:
                str_list += [self.dbs[typ].get_catalog(key).dumps(),
                             '\nFound in {0}'.format(typ)]

                # Only has one match, so quits when it is found
                break

        if len(str_list) > 0:
            fout.write('\n'.join(str_list) + '\n')


def _exclude_noncrit(in_list):
    """
    Exclude any items in input list containing
    `astropy.vo.validator.Conf.noncritical_warnings`.

    Parameters
    ----------
    in_list : list
        List of strings to process.

    Returns
    -------
    out_list : list
        List with only qualified strings.

    """
    from . import conf
    out_list = []
    for s in in_list:
        n = 0
        if s is not None:
            for w in conf.noncritical_warnings:
                n += s.count(w)
        if n == 0:  # pragma: no cover
            out_list.append(s)
    return out_list
