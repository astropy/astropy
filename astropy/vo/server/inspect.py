# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Inspect results from `astropy.vo.server.validate`.

Examples
--------
Load Cone Search validation results directly from
`astropy.vo.client.vos_baseurl`:

>>> from astropy.vo.server import inspect
>>> r = inspect.ConeSearchResults()

Print tally:

>>> r.tally()
good: 15 catalog(s)
warn: 15 catalog(s)
exception: 0 catalog(s)
error: 1 catalog(s)
total: 31 catalog(s)

Print a list of good Cone Search catalogs:

>>> r.list_cats('good')
2MASS All-Sky Catalog of Point Sources (Cutri+ 2003) 1
http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=II/246/out&
W22,W03
# ...

List Cone Search catalogs with warnings,
excluding warnings that were ignored in
`astropy.vo.server.noncrit_warnings`,
and writes the output to a file:

>>> with open('warn_cats.txt', 'w') as fout:
...     r.list_cats('warn', fout=fout, ignore_noncrit=True)

Print the first good Cone Search catalog:

>>> catkey = r.catkeys['good'][0]
>>> r.print_cat(catkey)
{
    \"validate_network_error\": null,
    \"capabilityClass\": \"ConeSearch\",
    \"updated\": \"2011-09-14T20:17:37\",
    \"capabilityValidationLevel\": \"\",
    # ...
    \"validate_xmllint\": true
}
Found in good

Load Cone Search validation results from local
validation in the current directory:

>>> import os
>>> from astropy.vo.client.vos_catalog import BASEURL
>>> BASEURL.set(os.curdir + os.sep)
>>> r = inspect.ConeSearchResults()

"""
from __future__ import print_function, division

# STDLIB
import sys

# THIRD PARTY
import numpy as np

# LOCAL
from ..client.vos_catalog import get_remote_catalog_db


__all__ = ['ConeSearchResults']


class ConeSearchResults(object):
    """
    A class to store Cone Search validation results.

    Attributes
    ----------
    dbtypes : list
        Conesearch database identifiers.

    dbs : dict
        Stores `astropy.vo.client.vos_catalog.VOSDatabase`
        for each `dbtypes`.

    catkeys : dict
        Stores sorted catalog keys for each `dbtypes`.

    Parameters
    ----------
    cache : bool
       Read from cache, if available.
       Default is `False` to ensure the latest data are read.

    """
    def __init__(self, cache=False):
        self.dbtypes = ['good', 'warn', 'exception', 'error']
        self.dbs = {}
        self.catkeys = {}

        for typ in self.dbtypes:
            self.dbs[typ] = get_remote_catalog_db(
                'conesearch_' + typ, cache=cache)
            self.catkeys[typ] = self.dbs[typ].list_catalogs(sort=True)

    def tally(self, fout=sys.stdout):
        """
        Tally databases.

        Parameters
        ----------
        fout : output stream
            Default is screen output.

        """
        str_list = []
        n_tot = 0

        for typ in self.dbtypes:
            n_cur = len(self.catkeys[typ])
            n_tot += n_cur
            str_list.append('{0}: {1} catalog(s)'.format(typ, n_cur))

        if len(str_list) > 0:
            str_list.append('total: {0} catalog(s)\n'.format(n_tot))
            fout.write('\n'.join(str_list))

    def list_cats(self, typ, fout=sys.stdout, ignore_noncrit=False):
        """
        List catalogs in given database.

        Listing contains:

            #. Catalog key
            #. Cone search access URL
            #. Warning codes
            #. Warning descriptions

        Parameters
        ----------
        typ : string
            Any value in `dbtypes`.

        fout : output stream
            Default is screen output.

        ignore_noncrit : bool
            Exclude warnings in `astropy.vo.server.noncrit_warnings`.
            This is useful to see why a catalog failed validation.

        """
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
            while None in out_wt:
                out_wt[out_wt.index(None)] = 'None'

            str_list += [cat, cat_db['url']]
            if len(out_wt) > 0:
                str_list.append(','.join(out_wt))
            if len(out_ws) > 0:
                str_list.append('\n'.join(out_ws))
            str_list[-1] += '\n'

        if len(str_list) > 0:
            fout.write('\n'.join(str_list))

    def print_cat(self, key, fout=sys.stdout):
        """
        Display a single catalog of given key.

        If not found, nothing is written out.

        Parameters
        -----------
        key : string
            Catalog key.

        fout : output stream
            Default is screen output.

        """
        str_list = []

        for typ in self.dbtypes:
            if key in self.catkeys[typ]:
                str_list += [repr(self.dbs[typ].get_catalog(key)),
                             '\nFound in {0}'.format(typ)]

                # Only has one match, so quits when it is found
                break

        if len(str_list) > 0:
            fout.write('\n'.join(str_list) + '\n')


def _exclude_noncrit(in_list):
    """
    Exclude any items in input list containing
    `astropy.vo.server.crit_warnings`.

    Parameters
    ----------
    in_list : list
        List of strings to process.

    Returns
    -------
    out_list : list
        List with only qualified strings.

    """
    from .validate import NONCRIT_WARNINGS
    out_list = []
    for s in in_list:
        n = 0
        if s is not None:
            for w in NONCRIT_WARNINGS():
                n += s.count(w)
        if n == 0:
            out_list.append(s)
    return out_list
