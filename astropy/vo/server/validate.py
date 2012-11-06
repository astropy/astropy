# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Validate VO Services.

*CONFIGURABLE PROPERTIES*

These properties are set via Astropy configuration system:

    * `astropy.vo.server.cs_mstr_list`
    * Also depends on properties in `astropy.vo.client`

.. note::

    This is not meant to be used by a typical AstroPy user.

"""
# STDLIB
from copy import deepcopy
import json
import os
import shutil
import time
import urllib2

# THIRD PARTY
import numpy

# LOCAL
from ..client import vos_catalog
from ...io import votable
from ...io.votable.validator import result
from ...logger import log
from ...utils.misc import NumpyScalarEncoder

# LOCAL CONFIG
from ...config.configuration import ConfigurationItem

__all__ = ['check_conesearch_sites']

CS_MSTR_LIST = ConfigurationItem(
    'cs_mstr_list',
    'http://vao.stsci.edu/directory/NVORegInt.asmx/VOTCapabilityPredOpt?'
    'predicate=1%3D1&capability=conesearch&VOTStyleOption=2',
    'Conesearch master list query from VAO.')

_TMP_ROOT = None  # Set by `check_conesearch_sites`


def check_conesearch_sites(destdir=os.curdir, verbose=True, multiproc=True):
    """
    Validate Cone Search Services.

    A master list of all available Cone Search sites is
    obtained from `astropy.vo.server.cs_mstr_list`, which
    is a URL query to an external VAO service.

    These sites are validated using `astropy.io.votable.validator`
    and separated into four groups below, each is stored as a JSON
    database in `destdir`. Existing files with same names will be
    deleted to avoid confusion:

        #. 'conesearch_good.json' - Passed validation without any
           warnings or exceptions. If this database is not empty,
           it will be used by `astropy.vo.client.conesearch`.
        #. 'conesearch_warn.json' - Has some warnings but no
           exceptions. Use this instead if database above is empty.
        #. 'conesearch_exception.json' - Has some exceptions.
           *Never* use this. For informational purpose only.
        #. 'conesearch_error.json' - Has network connection error.
           *Never* use this. For informational purpose only.

    For convenience, the database that `astropy.vo.client.conesearch`
    is supposed to use will be symbolically linked to
    'conesearch.json' (Unix only).

    .. note::

        At the time this was written, 'conesearch_good.json' did
        not yield any results.

    Parameters
    ----------
    destdir : string
        Directory to store output files. Will be created if does
        not exist.

    verbose : bool
        Print extra info to log.

    multiproc : bool
        Enable multiprocessing.

    Raises
    ------
    AssertionError
        Parameter failed assertion test.

    timeout
        URL request timed out.

    """
    global _TMP_ROOT

    assert (not os.path.exists(destdir) and len(destdir) > 0) or \
        (os.path.exists(destdir) and os.path.isdir(destdir)), \
        'Invalid destination directory'

    if not os.path.exists(destdir):
        os.mkdir(destdir)

    if destdir[-1] != os.sep:
        destdir += os.sep

    # Temp dir created by votable.validator
    _TMP_ROOT = destdir + 'results'

    # Output files
    db_file = {}
    db_file['good'] = destdir + 'conesearch_good.json'
    db_file['warn'] = destdir + 'conesearch_warn.json'
    db_file['excp'] = destdir + 'conesearch_exception.json'
    db_file['nerr'] = destdir + 'conesearch_error.json'
    db_to_use = destdir + 'conesearch.json'

    # JSON dictionaries for output files
    js_template = {'__version__': 1, 'catalogs':{}}
    js_mstr = deepcopy(js_template)
    js_tree = {}
    for key in db_file:
        js_tree[key] = deepcopy(js_template)

        # Delete existing files, if any, to be on the safe side.
        # Else can cause confusion if program exited prior to
        # new files being written but old files are still there.
        _do_rmfile(db_file[key], verbose=verbose)
    _do_rmfile(db_to_use, verbose=verbose)

    # Get all Cone Search sites
    tab_all = votable.parse_single_table(urllib2.urlopen(
        CS_MSTR_LIST(), timeout=vos_catalog.TIMEOUT()), pedantic=False)
    arr_cone = tab_all.array.data[numpy.where(
        tab_all.array['capabilityClass'] == 'ConeSearch')]

    assert arr_cone.size > 0, 'CS_MSTR_LIST yields no valid result'

    col_names = arr_cone.dtype.names
    col_to_rename = {'accessURL':'url'}
    uniq_title = sorted(set(arr_cone['title']))  # May have duplicates
    key_lookup_by_url = {}

    # Re-structure dictionary for JSON file
    for title in uniq_title:
        idx = numpy.where(arr_cone['title'] == title)
        for n,i in enumerate(idx[0], start=1):
            cat_key = '{} {}'.format(title, n)
            row_d = {}
            for col in col_names:
                if col in col_to_rename:
                    row_d[col_to_rename[col]] = arr_cone[i][col]
                else:
                    row_d[col] = arr_cone[i][col]

            cur_url = row_d['url']
            assert cur_url not in key_lookup_by_url, 'Duplicate access URL'
            key_lookup_by_url[cur_url] = cat_key

            js_mstr['catalogs'][cat_key] = row_d

    all_urls = arr_cone['accessURL'] + 'RA=0&DEC=0&SR=0'
    t_beg = time.time()

    # Validate URLs (with multiprocessing)
    if multiproc:
        import multiprocessing
        mp_list = []
        pool = multiprocessing.Pool()
        mp_proc = pool.map_async(_do_validation, all_urls,
                                 callback=mp_list.append)
        mp_proc.wait()
        mp_list = mp_list[0]

    # Validate URLs (no multiprocessing)
    else:
        mp_list = [_do_validation(cur_url) for cur_url in all_urls]

    t_end = time.time()

    if verbose:
        log.info('Validation of {} sites took {} s'.format(
            arr_cone.size, t_end - t_beg))

    # Categorize validation results
    for r in mp_list:
        if r['network_error'] is not None:
            success_key = 'nerr'
        elif r['nexceptions'] > 0:
            success_key = 'excp'
        elif r['nwarnings'] > 0:
            success_key = 'warn'
        else:
            success_key = 'good'

        cat_key = key_lookup_by_url[r.url]
        js_tree[success_key]['catalogs'][cat_key] = js_mstr['catalogs'][cat_key]

    # Write to JSON
    n = {}
    checksum = 0
    for key in db_file:
        n[key] = len(js_tree[key]['catalogs'])
        checksum += n[key]
        if verbose:
            log.info('{}: {} catalogs'.format(key, n[key]))
        with open(db_file[key], 'w') as f_json:
            f_json.write(json.dumps(js_tree[key], cls=NumpyScalarEncoder,
                                    sort_keys=True, indent=4))

    if checksum != arr_cone.size:
        log.warn('Expected {} sites in all databases combined '
                 'but only found {}'.format(arr_cone.size, checksum))

    # Make symbolic link
    if n['good'] > 0:
        _do_symlink(db_file['good'], db_to_use)
    elif n['warn'] > 0:
        _do_symlink(db_file['warn'], db_to_use)
        if verbose:
            log.info('No Cone Search sites cleanly passed validation.')
    else:
        log.warn('All sites have exceptions or errors. '
                 'No viable database for Cone Search.')

    # Delete temp dir
    if os.path.exists(_TMP_ROOT):
        shutil.rmtree(_TMP_ROOT)


def _do_validation(url):
    """Validation for multiprocessing support."""
    r = result.Result(url, root=_TMP_ROOT)
    r.download_xml_content()
    r.validate_vo()
    return r


def _do_rmfile(filename, verbose=True):
    """Delete a file."""
    if os.path.exists(filename):
        assert not os.path.isdir(filename), \
            '{} is a directory, cannot continue'.format(filename)
        os.remove(filename)
        if verbose:
            log.info('Existing copy of {} deleted'.format(filename))


def _do_symlink(src, dst):
    """Create symbolic link (Unix only)."""
    try:
        os.symlink(src, dst)
    except Exception as e:
        log.warn('_do_symlink {} -> {} failed with {}'.format(
            src, dst, e.message))
