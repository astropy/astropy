# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Validate VO Services."""
from __future__ import print_function, division

# STDLIB
from collections import defaultdict
from copy import deepcopy
import json
import os
import time
import warnings

# THIRD PARTY
import numpy as np

# LOCAL
from ..client import vos_catalog
from ...config.configuration import ConfigurationItem
from ...io import votable
from ...io.votable.exceptions import E19
from ...io.votable.validator import html, result
from ...logger import log
from ...utils.data import get_readable_fileobj, get_pkg_data_contents
from ...utils.data import REMOTE_TIMEOUT
from ...utils.misc import JsonCustomEncoder
from ...utils.xml.unescaper import unescape_all

# Temporary solution until STScI VAO registry formally provides
# <testQuery> tags
from .tstquery import parse_cs


__all__ = ['check_conesearch_sites']

CS_MSTR_LIST = ConfigurationItem(
    'cs_mstr_list',
    'http://vao.stsci.edu/directory/NVORegInt.asmx/VOTCapabilityPredOpt?'
    'predicate=1%3D1&capability=conesearch&VOTStyleOption=2',
    'Cone Search services master list for validation.')

CS_URLS = ConfigurationItem(
    'cs_urls',
    get_pkg_data_contents(
        os.path.join('data', 'conesearch_urls.txt')).split(),
    'Only check these Cone Search URLs.',
    'list')

NONCRIT_WARNINGS = ConfigurationItem(
    'noncrit_warnings',
    ['W03', 'W06', 'W07', 'W09', 'W10', 'W15', 'W17', 'W20', 'W21', 'W22',
     'W27', 'W28', 'W29', 'W41', 'W42', 'W48', 'W50'],
    'VO Table warning codes that are considered non-critical',
    'list')

_OUT_ROOT = None  # Set by `check_conesearch_sites`


def check_conesearch_sites(destdir=os.curdir, verbose=True, parallel=True,
                           url_list=CS_URLS()):
    """Validate Cone Search Services.

    .. note::

        URLs are unescaped prior to validation.

        Only check queries with real matched objects.
        Does not perform meta-data and erroneous queries.

    Parameters
    ----------
    destdir : string
        Directory to store output files. Will be created if does
        not exist. Existing files with these names will be deleted
        or replaced:
            * conesearch_good.json
            * conesearch_warn.json
            * conesearch_exception.json
            * conesearch_error.json

    verbose : bool
        Print extra info to log.

    parallel : bool
        Enable multiprocessing.

    url_list : list of string
        Only check these access URLs against
        `astropy.vo.validator.validate.CS_MSTR_LIST` and ignore the others,
        which will not appear in output files.
        By default, check those in `astropy.vo.validator.validate.CS_URLS`.
        If `None`, check everything.

    Raises
    ------
    AssertionError
        Parameter failed assertion test.

    IOError
        Invalid destination directory.

    timeout
        URL request timed out.

    """
    global _OUT_ROOT

    # Start timer
    t_beg = time.time()

    if (not isinstance(destdir, basestring) or len(destdir) == 0 or
            os.path.exists(destdir) and not os.path.isdir(destdir)):
        raise IOError('Invalid destination directory')  # pragma: no cover

    if not os.path.exists(destdir):
        os.mkdir(destdir)

    # Output dir created by votable.validator
    _OUT_ROOT = os.path.join(destdir, 'results')

    if not os.path.exists(_OUT_ROOT):
        os.mkdir(_OUT_ROOT)

    # Output files
    db_file = {}
    db_file['good'] = os.path.join(destdir, 'conesearch_good.json')
    db_file['warn'] = os.path.join(destdir, 'conesearch_warn.json')
    db_file['excp'] = os.path.join(destdir, 'conesearch_exception.json')
    db_file['nerr'] = os.path.join(destdir, 'conesearch_error.json')

    # JSON dictionaries for output files
    js_template = {'__version__': 1, 'catalogs': {}}
    js_mstr = deepcopy(js_template)
    js_tree = {}
    for key in db_file:
        js_tree[key] = deepcopy(js_template)

        # Delete existing files, if any, to be on the safe side.
        # Else can cause confusion if program exited prior to
        # new files being written but old files are still there.
        if os.path.exists(db_file[key]):  # pragma: no cover
            os.remove(db_file[key])
            if verbose:
                log.info('Existing file {0} deleted'.format(db_file[key]))

    # Get all Cone Search sites
    with get_readable_fileobj(CS_MSTR_LIST(), encoding='binary') as fd:
        tab_all = votable.parse_single_table(fd, pedantic=False)
    arr_cone = tab_all.array.data[np.where(
        tab_all.array['capabilityClass'] == b'ConeSearch')]

    assert arr_cone.size > 0, \
        'astropy.vo.validator.validate.CS_MSTR_LIST yields no valid result'

    fixed_urls = [unescape_all(cur_url) for cur_url in arr_cone['accessURL']]
    uniq_urls = set(fixed_urls)

    if url_list is None:
        url_list = uniq_urls
    else:
        tmp_list = [cur_url.encode('utf-8') if isinstance(cur_url, str)
                    else cur_url for cur_url in set(url_list)]
        url_list = [unescape_all(cur_url) for cur_url in tmp_list]

        if verbose:
            log.info('Only {0}/{1} site(s) are validated'.format(
                len(url_list), len(uniq_urls)))

    uniq_rows = len(url_list)

    # Re-structure dictionary for JSON file

    col_names = tab_all.array.dtype.names
    title_counter = defaultdict(int)
    key_lookup_by_url = {}

    for cur_url in url_list:
        num_match = fixed_urls.count(cur_url)

        if num_match == 0:
            log.warn(
                '{0} not found in cs_mstr_list! Skipping...'.format(cur_url))
            continue

        i = fixed_urls.index(cur_url)
        n_ignored = num_match - 1
        row_d = {'duplicatesIgnored': n_ignored}
        if verbose and n_ignored > 0:  # pragma: no cover
            log.info('{0} has {1} ignored duplicate entries in '
                     'cs_mstr_list'.format(cur_url, n_ignored))

        cur_title = arr_cone[i]['title']
        title_counter[cur_title] += 1

        if isinstance(cur_title, bytes):  # pragma: py3
            cat_key = '{0} {1}'.format(cur_title.decode('ascii'),
                                       title_counter[cur_title])
        else:  # pragma: py2
            cat_key = '{0} {1}'.format(cur_title, title_counter[cur_title])

        for col in col_names:
            if col == 'accessURL':
                row_d['url'] = fixed_urls[i]
            else:
                row_d[col] = arr_cone[i][col]

        # Use testQuery to return non-empty VO table
        testquery_pars = parse_cs(arr_cone[i]['resourceID'])
        cs_pars_arr = ['='.join([key, testquery_pars[key]]).encode('utf-8')
                       for key in testquery_pars]

        # Max verbosity
        cs_pars_arr += [b'VERB=3']
        js_mstr['catalogs'][cat_key] = row_d
        key_lookup_by_url[cur_url + b'&'.join(cs_pars_arr)] = cat_key

    # Validate URLs

    all_urls = key_lookup_by_url.keys()

    if parallel:
        from multiprocessing import Pool
        mp_list = []
        pool = Pool()
        mp_proc = pool.map_async(_do_validation, all_urls,
                                 callback=mp_list.append)
        mp_proc.wait()
        assert len(mp_list) > 0, \
            'Multiprocessing pool callback returned empty list'
        mp_list = mp_list[0]

    else:
        mp_list = [_do_validation(cur_url) for cur_url in all_urls]

    # Categorize validation results
    for r in mp_list:
        db_key = r['out_db_name']
        cat_key = key_lookup_by_url[r.url]
        js_tree[db_key]['catalogs'][cat_key] = js_mstr['catalogs'][cat_key]
        _copy_r_to_db(r, js_tree[db_key]['catalogs'][cat_key])

    # Write to HTML

    html_subsets = result.get_result_subsets(mp_list, _OUT_ROOT)
    html.write_index(html_subsets, all_urls, _OUT_ROOT)

    if parallel:
        html_subindex_args = [(html_subset, uniq_rows)
                              for html_subset in html_subsets]
        mp_proc = pool.map_async(_html_subindex, html_subindex_args)
        mp_proc.wait()

    else:
        for html_subset in html_subsets:
            _html_subindex((html_subset, uniq_rows))

    # Write to JSON
    n = {}
    n_tot = 0
    for key in db_file:
        n[key] = len(js_tree[key]['catalogs'])
        n_tot += n[key]
        if verbose:
            log.info('{0}: {1} catalog(s)'.format(key, n[key]))
        with open(db_file[key], 'w') as f_json:
            f_json.write(json.dumps(js_tree[key], cls=JsonCustomEncoder,
                                    sort_keys=True, indent=4))

    # End timer
    t_end = time.time()

    if verbose:
        log.info('total: {0} catalog(s)'.format(n_tot))
        log.info('Validation of {0} site(s) took {1:.3f} s'.format(
            uniq_rows, t_end - t_beg))

    if n['good'] == 0:  # pragma: no cover
        log.warn('No good sites available for Cone Search.')


def _do_validation(url):
    """Validation for multiprocessing support."""
    votable.table.reset_vo_warnings()

    r = result.Result(url, root=_OUT_ROOT, timeout=REMOTE_TIMEOUT())
    r.validate_vo()

    _categorize_result(r)

    # This was already checked above.
    # Calling this again to get VOTableFile object to catch
    # well-formed error responses in downloaded XML.
    #
    # 'incorrect' is also added in case user wants to use
    # 'conesearch_warn.json' anyway.
    #
    # If using cached data, it will not detect network error
    # like the first run, but will raise exception.
    #
    # When SR is not 0, VOSError is raised for empty table.
    #
    if r['expected'] in ('good', 'incorrect') and r['nexceptions'] == 0:
        nexceptions = 0
        nwarnings = 0
        lines = []

        with warnings.catch_warnings(record=True) as warning_lines:
            try:
                tab = vos_catalog.vo_tab_parse(votable.table.parse(
                    r.get_vo_xml_path(), pedantic=False), r.url, {})
            except (E19, IndexError, vos_catalog.VOSError) as e:
                lines.append(str(e))
                nexceptions += 1
        lines = [str(x.message) for x in warning_lines] + lines

        warning_types = set()
        for line in lines:
            w = votable.exceptions.parse_vowarning(line)
            if w['is_warning']:
                nwarnings += 1
            if w['is_exception']:
                nexceptions += 1
            warning_types.add(w['warning'])

        r['nwarnings'] += nwarnings
        r['nexceptions'] += nexceptions
        r['warnings'] += lines
        r['warning_types'] = r['warning_types'].union(warning_types)

        _categorize_result(r)

    html.write_result(r)
    return r


def _categorize_result(r):
    """Set success codes.

    Parameters
    ----------
    r : `astropy.io.votable.validator.result.Result` object

    """
    if 'network_error' in r and r['network_error'] is not None:
        r['out_db_name'] = 'nerr'
        r['expected'] = 'broken'
    elif ((r['nexceptions'] == 0 and r['nwarnings'] == 0) or
            r['warning_types'].issubset(NONCRIT_WARNINGS())):
        r['out_db_name'] = 'good'
        r['expected'] = 'good'
    elif r['nexceptions'] > 0:
        r['out_db_name'] = 'excp'
        r['expected'] = 'incorrect'
    elif r['nwarnings'] > 0:
        r['out_db_name'] = 'warn'
        r['expected'] = 'incorrect'
    else:  # pragma: no cover
        raise vos_catalog.VOSError(
            'Unhandled validation result attributes: '
            '{0}'.format(r._attributes))


def _html_subindex(args):
    """HTML writer for multiprocessing support."""
    subset, total = args
    html.write_index_table(_OUT_ROOT, *subset, total=total)


def _copy_r_to_db(r, db):
    """Copy validation result attributes to given JSON database entry.

    Parameters
    ----------
    r : `astropy.io.votable.validate.result.Result` object

    db : dict

    """
    for key in r._attributes:
        new_key = 'validate_' + key
        db[new_key] = r._attributes[key]
