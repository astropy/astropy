# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Validate VO Services."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six
from ...extern.six.moves import map

# STDLIB
import multiprocessing
import os
import warnings
from collections import OrderedDict

# LOCAL
from .exceptions import (ValidationMultiprocessingError,
                         InvalidValidationAttribute)
from ..client import vos_catalog
from ..client.exceptions import VOSError
from ...io import votable
from ...io.votable.exceptions import E19
from ...io.votable.validator import html, result
from ...logger import log
from ...utils import data
from ...utils.exceptions import AstropyUserWarning
from ...utils.timer import timefunc
from ...utils.xml.unescaper import unescape_all

# Temporary solution until STScI VAO registry formally provides
# <testQuery> tags
from .tstquery import parse_cs


__all__ = ['check_conesearch_sites']


@timefunc(1)
def check_conesearch_sites(destdir=os.curdir, verbose=True, parallel=True,
                           url_list='default'):
    """Validate Cone Search Services.

    .. note::

        URLs are unescaped prior to validation.

        Only check queries with ``<testQuery>`` parameters.
        Does not perform meta-data and erroneous queries.

    Parameters
    ----------
    destdir : str, optional
        Directory to store output files. Will be created if does
        not exist. Existing files with these names will be deleted
        or replaced:

            * conesearch_good.json
            * conesearch_warn.json
            * conesearch_exception.json
            * conesearch_error.json

    verbose : bool, optional
        Print extra info to log.

    parallel : bool, optional
        Enable multiprocessing.

    url_list : list of string, optional
        Only check these access URLs against
        `astropy.vo.validator.Conf.conesearch_master_list` and ignore
        the others, which will not appear in output files.  By
        default, check those in
        `astropy.vo.validator.Conf.conesearch_urls`.  If `None`, check
        everything.

    Raises
    ------
    IOError
        Invalid destination directory.

    timeout
        URL request timed out.

    ValidationMultiprocessingError
        Multiprocessing failed.

    """
    from . import conf

    if url_list == 'default':
        url_list = conf.conesearch_urls

    if (not isinstance(destdir, six.string_types) or len(destdir) == 0 or
            os.path.exists(destdir) and not os.path.isdir(destdir)):
        raise IOError('Invalid destination directory')  # pragma: no cover

    if not os.path.exists(destdir):
        os.mkdir(destdir)

    # Output dir created by votable.validator
    out_dir = os.path.join(destdir, 'results')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Output files
    db_file = OrderedDict()
    db_file['good'] = os.path.join(destdir, 'conesearch_good.json')
    db_file['warn'] = os.path.join(destdir, 'conesearch_warn.json')
    db_file['excp'] = os.path.join(destdir, 'conesearch_exception.json')
    db_file['nerr'] = os.path.join(destdir, 'conesearch_error.json')

    # JSON dictionaries for output files
    js_tree = {}
    for key in db_file:
        js_tree[key] = vos_catalog.VOSDatabase.create_empty()

        # Delete existing files, if any, to be on the safe side.
        # Else can cause confusion if program exited prior to
        # new files being written but old files are still there.
        if os.path.exists(db_file[key]):  # pragma: no cover
            os.remove(db_file[key])
            if verbose:
                log.info('Existing file {0} deleted'.format(db_file[key]))

    # Master VO database from registry. Silence all the warnings.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        js_mstr = vos_catalog.VOSDatabase.from_registry(
            conf.conesearch_master_list, encoding='binary',
            show_progress=verbose)

    # Validate only a subset of the services.
    if url_list is not None:
        # Make sure URL is unique and fixed.
        url_list = set(map(unescape_all,
            [cur_url.encode('utf-8') if isinstance(cur_url, str) else cur_url
             for cur_url in url_list]))
        uniq_rows = len(url_list)
        url_list_processed = []  # To track if given URL is valid in registry
        if verbose:
            log.info('Only {0}/{1} site(s) are validated'.format(
                uniq_rows, len(js_mstr)))
    # Validate all services.
    else:
        uniq_rows = len(js_mstr)

    key_lookup_by_url = {}

    # Process each catalog in the registry.
    for cur_key, cur_cat in js_mstr.get_catalogs():
        cur_url = cur_cat['url']

        # Skip if:
        #   a. not a Cone Search service
        #   b. not in given subset, if any
        if ((cur_cat['capabilityClass'] != b'ConeSearch') or
                (url_list is not None and cur_url not in url_list)):
            continue

        # Use testQuery to return non-empty VO table with max verbosity.
        testquery_pars = parse_cs(cur_cat['resourceID'])
        cs_pars_arr = ['='.join([key, testquery_pars[key]]).encode('utf-8')
                       for key in testquery_pars]
        cs_pars_arr += [b'VERB=3']

        # Track the service.
        key_lookup_by_url[cur_url + b'&'.join(cs_pars_arr)] = cur_key
        if url_list is not None:
            url_list_processed.append(cur_url)

    # Give warning if any of the user given subset is not in the registry.
    if url_list is not None:
        url_list_skipped = url_list - set(url_list_processed)
        n_skipped = len(url_list_skipped)
        if n_skipped > 0:
            warn_str = '{0} not found in registry! Skipped:\n'.format(n_skipped)
            for cur_url in url_list_skipped:
                warn_str += '\t{0}\n'.format(cur_url)
            warnings.warn(warn_str, AstropyUserWarning)

    all_urls = list(key_lookup_by_url)
    timeout = data.conf.remote_timeout
    map_args = [(out_dir, url, timeout) for url in all_urls]

    # Validate URLs
    if parallel:
        pool = multiprocessing.Pool()
        try:
            mp_list = pool.map(_do_validation, map_args)
        except Exception as exc:  # pragma: no cover
            raise ValidationMultiprocessingError(
                'An exception occurred during parallel processing '
                'of validation results: {0}'.format(exc))
    else:
        mp_list = map(_do_validation, map_args)

    # Categorize validation results
    for r in mp_list:
        db_key = r['out_db_name']
        cat_key = key_lookup_by_url[r.url]
        cur_cat = js_mstr.get_catalog(cat_key)
        _copy_r_to_cat(r, cur_cat)
        js_tree[db_key].add_catalog(cat_key, cur_cat)

    # Write to HTML
    html_subsets = result.get_result_subsets(mp_list, out_dir)
    html.write_index(html_subsets, all_urls, out_dir)
    if parallel:
        html_subindex_args = [(out_dir, html_subset, uniq_rows)
                              for html_subset in html_subsets]
        pool.map(_html_subindex, html_subindex_args)
    else:
        for html_subset in html_subsets:
            _html_subindex((out_dir, html_subset, uniq_rows))

    # Write to JSON
    n = {}
    n_tot = 0
    for key in db_file:
        n[key] = len(js_tree[key])
        n_tot += n[key]
        js_tree[key].to_json(db_file[key], overwrite=True)
        if verbose:
            log.info('{0}: {1} catalog(s)'.format(key, n[key]))

    # Checksum
    if verbose:
        log.info('total: {0} out of {1} catalog(s)'.format(n_tot, uniq_rows))

    if n['good'] == 0:  # pragma: no cover
        warnings.warn(
            'No good sites available for Cone Search.', AstropyUserWarning)


def _do_validation(args):
    """Validation for multiprocessing support."""

    root, url, timeout = args

    votable.table.reset_vo_warnings()

    r = result.Result(url, root=root, timeout=timeout)
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
            except (E19, IndexError, VOSError) as e:  # pragma: no cover
                lines.append(str(e))
                nexceptions += 1
        lines = [str(x.message) for x in warning_lines] + lines

        warning_types = set()
        for line in lines:  # pragma: no cover
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
    r : `astropy.io.votable.validator.result.Result`

    Raises
    ------
    InvalidValidationAttribute
        Unhandled validation result attributes.

    """
    from . import conf

    if 'network_error' in r and r['network_error'] is not None:  # pragma: no cover
        r['out_db_name'] = 'nerr'
        r['expected'] = 'broken'
    elif ((r['nexceptions'] == 0 and r['nwarnings'] == 0) or
            r['warning_types'].issubset(conf.noncritical_warnings)):
        r['out_db_name'] = 'good'
        r['expected'] = 'good'
    elif r['nexceptions'] > 0:  # pragma: no cover
        r['out_db_name'] = 'excp'
        r['expected'] = 'incorrect'
    elif r['nwarnings'] > 0:  # pragma: no cover
        r['out_db_name'] = 'warn'
        r['expected'] = 'incorrect'
    else:  # pragma: no cover
        raise InvalidValidationAttribute(
            'Unhandled validation result attributes: {0}'.format(r._attributes))


def _html_subindex(args):
    """HTML writer for multiprocessing support."""
    out_dir, subset, total = args
    html.write_index_table(out_dir, *subset, total=total)


def _copy_r_to_cat(r, cat):
    """Copy validation result attributes to given VO catalog.

    Parameters
    ----------
    r : `astropy.io.votable.validate.result.Result`

    cat : `astropy.vo.client.vos_catalog.VOSCatalog`

    """
    for key in r._attributes:
        new_key = 'validate_' + key
        cat[new_key] = r._attributes[key]
