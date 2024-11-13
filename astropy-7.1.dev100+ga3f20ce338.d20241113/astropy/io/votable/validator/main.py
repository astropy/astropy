# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Validates a large collection of web-accessible VOTable files,
and generates a report as a directory tree of HTML files.
"""

# STDLIB
import os

# LOCAL
from astropy.utils.data import get_pkg_data_filename

from . import html, result

__all__ = ["make_validation_report"]


def get_srcdir():
    return os.path.dirname(__file__)


def get_urls(destdir, s):
    import gzip

    types = ["good", "broken", "incorrect"]

    seen = set()
    urls = []
    for type in types:
        filename = get_pkg_data_filename(f"data/urls/cone.{type}.dat.gz")
        with gzip.open(filename, "rb") as fd:
            for url in fd.readlines():
                next(s)
                url = url.strip()
                if url not in seen:
                    with result.Result(url, root=destdir) as r:
                        r["expected"] = type
                    urls.append(url)
                seen.add(url)

    return urls


def download(args):
    url, destdir = args
    with result.Result(url, root=destdir) as r:
        r.download_xml_content()


def validate_vo(args):
    url, destdir = args
    with result.Result(url, root=destdir) as r:
        r.validate_vo()


def votlint_validate(args):
    path_to_stilts_jar, url, destdir = args
    with result.Result(url, root=destdir) as r:
        if r["network_error"] is None:
            r.validate_with_votlint(path_to_stilts_jar)


def write_html_result(args):
    url, destdir = args
    with result.Result(url, root=destdir) as r:
        html.write_result(r)


def write_subindex(args):
    subset, destdir, total = args
    html.write_index_table(destdir, *subset, total=total)


def make_validation_report(
    urls=None,
    destdir="astropy.io.votable.validator.results",
    multiprocess=True,
    stilts=None,
):
    """
    Validates a large collection of web-accessible VOTable files.

    Generates a report as a directory tree of HTML files.

    Parameters
    ----------
    urls : list of str, optional
        If provided, is a list of HTTP urls to download VOTable files
        from.  If not provided, a built-in set of ~22,000 urls
        compiled by HEASARC will be used.

    destdir : path-like, optional
        The directory to write the report to.  By default, this is a
        directory called ``'results'`` in the current directory. If the
        directory does not exist, it will be created.

    multiprocess : bool, optional
        If `True` (default), perform validations in parallel using all
        of the cores on this machine.

    stilts : path-like, optional
        To perform validation with ``votlint`` from the Java-based |STILTS|
        VOTable parser, in addition to `astropy.io.votable`, set this to the
        path of the ``'stilts.jar'`` file.  ``java`` on the system shell
        path will be used to run it.

    Notes
    -----
    Downloads of each given URL will be performed only once and cached
    locally in *destdir*.  To refresh the cache, remove *destdir*
    first.
    """
    from astropy.utils.console import ProgressBar, Spinner, color_print

    if stilts is not None:
        if not os.path.exists(stilts):
            raise ValueError(f"{stilts} does not exist.")

    destdir = os.path.expanduser(destdir)
    destdir = os.path.abspath(destdir)

    if urls is None:
        with Spinner("Loading URLs", "green") as s:
            urls = get_urls(destdir, s)
    else:
        urls = [url.encode() for url in urls if isinstance(url, str)]

        color_print("Marking URLs", "green")
        for url in ProgressBar(urls):
            with result.Result(url, root=destdir) as r:
                r["expected"] = type

    args = [(url, destdir) for url in urls]

    color_print("Downloading VO files", "green")
    ProgressBar.map(download, args, multiprocess=multiprocess)

    color_print("Validating VO files", "green")
    ProgressBar.map(validate_vo, args, multiprocess=multiprocess)

    if stilts is not None:
        color_print("Validating with votlint", "green")
        votlint_args = [(stilts, x, destdir) for x in urls]
        ProgressBar.map(votlint_validate, votlint_args, multiprocess=multiprocess)

    color_print("Generating HTML files", "green")
    ProgressBar.map(write_html_result, args, multiprocess=multiprocess)

    with Spinner("Grouping results", "green") as s:
        subsets = result.get_result_subsets(urls, destdir, s)

    color_print("Generating index", "green")
    html.write_index(subsets, urls, destdir)

    color_print("Generating subindices", "green")
    subindex_args = [(subset, destdir, len(urls)) for subset in subsets]
    ProgressBar.map(write_subindex, subindex_args, multiprocess=multiprocess)
