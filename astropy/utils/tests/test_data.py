# Licensed under a 3-clause BSD style license - see LICENSE.rst

import base64
import contextlib
import errno
import hashlib
import io
import itertools
import os
import pathlib
import platform
import random
import shutil
import stat
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import warnings
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from itertools import islice
from tempfile import NamedTemporaryFile, TemporaryDirectory

import pytest

import astropy.utils.data
from astropy import units as _u  # u is taken
from astropy.config import paths
from astropy.tests.helper import CI, IS_CRON, PYTEST_LT_8_0
from astropy.utils.compat.optional_deps import HAS_BZ2, HAS_LZMA
from astropy.utils.data import (
    CacheDamaged,
    CacheMissingWarning,
    _deltemps,
    _get_download_cache_loc,
    _tempfilestodel,
    cache_contents,
    cache_total_size,
    check_download_cache,
    check_free_space_in_dir,
    clear_download_cache,
    compute_hash,
    conf,
    download_file,
    download_files_in_parallel,
    export_download_cache,
    get_cached_urls,
    get_file_contents,
    get_free_space_in_dir,
    get_pkg_data_contents,
    get_pkg_data_filename,
    get_pkg_data_fileobj,
    get_pkg_data_path,
    get_readable_fileobj,
    import_download_cache,
    import_file_to_cache,
    is_url,
    is_url_in_cache,
)
from astropy.utils.exceptions import AstropyWarning

TESTURL = "http://www.astropy.org"
TESTURL2 = "http://www.astropy.org/about.html"
TESTURL_SSL = "https://www.astropy.org"
TESTLOCAL = get_pkg_data_filename(os.path.join("data", "local.dat"))

# For when we need "some" test URLs
FEW = 5

# For stress testing the locking system using multiprocessing
N_PARALLEL_HAMMER = 5  # as high as 500 to replicate a bug

# For stress testing the locking system using threads
# (cheaper, works with coverage)
N_THREAD_HAMMER = 10  # as high as 1000 to replicate a bug


def can_rename_directory_in_use():
    with TemporaryDirectory() as d:
        d1 = os.path.join(d, "a")
        d2 = os.path.join(d, "b")
        f1 = os.path.join(d1, "file")
        os.mkdir(d1)
        with open(f1, "w") as f:
            f.write("some contents\n")
        try:
            with open(f1):
                os.rename(d1, d2)
        except PermissionError:
            return False
        else:
            return True


CAN_RENAME_DIRECTORY_IN_USE = can_rename_directory_in_use()


def url_to(path):
    return pathlib.Path(path).resolve().as_uri()


@pytest.fixture
def valid_urls(tmp_path):
    def _valid_urls(tmp_path):
        for i in itertools.count():
            c = os.urandom(16).hex()
            fn = tmp_path / f"valid_{i}"
            with open(fn, "w") as f:
                f.write(c)
            u = url_to(fn)
            yield u, c

    return _valid_urls(tmp_path)


@pytest.fixture
def invalid_urls(tmp_path):
    def _invalid_urls(tmp_path):
        for i in itertools.count():
            fn = tmp_path / f"invalid_{i}"
            if not os.path.exists(fn):
                yield url_to(fn)

    return _invalid_urls(tmp_path)


@pytest.fixture
def temp_cache(tmp_path):
    with paths.set_temp_cache(tmp_path):
        yield None
        check_download_cache()


def change_tree_permission(d, writable=False):
    if writable:
        dirperm = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
        fileperm = stat.S_IRUSR | stat.S_IWUSR
    else:
        dirperm = stat.S_IRUSR | stat.S_IXUSR
        fileperm = stat.S_IRUSR
    for dirpath, dirnames, filenames in os.walk(d):
        os.chmod(dirpath, dirperm)
        for f in filenames:
            os.chmod(os.path.join(dirpath, f), fileperm)


def is_dir_readonly(d):
    try:
        with NamedTemporaryFile(dir=d):
            return False
    except PermissionError:
        return True


@contextlib.contextmanager
def readonly_dir(d):
    try:
        change_tree_permission(d, writable=False)
        yield
    finally:
        change_tree_permission(d, writable=True)


@pytest.fixture
def readonly_cache(tmp_path, valid_urls):
    with TemporaryDirectory(dir=tmp_path) as d:
        # other fixtures use the same tmp_path so we need a subdirectory
        # to make into the cache
        d = pathlib.Path(d)
        with paths.set_temp_cache(d):
            us = {u for u, c in islice(valid_urls, FEW)}
            urls = {u: download_file(u, cache=True) for u in us}
            files = set(d.iterdir())
            with readonly_dir(d):
                if not is_dir_readonly(d):
                    pytest.skip("Unable to make directory readonly")
                yield urls
            assert set(d.iterdir()) == files
            check_download_cache()


@pytest.fixture
def fake_readonly_cache(tmp_path, valid_urls, monkeypatch):
    def no_mkdir(path, mode=None):
        raise OSError(errno.EPERM, "os.mkdir monkeypatched out")

    def no_mkdtemp(*args, **kwargs):
        """On Windows, mkdtemp uses mkdir in a loop and therefore hangs
        with it monkeypatched out.
        """
        raise OSError(errno.EPERM, "os.mkdtemp monkeypatched out")

    def no_TemporaryDirectory(*args, **kwargs):
        raise OSError(errno.EPERM, "_SafeTemporaryDirectory monkeypatched out")

    with TemporaryDirectory(dir=tmp_path) as d:
        # other fixtures use the same tmp_path so we need a subdirectory
        # to make into the cache
        d = pathlib.Path(d)
        with paths.set_temp_cache(d):
            us = {u for u, c in islice(valid_urls, FEW)}
            urls = {u: download_file(u, cache=True) for u in us}
            files = set(d.iterdir())
            monkeypatch.setattr(os, "mkdir", no_mkdir)
            monkeypatch.setattr(tempfile, "mkdtemp", no_mkdtemp)
            monkeypatch.setattr(
                astropy.utils.data, "_SafeTemporaryDirectory", no_TemporaryDirectory
            )
            yield urls
            assert set(d.iterdir()) == files
            check_download_cache()


def test_download_file_basic(valid_urls, temp_cache):
    u, c = next(valid_urls)
    assert get_file_contents(download_file(u, cache=False)) == c
    assert not is_url_in_cache(u)
    assert get_file_contents(download_file(u, cache=True)) == c  # Cache miss
    assert is_url_in_cache(u)
    assert get_file_contents(download_file(u, cache=True)) == c  # Cache hit
    assert get_file_contents(download_file(u, cache=True, sources=[])) == c


def test_download_file_absolute_path(valid_urls, temp_cache):
    def is_abs(p):
        return p == os.path.abspath(p)

    u, c = next(valid_urls)
    assert is_abs(download_file(u, cache=False))  # no cache
    assert is_abs(download_file(u, cache=True))  # not in cache
    assert is_abs(download_file(u, cache=True))  # in cache
    for v in cache_contents().values():
        assert is_abs(v)


def test_unicode_url(valid_urls, temp_cache):
    u, c = next(valid_urls)
    unicode_url = "http://é—☃—è.com"
    download_file(unicode_url, cache=False, sources=[u])
    download_file(unicode_url, cache=True, sources=[u])
    download_file(unicode_url, cache=True, sources=[])
    assert is_url_in_cache(unicode_url)
    assert unicode_url in cache_contents()


def test_too_long_url(valid_urls, temp_cache):
    u, c = next(valid_urls)
    long_url = "http://" + "a" * 256 + ".com"
    download_file(long_url, cache=False, sources=[u])
    download_file(long_url, cache=True, sources=[u])
    download_file(long_url, cache=True, sources=[])


def test_case_collision(valid_urls, temp_cache):
    u, c = next(valid_urls)
    u2, c2 = next(valid_urls)
    f1 = download_file("http://example.com/thing", cache=True, sources=[u])
    f2 = download_file("http://example.com/THING", cache=True, sources=[u2])
    assert f1 != f2
    assert get_file_contents(f1) != get_file_contents(f2)


def test_domain_name_case(valid_urls, temp_cache):
    u, c = next(valid_urls)
    download_file("http://Example.com/thing", cache=True, sources=[u])
    assert is_url_in_cache("http://EXAMPLE.com/thing")
    download_file("http://EXAMPLE.com/thing", cache=True, sources=[])
    assert is_url_in_cache("Http://example.com/thing")
    download_file("Http://example.com/thing", cache=True, sources=[])


@pytest.mark.remote_data(source="astropy")
def test_download_nocache_from_internet():
    fnout = download_file(TESTURL, cache=False)
    assert os.path.isfile(fnout)


@pytest.fixture
def a_binary_file(tmp_path):
    fn = tmp_path / "file"
    b_contents = b"\xde\xad\xbe\xef"
    with open(fn, "wb") as f:
        f.write(b_contents)
    yield fn, b_contents


@pytest.fixture
def a_file(tmp_path):
    fn = tmp_path / "file.txt"
    contents = "contents\n"
    with open(fn, "w") as f:
        f.write(contents)
    yield fn, contents


def test_temp_cache(tmp_path):
    dldir0 = _get_download_cache_loc()
    check_download_cache()

    with paths.set_temp_cache(tmp_path):
        dldir1 = _get_download_cache_loc()
        check_download_cache()
        assert dldir1 != dldir0

    dldir2 = _get_download_cache_loc()
    check_download_cache()
    assert dldir2 != dldir1
    assert dldir2 == dldir0

    # Check that things are okay even if we exit via an exception
    class Special(Exception):
        pass

    try:
        with paths.set_temp_cache(tmp_path):
            dldir3 = _get_download_cache_loc()
            check_download_cache()
            assert dldir3 == dldir1
            raise Special
    except Special:
        pass

    dldir4 = _get_download_cache_loc()
    check_download_cache()
    assert dldir4 != dldir3
    assert dldir4 == dldir0


@pytest.mark.parametrize("parallel", [False, True])
def test_download_with_sources_and_bogus_original(
    valid_urls, invalid_urls, temp_cache, parallel
):
    # This is a combined test because the parallel version triggered a nasty
    # bug and I was trying to track it down by comparing with the non-parallel
    # version. I think the bug was that the parallel downloader didn't respect
    # temporary cache settings.

    # Make a big list of test URLs
    u, c = next(valid_urls)
    # as tuples (URL, right_content, wrong_content)
    urls = [(u, c, None)]
    # where to download the contents
    sources = {}
    # Set up some URLs to download where the "true" URL is not in the sources
    # list; make the true URL valid with different contents so we can tell if
    # it was loaded by mistake.
    for i, (um, c_bad) in enumerate(islice(valid_urls, FEW)):
        assert not is_url_in_cache(um)
        # For many of them the sources list starts with invalid URLs
        sources[um] = list(islice(invalid_urls, i))
        u, c = next(valid_urls)
        sources[um].append(u)
        urls.append((um, c, c_bad))

    # Now fetch them all
    if parallel:
        rs = download_files_in_parallel(
            [u for (u, c, c_bad) in urls], cache=True, sources=sources
        )
    else:
        rs = [
            download_file(u, cache=True, sources=sources.get(u))
            for (u, c, c_bad) in urls
        ]
    assert len(rs) == len(urls)
    for r, (u, c, c_bad) in zip(rs, urls):
        assert get_file_contents(r) == c
        assert get_file_contents(r) != c_bad
        assert is_url_in_cache(u)


@pytest.mark.skipif(
    (sys.platform.startswith("win") and CI), reason="flaky cache error on Windows CI"
)
def test_download_file_threaded_many(temp_cache, valid_urls):
    """Hammer download_file with multiple threaded requests.

    The goal is to stress-test the locking system. Normal parallel downloading
    also does this but coverage tools lose track of which paths are explored.

    """
    urls = list(islice(valid_urls, N_THREAD_HAMMER))
    with ThreadPoolExecutor(max_workers=len(urls)) as P:
        r = list(P.map(lambda u: download_file(u, cache=True), [u for (u, c) in urls]))
    check_download_cache()
    assert len(r) == len(urls)
    for r_, (u, c) in zip(r, urls):
        assert get_file_contents(r_) == c


@pytest.mark.skipif(
    (sys.platform.startswith("win") and CI), reason="flaky cache error on Windows CI"
)
def test_threaded_segfault(valid_urls):
    """Demonstrate urllib's segfault."""

    def slurp_url(u):
        with urllib.request.urlopen(u) as remote:
            block = True
            while block:
                block = remote.read(1024)

    urls = list(islice(valid_urls, N_THREAD_HAMMER))
    with ThreadPoolExecutor(max_workers=len(urls)) as P:
        list(P.map(lambda u: slurp_url(u), [u for (u, c) in urls]))


@pytest.mark.skipif(
    (sys.platform.startswith("win") and CI), reason="flaky cache error on Windows CI"
)
def test_download_file_threaded_many_partial_success(
    temp_cache, valid_urls, invalid_urls
):
    """Hammer download_file with multiple threaded requests.

    Because some of these requests fail, the locking context manager is
    exercised with exceptions as well as success returns. I do not expect many
    surprises from the threaded version, but the process version gave trouble
    here.

    """
    urls = []
    contents = {}
    for (u, c), i in islice(zip(valid_urls, invalid_urls), N_THREAD_HAMMER):
        urls.append(u)
        contents[u] = c
        urls.append(i)

    def get(u):
        try:
            return download_file(u, cache=True)
        except OSError:
            return None

    with ThreadPoolExecutor(max_workers=len(urls)) as P:
        r = list(P.map(get, urls))
    check_download_cache()
    assert len(r) == len(urls)
    for r_, u in zip(r, urls):
        if u in contents:
            assert get_file_contents(r_) == contents[u]
        else:
            assert r_ is None


def test_clear_download_cache(valid_urls):
    u1, c1 = next(valid_urls)
    download_file(u1, cache=True)

    u2, c2 = next(valid_urls)
    download_file(u2, cache=True)

    assert is_url_in_cache(u2)
    clear_download_cache(u2)
    assert not is_url_in_cache(u2)
    assert is_url_in_cache(u1)

    u3, c3 = next(valid_urls)
    f3 = download_file(u3, cache=True)

    assert is_url_in_cache(u3)
    clear_download_cache(f3)
    assert not is_url_in_cache(u3)
    assert is_url_in_cache(u1)

    u4, c4 = next(valid_urls)
    f4 = download_file(u4, cache=True)

    assert is_url_in_cache(u4)
    clear_download_cache(compute_hash(f4))
    assert not is_url_in_cache(u4)
    assert is_url_in_cache(u1)


def test_clear_download_multiple_references_doesnt_corrupt_storage(
    temp_cache, tmp_path
):
    """Check that files with the same hash don't confuse the storage."""
    content = "Test data; doesn't matter much.\n"

    def make_url():
        with NamedTemporaryFile("w", dir=tmp_path, delete=False) as f:
            f.write(content)
        url = url_to(f.name)
        clear_download_cache(url)
        filename = download_file(url, cache=True)
        return url, filename

    a_url, a_filename = make_url()
    clear_download_cache(a_filename)
    assert not is_url_in_cache(a_url)

    f_url, f_filename = make_url()
    g_url, g_filename = make_url()

    assert f_url != g_url
    assert is_url_in_cache(f_url)
    assert is_url_in_cache(g_url)

    clear_download_cache(f_url)
    assert not is_url_in_cache(f_url)
    assert is_url_in_cache(g_url)
    assert os.path.exists(
        g_filename
    ), "Contents should not be deleted while a reference exists"

    clear_download_cache(g_url)
    assert not os.path.exists(
        g_filename
    ), "No reference exists any more, file should be deleted"


@pytest.mark.parametrize("use_cache", [False, True])
def test_download_file_local_cache_survives(tmp_path, temp_cache, use_cache):
    """Confirm that downloading a local file does not delete it.

    When implemented with urlretrieve (rather than urlopen) local files are
    not copied to create temporaries, so importing them to the cache deleted
    the original from wherever it was in the filesystem. I lost some built-in
    astropy data.

    """
    fn = tmp_path / "file"
    contents = "some text"
    with open(fn, "w") as f:
        f.write(contents)
    u = url_to(fn)
    f = download_file(u, cache=use_cache)
    assert fn not in _tempfilestodel, "File should not be deleted!"
    assert os.path.isfile(fn), "File should not be deleted!"
    assert get_file_contents(f) == contents


def test_sources_normal(temp_cache, valid_urls, invalid_urls):
    primary, contents = next(valid_urls)
    fallback1 = next(invalid_urls)
    f = download_file(primary, cache=True, sources=[primary, fallback1])
    assert get_file_contents(f) == contents
    assert is_url_in_cache(primary)
    assert not is_url_in_cache(fallback1)


def test_sources_fallback(temp_cache, valid_urls, invalid_urls):
    primary = next(invalid_urls)
    fallback1, contents = next(valid_urls)
    f = download_file(primary, cache=True, sources=[primary, fallback1])
    assert get_file_contents(f) == contents
    assert is_url_in_cache(primary)
    assert not is_url_in_cache(fallback1)


def test_sources_ignore_primary(temp_cache, valid_urls, invalid_urls):
    primary, bogus = next(valid_urls)
    fallback1, contents = next(valid_urls)
    f = download_file(primary, cache=True, sources=[fallback1])
    assert get_file_contents(f) == contents
    assert is_url_in_cache(primary)
    assert not is_url_in_cache(fallback1)


def test_sources_multiple(temp_cache, valid_urls, invalid_urls):
    primary = next(invalid_urls)
    fallback1 = next(invalid_urls)
    fallback2, contents = next(valid_urls)
    f = download_file(primary, cache=True, sources=[primary, fallback1, fallback2])
    assert get_file_contents(f) == contents
    assert is_url_in_cache(primary)
    assert not is_url_in_cache(fallback1)
    assert not is_url_in_cache(fallback2)


def test_sources_multiple_missing(temp_cache, valid_urls, invalid_urls):
    primary = next(invalid_urls)
    fallback1 = next(invalid_urls)
    fallback2 = next(invalid_urls)
    with pytest.raises(urllib.error.URLError):
        download_file(primary, cache=True, sources=[primary, fallback1, fallback2])
    assert not is_url_in_cache(primary)
    assert not is_url_in_cache(fallback1)
    assert not is_url_in_cache(fallback2)


def test_update_url(tmp_path, temp_cache):
    with TemporaryDirectory(dir=tmp_path) as d:
        f_name = os.path.join(d, "f")
        with open(f_name, "w") as f:
            f.write("old")
        f_url = url_to(f.name)
        assert get_file_contents(download_file(f_url, cache=True)) == "old"
        with open(f_name, "w") as f:
            f.write("new")
        assert get_file_contents(download_file(f_url, cache=True)) == "old"
        assert get_file_contents(download_file(f_url, cache="update")) == "new"
    # Now the URL doesn't exist any more.
    assert not os.path.exists(f_name)
    with pytest.raises(urllib.error.URLError):
        # Direct download should fail
        download_file(f_url, cache=False)
    assert (
        get_file_contents(download_file(f_url, cache=True)) == "new"
    ), "Cached version should still exist"
    with pytest.raises(urllib.error.URLError):
        # cannot download new version to check for updates
        download_file(f_url, cache="update")
    assert (
        get_file_contents(download_file(f_url, cache=True)) == "new"
    ), "Failed update should not remove the current version"


@pytest.mark.remote_data(source="astropy")
def test_download_noprogress():
    fnout = download_file(TESTURL, cache=False, show_progress=False)
    assert os.path.isfile(fnout)


@pytest.mark.remote_data(source="astropy")
def test_download_cache():
    download_dir = _get_download_cache_loc()

    # Download the test URL and make sure it exists, then clear just that
    # URL and make sure it got deleted.
    fnout = download_file(TESTURL, cache=True)
    assert os.path.isdir(download_dir)
    assert os.path.isfile(fnout)
    clear_download_cache(TESTURL)
    assert not os.path.exists(fnout)

    # Clearing download cache succeeds even if the URL does not exist.
    clear_download_cache("http://this_was_never_downloaded_before.com")

    # Make sure lockdir was released
    lockdir = os.path.join(download_dir, "lock")
    assert not os.path.isdir(lockdir), "Cache dir lock was not released!"


@pytest.mark.remote_data(source="astropy")
def test_download_certificate_verification_failed():
    """Tests for https://github.com/astropy/astropy/pull/10434"""

    # First test the expected exception when download fails due to a
    # certificate verification error; we simulate this by passing a bogus
    # CA directory to the ssl_context argument
    ssl_context = {"cafile": None, "capath": "/does/not/exist"}
    msg = f"Verification of TLS/SSL certificate at {TESTURL_SSL} failed"
    with pytest.raises(urllib.error.URLError, match=msg):
        download_file(TESTURL_SSL, cache=False, ssl_context=ssl_context)

    with pytest.warns(AstropyWarning, match=msg) as warning_lines:
        fnout = download_file(
            TESTURL_SSL, cache=False, ssl_context=ssl_context, allow_insecure=True
        )

    assert len(warning_lines) == 1
    assert os.path.isfile(fnout)


def test_download_cache_after_clear(tmp_path, temp_cache, valid_urls):
    testurl, contents = next(valid_urls)
    # Test issues raised in #4427 with clear_download_cache() without a URL,
    # followed by subsequent download.
    download_dir = _get_download_cache_loc()

    fnout = download_file(testurl, cache=True)
    assert os.path.isfile(fnout)

    clear_download_cache()
    assert not os.path.exists(fnout)
    assert not os.path.exists(download_dir)

    fnout = download_file(testurl, cache=True)
    assert os.path.isfile(fnout)


@pytest.mark.remote_data(source="astropy")
def test_download_parallel_from_internet_works(temp_cache):
    main_url = conf.dataurl
    mirror_url = conf.dataurl_mirror
    fileloc = "intersphinx/README"
    urls = []
    sources = {}
    for s in ["", fileloc]:
        urls.append(main_url + s)
        sources[urls[-1]] = [urls[-1], mirror_url + s]
    fnout = download_files_in_parallel(urls, sources=sources)
    assert all(os.path.isfile(f) for f in fnout), fnout


@pytest.mark.parametrize("method", [None, "spawn"])
def test_download_parallel_fills_cache(tmp_path, valid_urls, method):
    urls = []
    # tmp_path is shared between many tests, and that can cause weird
    # interactions if we set the temporary cache too directly
    with paths.set_temp_cache(tmp_path):
        for um, c in islice(valid_urls, FEW):
            assert not is_url_in_cache(um)
            urls.append((um, c))
        rs = download_files_in_parallel(
            [u for (u, c) in urls], multiprocessing_start_method=method
        )
        assert len(rs) == len(urls)
        url_set = {u for (u, c) in urls}
        assert url_set <= set(get_cached_urls())
        for r, (u, c) in zip(rs, urls):
            assert get_file_contents(r) == c
        check_download_cache()
    assert not url_set.intersection(get_cached_urls())
    check_download_cache()


def test_download_parallel_with_empty_sources(valid_urls, temp_cache):
    urls = []
    sources = {}
    for um, c in islice(valid_urls, FEW):
        assert not is_url_in_cache(um)
        urls.append((um, c))
    rs = download_files_in_parallel([u for (u, c) in urls], sources=sources)
    assert len(rs) == len(urls)
    # u = set(u for (u, c) in urls)
    # assert u <= set(get_cached_urls())
    check_download_cache()
    for r, (u, c) in zip(rs, urls):
        assert get_file_contents(r) == c


def test_download_parallel_with_sources_and_bogus_original(
    valid_urls, invalid_urls, temp_cache
):
    u, c = next(valid_urls)
    urls = [(u, c, None)]
    sources = {}
    for i, (um, c_bad) in enumerate(islice(valid_urls, FEW)):
        assert not is_url_in_cache(um)
        sources[um] = list(islice(invalid_urls, i))
        u, c = next(valid_urls)
        sources[um].append(u)
        urls.append((um, c, c_bad))
    rs = download_files_in_parallel([u for (u, c, c_bad) in urls], sources=sources)
    assert len(rs) == len(urls)
    # u = set(u for (u, c, c_bad) in urls)
    # assert u <= set(get_cached_urls())
    for r, (u, c, c_bad) in zip(rs, urls):
        assert get_file_contents(r) == c
        assert get_file_contents(r) != c_bad


def test_download_parallel_many(temp_cache, valid_urls):
    td = list(islice(valid_urls, N_PARALLEL_HAMMER))

    r = download_files_in_parallel([u for (u, c) in td])
    assert len(r) == len(td)
    for r_, (_, c) in zip(r, td):
        assert get_file_contents(r_) == c


def test_download_parallel_partial_success(temp_cache, valid_urls, invalid_urls):
    """Check that a partially successful download works.

    Even in the presence of many requested URLs, presumably hitting all the
    parallelism this system can manage, a download failure leads to a tidy
    shutdown.

    """
    td = list(islice(valid_urls, N_PARALLEL_HAMMER))

    u_bad = next(invalid_urls)

    with pytest.raises(urllib.request.URLError):
        download_files_in_parallel([u_bad] + [u for (u, c) in td])
    # Actually some files may get downloaded, others not.
    # Is this good? Should we stubbornly keep trying?
    # assert not any([is_url_in_cache(u) for (u, c) in td])


@pytest.mark.slow
def test_download_parallel_partial_success_lock_safe(
    temp_cache, valid_urls, invalid_urls
):
    """Check that a partially successful parallel download leaves the cache unlocked.

    This needs to be repeated many times because race conditions are what cause
    this sort of thing, especially situations where a process might be forcibly
    shut down while it holds the lock.

    """
    s = random.getstate()
    try:
        random.seed(0)
        for _ in range(N_PARALLEL_HAMMER):
            td = list(islice(valid_urls, FEW))

            u_bad = next(invalid_urls)
            urls = [u_bad] + [u for (u, c) in td]
            random.shuffle(urls)
            with pytest.raises(urllib.request.URLError):
                download_files_in_parallel(urls)
    finally:
        random.setstate(s)


def test_download_parallel_update(temp_cache, tmp_path):
    td = []
    for i in range(N_PARALLEL_HAMMER):
        c = f"{i:04d}"
        fn = tmp_path / c
        with open(fn, "w") as f:
            f.write(c)
        u = url_to(fn)
        clear_download_cache(u)
        td.append((fn, u, c))

    r1 = download_files_in_parallel([u for (fn, u, c) in td])
    assert len(r1) == len(td)
    for r_1, (fn, u, c) in zip(r1, td):
        assert get_file_contents(r_1) == c

    td2 = []
    for fn, u, c in td:
        c_plus = f"{c} updated"
        fn = tmp_path / c
        with open(fn, "w") as f:
            f.write(c_plus)
        td2.append((fn, u, c, c_plus))

    r2 = download_files_in_parallel([u for (fn, u, c) in td], cache=True)
    assert len(r2) == len(td)
    for r_2, (fn, u, c, c_plus) in zip(r2, td2):
        assert get_file_contents(r_2) == c
        assert c != c_plus
    r3 = download_files_in_parallel([u for (fn, u, c) in td], cache="update")

    assert len(r3) == len(td)
    for r_3, (fn, u, c, c_plus) in zip(r3, td2):
        assert get_file_contents(r_3) != c
        assert get_file_contents(r_3) == c_plus


@pytest.mark.skipif(
    (sys.platform.startswith("win") and CI), reason="flaky cache error on Windows CI"
)
def test_update_parallel(temp_cache, valid_urls):
    u, c = next(valid_urls)
    u2, c2 = next(valid_urls)

    f = download_file(u, cache=True)
    assert get_file_contents(f) == c

    def update(i):
        return download_file(u, cache="update", sources=[u2])

    with ThreadPoolExecutor(max_workers=N_THREAD_HAMMER) as P:
        r = set(P.map(update, range(N_THREAD_HAMMER)))

    check_download_cache()
    for f in r:
        assert get_file_contents(f) == c2


@pytest.mark.skipif(
    (sys.platform.startswith("win") and CI), reason="flaky cache error on Windows CI"
)
def test_update_parallel_multi(temp_cache, valid_urls):
    u, c = next(valid_urls)
    iucs = list(islice(valid_urls, N_THREAD_HAMMER))

    f = download_file(u, cache=True)
    assert get_file_contents(f) == c

    def update(uc):
        u2, c2 = uc
        return download_file(u, cache="update", sources=[u2]), c2

    with ThreadPoolExecutor(max_workers=len(iucs)) as P:
        r = list(P.map(update, iucs))

    check_download_cache()
    assert any(get_file_contents(f) == c for (f, c) in r)


@pytest.mark.remote_data(source="astropy")
def test_url_nocache():
    with get_readable_fileobj(TESTURL, cache=False, encoding="utf-8") as page:
        assert page.read().find("Astropy") > -1


def test_find_by_hash(valid_urls, temp_cache):
    testurl, contents = next(valid_urls)
    p = download_file(testurl, cache=True)
    hash = compute_hash(p)

    hashstr = "hash/" + hash

    fnout = get_pkg_data_filename(hashstr)
    assert os.path.isfile(fnout)
    clear_download_cache(fnout)
    assert not os.path.isfile(fnout)


@pytest.mark.remote_data(source="astropy")
def test_find_invalid():
    # this is of course not a real data file and not on any remote server, but
    # it should *try* to go to the remote server
    with pytest.raises((urllib.error.URLError, TimeoutError)):
        get_pkg_data_filename(
            "kjfrhgjklahgiulrhgiuraehgiurhgiuhreglhurieghruelighiuerahiulruli"
        )


@pytest.mark.parametrize("package", [None, "astropy", "numpy"])
def test_get_invalid(package):
    """Test can create a file path to an invalid file."""
    path = get_pkg_data_path("kjfrhgjkla", "hgiulrhgiu", package=package)
    assert not os.path.isfile(path)
    assert not os.path.isdir(path)


# Package data functions
@pytest.mark.parametrize(
    "filename", ["local.dat", "local.dat.gz", "local.dat.bz2", "local.dat.xz"]
)
def test_local_data_obj(filename):
    if (not HAS_BZ2 and "bz2" in filename) or (not HAS_LZMA and "xz" in filename):
        with pytest.raises(ValueError, match=r" format files are not supported"):
            with get_pkg_data_fileobj(
                os.path.join("data", filename), encoding="binary"
            ) as f:
                f.readline()
                # assert f.read().rstrip() == b'CONTENT'
    else:
        with get_pkg_data_fileobj(
            os.path.join("data", filename), encoding="binary"
        ) as f:
            f.readline()
            assert f.read().rstrip() == b"CONTENT"


@pytest.fixture(params=["invalid.dat.bz2", "invalid.dat.gz"])
def bad_compressed(request, tmp_path):
    # These contents have valid headers for their respective file formats, but
    # are otherwise malformed and invalid.
    bz_content = b"BZhinvalid"
    gz_content = b"\x1f\x8b\x08invalid"

    datafile = tmp_path / request.param
    filename = str(datafile)

    if filename.endswith(".bz2"):
        contents = bz_content
    elif filename.endswith(".gz"):
        contents = gz_content
    else:
        contents = "invalid"

    datafile.write_bytes(contents)

    return filename


def test_local_data_obj_invalid(bad_compressed):
    is_bz2 = bad_compressed.endswith(".bz2")
    is_xz = bad_compressed.endswith(".xz")

    # Note, since these invalid files are created on the fly in order to avoid
    # problems with detection by antivirus software
    # (see https://github.com/astropy/astropy/issues/6520), it is no longer
    # possible to use ``get_pkg_data_fileobj`` to read the files. Technically,
    # they're not local anymore: they just live in a temporary directory
    # created by pytest. However, we can still use get_readable_fileobj for the
    # test.
    if (not HAS_BZ2 and is_bz2) or (not HAS_LZMA and is_xz):
        with pytest.raises(
            ModuleNotFoundError, match=r"does not provide the [lb]z[2m]a? module\."
        ):
            with get_readable_fileobj(bad_compressed, encoding="binary") as f:
                f.read()
    else:
        with get_readable_fileobj(bad_compressed, encoding="binary") as f:
            assert f.read().rstrip().endswith(b"invalid")


def test_local_data_name():
    assert os.path.isfile(TESTLOCAL) and TESTLOCAL.endswith("local.dat")

    # TODO: if in the future, the root data/ directory is added in, the below
    # test should be uncommented and the README.rst should be replaced with
    # whatever file is there

    # get something in the astropy root
    # fnout2 = get_pkg_data_filename('../../data/README.rst')
    # assert os.path.isfile(fnout2) and fnout2.endswith('README.rst')


def test_data_name_third_party_package():
    """Regression test for issue #1256

    Tests that `get_pkg_data_filename` works in a third-party package that
    doesn't make any relative imports from the module it's used from.

    Uses a test package under ``data/test_package``.
    """

    # Get the actual data dir:
    data_dir = os.path.join(os.path.dirname(__file__), "data")

    sys.path.insert(0, data_dir)
    try:
        import test_package

        filename = test_package.get_data_filename()
        assert os.path.normcase(filename) == (
            os.path.normcase(os.path.join(data_dir, "test_package", "data", "foo.txt"))
        )
    finally:
        sys.path.pop(0)


def test_local_data_nonlocalfail():
    # this would go *outside* the astropy tree
    with pytest.raises(RuntimeError):
        get_pkg_data_filename("../../../data/README.rst")


def test_compute_hash(tmp_path):
    rands = b"1234567890abcdefghijklmnopqrstuvwxyz"

    filename = tmp_path / "tmp.dat"

    with open(filename, "wb") as ntf:
        ntf.write(rands)
        ntf.flush()

    chhash = compute_hash(filename)
    shash = hashlib.md5(rands).hexdigest()

    assert chhash == shash


def test_get_pkg_data_contents():
    with get_pkg_data_fileobj("data/local.dat") as f:
        contents1 = f.read()

    contents2 = get_pkg_data_contents("data/local.dat")

    assert contents1 == contents2


@pytest.mark.remote_data(source="astropy")
def test_data_noastropy_fallback(monkeypatch):
    """
    Tests to make sure the default behavior when the cache directory can't
    be located is correct
    """

    # better yet, set the configuration to make sure the temp files are deleted
    conf.delete_temporary_downloads_at_exit = True

    # make sure the config and cache directories are not searched
    monkeypatch.setenv("XDG_CONFIG_HOME", "foo")
    monkeypatch.delenv("XDG_CONFIG_HOME")
    monkeypatch.setenv("XDG_CACHE_HOME", "bar")
    monkeypatch.delenv("XDG_CACHE_HOME")

    monkeypatch.setattr(paths.set_temp_config, "_temp_path", None)
    monkeypatch.setattr(paths.set_temp_cache, "_temp_path", None)

    # make sure the _find_or_create_astropy_dir function fails as though the
    # astropy dir could not be accessed
    def osraiser(dirnm, linkto, pkgname=None):
        raise OSError()

    monkeypatch.setattr(paths, "_find_or_create_root_dir", osraiser)

    with pytest.raises(OSError):
        # make sure the config dir search fails
        paths.get_cache_dir(rootname="astropy")

    with pytest.warns(CacheMissingWarning) as warning_lines:
        fnout = download_file(TESTURL, cache=True)
    n_warns = len(warning_lines)

    partial_warn_msgs = ["remote data cache could not be accessed", "temporary file"]
    if n_warns == 4:
        partial_warn_msgs.extend(["socket", "socket"])

    for wl in warning_lines:
        cur_w = str(wl).lower()
        for i, partial_msg in enumerate(partial_warn_msgs):
            if partial_msg in cur_w:
                del partial_warn_msgs[i]
                break
    assert (
        len(partial_warn_msgs) == 0
    ), f"Got some unexpected warnings: {partial_warn_msgs}"

    assert n_warns in (2, 4), f"Expected 2 or 4 warnings, got {n_warns}"

    assert os.path.isfile(fnout)

    # clearing the cache should be a no-up that doesn't affect fnout
    with pytest.warns(CacheMissingWarning) as record:
        clear_download_cache(TESTURL)
    assert len(record) == 2
    assert (
        record[0].message.args[0]
        == "Remote data cache could not be accessed due to OSError"
    )
    assert "Not clearing data cache - cache inaccessible" in record[1].message.args[0]
    assert os.path.isfile(fnout)

    # now remove it so tests don't clutter up the temp dir this should get
    # called at exit, anyway, but we do it here just to make sure it's working
    # correctly
    _deltemps()
    assert not os.path.isfile(fnout)

    # now try with no cache
    fnnocache = download_file(TESTURL, cache=False)
    with open(fnnocache, "rb") as page:
        assert page.read().decode("utf-8").find("Astropy") > -1

    # no warnings should be raise in fileobj because cache is unnecessary


@pytest.mark.parametrize(
    "filename",
    [
        "unicode.txt",
        "unicode.txt.gz",
        pytest.param(
            "unicode.txt.bz2",
            marks=pytest.mark.xfail(not HAS_BZ2, reason="no bz2 support"),
        ),
        pytest.param(
            "unicode.txt.xz",
            marks=pytest.mark.xfail(not HAS_LZMA, reason="no lzma support"),
        ),
    ],
)
def test_read_unicode(filename):
    contents = get_pkg_data_contents(os.path.join("data", filename), encoding="utf-8")
    assert isinstance(contents, str)
    contents = contents.splitlines()[1]
    assert contents == "האסטרונומי פייתון"

    contents = get_pkg_data_contents(os.path.join("data", filename), encoding="binary")
    assert isinstance(contents, bytes)
    x = contents.splitlines()[1]
    expected = (
        b"\xff\xd7\x94\xd7\x90\xd7\xa1\xd7\x98\xd7\xa8\xd7\x95\xd7\xa0\xd7\x95"
        b"\xd7\x9e\xd7\x99 \xd7\xa4\xd7\x99\xd7\x99\xd7\xaa\xd7\x95\xd7\x9f"[1:]
    )
    assert x == expected


def test_compressed_stream():
    gzipped_data = (
        b"H4sICIxwG1AAA2xvY2FsLmRhdAALycgsVkjLzElVANKlxakpCpl5CiUZqQ"
        b"olqcUl8Tn5yYk58SmJJYnxWmCRzLx0hbTSvOSSzPy8Yi5nf78QV78QLgAlLytnRQAAAA=="
    )
    gzipped_data = base64.b64decode(gzipped_data)
    assert isinstance(gzipped_data, bytes)

    class FakeStream:
        """
        A fake stream that has `read`, but no `seek`.
        """

        def __init__(self, data):
            self.data = data

        def read(self, nbytes=None):
            if nbytes is None:
                result = self.data
                self.data = b""
            else:
                result = self.data[:nbytes]
                self.data = self.data[nbytes:]
            return result

    stream = FakeStream(gzipped_data)
    with get_readable_fileobj(stream, encoding="binary") as f:
        f.readline()
        assert f.read().rstrip() == b"CONTENT"


@pytest.mark.remote_data(source="astropy")
def test_invalid_location_download_raises_urlerror():
    """
    checks that download_file gives a URLError and not an AttributeError,
    as its code pathway involves some fiddling with the exception.
    """

    with pytest.raises(urllib.error.URLError):
        download_file("http://www.astropy.org/nonexistentfile")


def test_invalid_location_download_noconnect():
    """
    checks that download_file gives an OSError if the socket is blocked
    """

    # This should invoke socket's monkeypatched failure
    with pytest.raises(OSError):
        download_file("http://astropy.org/nonexistentfile")


@pytest.mark.remote_data(source="astropy")
def test_is_url_in_cache_remote():
    assert not is_url_in_cache("http://astropy.org/nonexistentfile")

    download_file(TESTURL, cache=True, show_progress=False)
    assert is_url_in_cache(TESTURL)


def test_is_url_in_cache_local(temp_cache, valid_urls, invalid_urls):
    testurl, contents = next(valid_urls)
    nonexistent = next(invalid_urls)

    assert not is_url_in_cache(testurl)
    assert not is_url_in_cache(nonexistent)

    download_file(testurl, cache=True, show_progress=False)
    assert is_url_in_cache(testurl)
    assert not is_url_in_cache(nonexistent)


# If non-deterministic failure happens see
# https://github.com/astropy/astropy/issues/9765
def test_check_download_cache(tmp_path, temp_cache, valid_urls, invalid_urls):
    testurl, testurl_contents = next(valid_urls)
    testurl2, testurl2_contents = next(valid_urls)

    zip_file_name = tmp_path / "the.zip"
    clear_download_cache()
    assert not check_download_cache()

    download_file(testurl, cache=True)
    check_download_cache()
    download_file(testurl2, cache=True)
    check_download_cache()

    export_download_cache(zip_file_name, [testurl, testurl2])
    check_download_cache()

    clear_download_cache(testurl2)
    check_download_cache()

    import_download_cache(zip_file_name, [testurl])
    check_download_cache()


def test_export_import_roundtrip_one(tmp_path, temp_cache, valid_urls):
    testurl, contents = next(valid_urls)
    f = download_file(testurl, cache=True, show_progress=False)
    assert get_file_contents(f) == contents

    initial_urls_in_cache = set(get_cached_urls())
    zip_file_name = tmp_path / "the.zip"

    export_download_cache(zip_file_name, [testurl])
    clear_download_cache(testurl)
    import_download_cache(zip_file_name)
    assert is_url_in_cache(testurl)
    assert set(get_cached_urls()) == initial_urls_in_cache
    assert (
        get_file_contents(download_file(testurl, cache=True, show_progress=False))
        == contents
    )


def test_export_url_not_present(temp_cache, valid_urls):
    testurl, contents = next(valid_urls)
    with NamedTemporaryFile("wb") as zip_file:
        assert not is_url_in_cache(testurl)
        with pytest.raises(KeyError):
            export_download_cache(zip_file, [testurl])


def test_import_one(tmp_path, temp_cache, valid_urls):
    testurl, testurl_contents = next(valid_urls)
    testurl2, testurl2_contents = next(valid_urls)
    zip_file_name = tmp_path / "the.zip"

    download_file(testurl, cache=True)
    download_file(testurl2, cache=True)
    assert is_url_in_cache(testurl2)
    export_download_cache(zip_file_name, [testurl, testurl2])
    clear_download_cache(testurl)
    clear_download_cache(testurl2)
    import_download_cache(zip_file_name, [testurl])
    assert is_url_in_cache(testurl)
    assert not is_url_in_cache(testurl2)


def test_export_import_roundtrip(tmp_path, temp_cache, valid_urls):
    zip_file_name = tmp_path / "the.zip"
    for u, c in islice(valid_urls, FEW):
        download_file(u, cache=True)

    initial_urls_in_cache = set(get_cached_urls())

    export_download_cache(zip_file_name)
    clear_download_cache()
    import_download_cache(zip_file_name)

    assert set(get_cached_urls()) == initial_urls_in_cache


def test_export_import_roundtrip_stream(temp_cache, valid_urls):
    for u, c in islice(valid_urls, FEW):
        download_file(u, cache=True)
    initial_urls_in_cache = set(get_cached_urls())

    with io.BytesIO() as f:
        export_download_cache(f)
        b = f.getvalue()
    clear_download_cache()
    with io.BytesIO(b) as f:
        import_download_cache(f)

    assert set(get_cached_urls()) == initial_urls_in_cache


def test_export_overwrite_flag_works(temp_cache, valid_urls, tmp_path):
    fn = tmp_path / "f.zip"
    c = b"Some contents\nto check later"
    with open(fn, "wb") as f:
        f.write(c)
    for u, _ in islice(valid_urls, FEW):
        download_file(u, cache=True)

    with pytest.raises(FileExistsError):
        export_download_cache(fn)
    assert get_file_contents(fn, encoding="binary") == c

    export_download_cache(fn, overwrite=True)
    assert get_file_contents(fn, encoding="binary") != c


def test_export_import_roundtrip_different_location(tmp_path, valid_urls):
    original_cache = tmp_path / "original"
    original_cache.mkdir()
    zip_file_name = tmp_path / "the.zip"

    urls = list(islice(valid_urls, FEW))
    initial_urls_in_cache = {u for (u, c) in urls}
    with paths.set_temp_cache(original_cache):
        for u, c in urls:
            download_file(u, cache=True)
        assert set(get_cached_urls()) == initial_urls_in_cache
        export_download_cache(zip_file_name)

    new_cache = tmp_path / "new"
    new_cache.mkdir()
    with paths.set_temp_cache(new_cache):
        import_download_cache(zip_file_name)
        check_download_cache()
        assert set(get_cached_urls()) == initial_urls_in_cache
        for u, c in urls:
            assert get_file_contents(download_file(u, cache=True)) == c


def test_cache_size_is_zero_when_empty(temp_cache):
    assert not get_cached_urls()
    assert cache_total_size() == 0


def test_cache_size_changes_correctly_when_files_are_added_and_removed(
    temp_cache, valid_urls
):
    u, c = next(valid_urls)
    clear_download_cache(u)
    s_i = cache_total_size()
    download_file(u, cache=True)
    assert cache_total_size() == s_i + len(c) + len(u.encode("utf-8"))
    clear_download_cache(u)
    assert cache_total_size() == s_i


def test_cache_contents_agrees_with_get_urls(temp_cache, valid_urls):
    r = []
    for a, a_c in islice(valid_urls, FEW):
        a_f = download_file(a, cache=True)
        r.append((a, a_c, a_f))
    assert set(cache_contents().keys()) == set(get_cached_urls())
    for u, c, h in r:
        assert cache_contents()[u] == h


@pytest.mark.parametrize("desired_size", [1_000_000_000_000_000_000, 1 * _u.Ebyte])
def test_free_space_checker_huge(tmp_path, desired_size):
    with pytest.raises(OSError):
        check_free_space_in_dir(tmp_path, desired_size)


def test_get_free_space_file_directory(tmp_path):
    fn = tmp_path / "file"
    with open(fn, "w"):
        pass
    with pytest.raises(OSError):
        get_free_space_in_dir(fn)

    free_space = get_free_space_in_dir(tmp_path)
    assert free_space > 0 and not hasattr(free_space, "unit")

    # TODO: If unit=True starts to auto-guess prefix, this needs updating.
    free_space = get_free_space_in_dir(tmp_path, unit=True)
    assert free_space > 0 and free_space.unit == _u.byte

    free_space = get_free_space_in_dir(tmp_path, unit=_u.Mbit)
    assert free_space > 0 and free_space.unit == _u.Mbit


def test_download_file_bogus_settings(invalid_urls, temp_cache):
    u = next(invalid_urls)
    with pytest.raises(KeyError):
        download_file(u, sources=[])


def test_download_file_local_directory(tmp_path):
    """Make sure we get a URLError rather than OSError even if it's a
    local directory."""
    with pytest.raises(urllib.request.URLError):
        download_file(url_to(tmp_path))


def test_download_file_schedules_deletion(valid_urls):
    u, c = next(valid_urls)
    f = download_file(u)
    assert f in _tempfilestodel
    # how to test deletion actually occurs?


def test_clear_download_cache_refuses_to_delete_outside_the_cache(tmp_path):
    fn = str(tmp_path / "file")
    with open(fn, "w") as f:
        f.write("content")
    assert os.path.exists(fn)
    with pytest.raises(RuntimeError):
        clear_download_cache(fn)
    assert os.path.exists(fn)


def test_check_download_cache_finds_bogus_entries(temp_cache, valid_urls):
    u, c = next(valid_urls)
    download_file(u, cache=True)
    dldir = _get_download_cache_loc()
    bf = os.path.abspath(os.path.join(dldir, "bogus"))
    with open(bf, "w") as f:
        f.write("bogus file that exists")
    with pytest.raises(CacheDamaged) as e:
        check_download_cache()
    assert bf in e.value.bad_files
    clear_download_cache()


def test_check_download_cache_finds_bogus_subentries(temp_cache, valid_urls):
    u, c = next(valid_urls)
    f = download_file(u, cache=True)
    bf = os.path.abspath(os.path.join(os.path.dirname(f), "bogus"))
    with open(bf, "w") as f:
        f.write("bogus file that exists")
    with pytest.raises(CacheDamaged) as e:
        check_download_cache()
    assert bf in e.value.bad_files
    clear_download_cache()


def test_check_download_cache_cleanup(temp_cache, valid_urls):
    u, c = next(valid_urls)
    fn = download_file(u, cache=True)
    dldir = _get_download_cache_loc()

    bf1 = os.path.abspath(os.path.join(dldir, "bogus1"))
    with open(bf1, "w") as f:
        f.write("bogus file that exists")

    bf2 = os.path.abspath(os.path.join(os.path.dirname(fn), "bogus2"))
    with open(bf2, "w") as f:
        f.write("other bogus file that exists")

    bf3 = os.path.abspath(os.path.join(dldir, "contents"))
    with open(bf3, "w") as f:
        f.write("awkwardly-named bogus file that exists")

    u2, c2 = next(valid_urls)
    f2 = download_file(u, cache=True)
    os.unlink(f2)
    bf4 = os.path.dirname(f2)

    with pytest.raises(CacheDamaged) as e:
        check_download_cache()
    assert set(e.value.bad_files) == {bf1, bf2, bf3, bf4}
    for bf in e.value.bad_files:
        clear_download_cache(bf)
    # download cache will be checked on exit


def test_download_cache_update_doesnt_damage_cache(temp_cache, valid_urls):
    u, _ = next(valid_urls)
    download_file(u, cache=True)
    download_file(u, cache="update")


def test_cache_dir_is_actually_a_file(tmp_path, valid_urls):
    """Ensure that bogus cache settings are handled sensibly.

    Because the user can specify the cache location in a config file, and
    because they might try to deduce the location by looking around at what's
    in their directory tree, and because the cache directory is actual several
    tree levels down from the directory set in the config file, it's important
    to check what happens if each of the steps in the path is wrong somehow.
    """

    def check_quietly_ignores_bogus_cache():
        """We want a broken cache to produce a warning but then astropy should
        act like there isn't a cache.
        """
        with pytest.warns(CacheMissingWarning):
            assert not get_cached_urls()
        with pytest.warns(CacheMissingWarning):
            assert not is_url_in_cache("http://www.example.com/")
        with pytest.warns(CacheMissingWarning):
            assert not cache_contents()
        with pytest.warns(CacheMissingWarning):
            u, c = next(valid_urls)
            r = download_file(u, cache=True)
            assert get_file_contents(r) == c
            # check the filename r appears in a warning message?
            # check r is added to the delete_at_exit list?
            # in fact should there be testing of the delete_at_exit mechanism,
            # as far as that is possible?
        with pytest.warns(CacheMissingWarning):
            assert not is_url_in_cache(u)
        with pytest.warns(CacheMissingWarning):
            with pytest.raises(OSError):
                check_download_cache()

    dldir = _get_download_cache_loc()
    # set_temp_cache acts weird if it is pointed at a file (see below)
    # but we want to see what happens when the cache is pointed
    # at a file instead of a directory, so make a directory we can
    # replace later.
    fn = tmp_path / "file"
    ct = "contents\n"
    os.mkdir(fn)
    with paths.set_temp_cache(fn):
        shutil.rmtree(fn)
        with open(fn, "w") as f:
            f.write(ct)
        with pytest.raises(OSError):
            paths.get_cache_dir()
        check_quietly_ignores_bogus_cache()
    assert dldir == _get_download_cache_loc()
    assert get_file_contents(fn) == ct, "File should not be harmed."

    # See what happens when set_temp_cache is pointed at a file
    with pytest.raises(OSError):
        with paths.set_temp_cache(fn):
            pass
    assert dldir == _get_download_cache_loc()
    assert get_file_contents(str(fn)) == ct

    # Now the cache directory is normal but the subdirectory it wants
    # to make is a file
    cd = tmp_path / "astropy"
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmp_path):
        check_quietly_ignores_bogus_cache()
    assert dldir == _get_download_cache_loc()
    assert get_file_contents(cd) == ct
    os.remove(cd)

    # Ditto one level deeper
    os.makedirs(cd)
    cd = tmp_path / "astropy" / "download"
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmp_path):
        check_quietly_ignores_bogus_cache()
    assert dldir == _get_download_cache_loc()
    assert get_file_contents(cd) == ct
    os.remove(cd)

    # Ditto another level deeper
    os.makedirs(cd)
    cd = tmp_path / "astropy" / "download" / "url"
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmp_path):
        check_quietly_ignores_bogus_cache()
    assert dldir == _get_download_cache_loc()
    assert get_file_contents(cd) == ct
    os.remove(cd)


def test_get_fileobj_str(a_file):
    fn, c = a_file
    with get_readable_fileobj(str(fn)) as rf:
        assert rf.read() == c


def test_get_fileobj_pathlib(a_file):
    fn, c = a_file
    with get_readable_fileobj(pathlib.Path(fn)) as rf:
        assert rf.read() == c


def test_get_fileobj_binary(a_binary_file):
    fn, c = a_binary_file
    with get_readable_fileobj(fn, encoding="binary") as rf:
        assert rf.read() == c


def test_get_fileobj_already_open_text(a_file):
    fn, c = a_file
    with open(fn) as f:
        with get_readable_fileobj(f) as rf:
            with pytest.raises(TypeError):
                rf.read()


def test_get_fileobj_already_open_binary(a_file):
    fn, c = a_file
    with open(fn, "rb") as f:
        with get_readable_fileobj(f) as rf:
            assert rf.read() == c


def test_get_fileobj_binary_already_open_binary(a_binary_file):
    fn, c = a_binary_file
    with open(fn, "rb") as f:
        with get_readable_fileobj(f, encoding="binary") as rf:
            assert rf.read() == c


def test_cache_contents_not_writable(temp_cache, valid_urls):
    c = cache_contents()
    with pytest.raises(TypeError):
        c["foo"] = 7
    u, _ = next(valid_urls)
    download_file(u, cache=True)
    c = cache_contents()
    assert u in c
    with pytest.raises(TypeError):
        c["foo"] = 7


def test_cache_relocatable(tmp_path, valid_urls):
    u, c = next(valid_urls)
    d1 = tmp_path / "1"
    d2 = tmp_path / "2"
    os.mkdir(d1)
    with paths.set_temp_cache(d1):
        p1 = download_file(u, cache=True)
        assert is_url_in_cache(u)
        assert get_file_contents(p1) == c
        shutil.copytree(d1, d2)
        clear_download_cache()
    with paths.set_temp_cache(d2):
        assert is_url_in_cache(u)
        p2 = download_file(u, cache=True)
        assert p1 != p2
        assert os.path.exists(p2)
        clear_download_cache(p2)
        check_download_cache()


def test_get_readable_fileobj_cleans_up_temporary_files(tmp_path, monkeypatch):
    """checks that get_readable_fileobj leaves no temporary files behind"""
    # Create a 'file://' URL pointing to a path on the local filesystem
    url = url_to(TESTLOCAL)

    # Save temporary files to a known location
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))

    # Call get_readable_fileobj() as a context manager
    with get_readable_fileobj(url) as f:
        f.read()

    # Get listing of files in temporary directory
    tempdir_listing = list(tmp_path.iterdir())

    # Assert that the temporary file was empty after get_readable_fileobj()
    # context manager finished running
    assert len(tempdir_listing) == 0


def test_path_objects_get_readable_fileobj():
    fpath = pathlib.Path(TESTLOCAL)
    with get_readable_fileobj(fpath) as f:
        assert (
            f.read().rstrip()
            == "This file is used in the test_local_data_* testing functions\nCONTENT"
        )


def test_nested_get_readable_fileobj():
    """Ensure fileobj state is as expected when get_readable_fileobj()
    is called inside another get_readable_fileobj().
    """
    with get_readable_fileobj(TESTLOCAL, encoding="binary") as fileobj:
        with get_readable_fileobj(fileobj, encoding="UTF-8") as fileobj2:
            fileobj2.seek(1)
        fileobj.seek(1)

        # Theoretically, fileobj2 should be closed already here but it is not.
        # See https://github.com/astropy/astropy/pull/8675.
        # UNCOMMENT THIS WHEN PYTHON FINALLY LETS IT HAPPEN.
        # assert fileobj2.closed

    assert fileobj.closed and fileobj2.closed


def test_download_file_wrong_size(monkeypatch):
    @contextlib.contextmanager
    def mockurl(remote_url, timeout=None):
        yield MockURL()

    def mockurl_builder(*args, tlscontext=None, **kwargs):
        mock_opener = type("MockOpener", (object,), {})()
        mock_opener.open = mockurl
        return mock_opener

    class MockURL:
        def __init__(self):
            self.reader = io.BytesIO(b"a" * real_length)

        def info(self):
            return {"Content-Length": str(report_length)}

        def read(self, length=None):
            return self.reader.read(length)

    monkeypatch.setattr(astropy.utils.data, "_build_urlopener", mockurl_builder)

    with pytest.raises(urllib.error.ContentTooShortError):
        report_length = 1024
        real_length = 1023
        download_file(TESTURL, cache=False)

    with pytest.raises(urllib.error.URLError):
        report_length = 1023
        real_length = 1024
        download_file(TESTURL, cache=False)

    report_length = 1023
    real_length = 1023
    fn = download_file(TESTURL, cache=False)
    with open(fn, "rb") as f:
        assert f.read() == b"a" * real_length

    report_length = None
    real_length = 1023
    fn = download_file(TESTURL, cache=False)
    with open(fn, "rb") as f:
        assert f.read() == b"a" * real_length


def test_can_make_directories_readonly(tmp_path):
    try:
        with readonly_dir(tmp_path):
            assert is_dir_readonly(tmp_path)
    except AssertionError:
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            pytest.skip(
                "We are root, we can't make a directory un-writable with chmod."
            )
        elif platform.system() == "Windows":
            pytest.skip(
                "It seems we can't make a directory un-writable under Windows "
                "with chmod, in spite of the documentation."
            )
        else:
            raise


def test_can_make_files_readonly(tmp_path):
    fn = tmp_path / "test"
    c = "contents\n"
    with open(fn, "w") as f:
        f.write(c)
    with readonly_dir(tmp_path):
        try:
            with open(fn, "w+") as f:
                f.write("more contents\n")
        except PermissionError:
            pass
        else:
            if hasattr(os, "geteuid") and os.geteuid() == 0:
                pytest.skip("We are root, we can't make a file un-writable with chmod.")
    assert get_file_contents(fn) == c


def test_read_cache_readonly(readonly_cache):
    assert cache_contents() == readonly_cache


def test_download_file_cache_readonly(readonly_cache):
    for u in readonly_cache:
        f = download_file(u, cache=True)
        assert f == readonly_cache[u]


def test_import_file_cache_readonly(readonly_cache, tmp_path):
    filename = tmp_path / "test-file"
    content = "Some text or other"
    url = "http://example.com/"
    with open(filename, "w") as f:
        f.write(content)

    with pytest.raises(OSError):
        import_file_to_cache(url, filename, remove_original=True)
    assert not is_url_in_cache(url)


def test_import_file_cache_invalid_cross_device_link(tmp_path, monkeypatch):
    def no_rename(path, mode=None):
        if os.path.exists(path):
            raise OSError(errno.EXDEV, "os.rename monkeypatched out")
        else:
            raise FileNotFoundError(f"File {path} does not exist.")

    monkeypatch.setattr(os, "rename", no_rename)

    filename = tmp_path / "test-file"
    content = "Some text or other"
    url = "http://example.com/"
    with open(filename, "w") as f:
        f.write(content)

    with pytest.warns(AstropyWarning, match="os.rename monkeypatched out"):
        import_file_to_cache(url, filename, remove_original=True, replace=True)
    assert is_url_in_cache(url)


def test_download_file_cache_readonly_cache_miss(readonly_cache, valid_urls):
    u, c = next(valid_urls)
    with pytest.warns(CacheMissingWarning):
        f = download_file(u, cache=True)
    assert get_file_contents(f) == c
    assert not is_url_in_cache(u)


def test_download_file_cache_readonly_update(readonly_cache):
    for u in readonly_cache:
        with pytest.warns(CacheMissingWarning):
            f = download_file(u, cache="update")
        assert f != readonly_cache[u]
        assert compute_hash(f) == compute_hash(readonly_cache[u])


def test_check_download_cache_works_if_readonly(readonly_cache):
    check_download_cache()


# On Windows I can't make directories readonly. On CircleCI I can't make
# anything readonly because the test suite runs as root. So on those platforms
# none of the "real" tests above can be run. I can use monkeypatch to trigger
# the readonly code paths, see the "fake" versions of the tests below, but I
# don't totally trust those to completely explore what happens either, so we
# have both. I couldn't see an easy way to parameterize over fixtures and share
# tests.


def test_read_cache_fake_readonly(fake_readonly_cache):
    assert cache_contents() == fake_readonly_cache


def test_download_file_cache_fake_readonly(fake_readonly_cache):
    for u in fake_readonly_cache:
        f = download_file(u, cache=True)
        assert f == fake_readonly_cache[u]


def test_mkdtemp_cache_fake_readonly(fake_readonly_cache):
    with pytest.raises(OSError):
        tempfile.mkdtemp()


def test_TD_cache_fake_readonly(fake_readonly_cache):
    with pytest.raises(OSError):
        with TemporaryDirectory():
            pass


def test_import_file_cache_fake_readonly(fake_readonly_cache, tmp_path):
    filename = tmp_path / "test-file"
    content = "Some text or other"
    url = "http://example.com/"
    with open(filename, "w") as f:
        f.write(content)

    with pytest.raises(OSError):
        import_file_to_cache(url, filename, remove_original=True)
    assert not is_url_in_cache(url)


def test_download_file_cache_fake_readonly_cache_miss(fake_readonly_cache, valid_urls):
    u, c = next(valid_urls)
    with pytest.warns(CacheMissingWarning):
        f = download_file(u, cache=True)
    assert not is_url_in_cache(u)
    assert get_file_contents(f) == c


def test_download_file_cache_fake_readonly_update(fake_readonly_cache):
    for u in fake_readonly_cache:
        with pytest.warns(CacheMissingWarning):
            f = download_file(u, cache="update")
        assert f != fake_readonly_cache[u]
        assert compute_hash(f) == compute_hash(fake_readonly_cache[u])


def test_check_download_cache_works_if_fake_readonly(fake_readonly_cache):
    check_download_cache()


def test_pkgname_isolation(temp_cache, valid_urls):
    a = "bogus_cache_name"

    assert not get_cached_urls()
    assert not get_cached_urls(pkgname=a)

    for u, _ in islice(valid_urls, FEW):
        download_file(u, cache=True, pkgname=a)
    assert not get_cached_urls()
    assert len(get_cached_urls(pkgname=a)) == FEW
    assert cache_total_size() < cache_total_size(pkgname=a)

    for u, _ in islice(valid_urls, FEW + 1):
        download_file(u, cache=True)
    assert len(get_cached_urls()) == FEW + 1
    assert len(get_cached_urls(pkgname=a)) == FEW
    assert cache_total_size() > cache_total_size(pkgname=a)

    assert set(get_cached_urls()) == set(cache_contents().keys())
    assert set(get_cached_urls(pkgname=a)) == set(cache_contents(pkgname=a).keys())
    for i in get_cached_urls():
        assert is_url_in_cache(i)
        assert not is_url_in_cache(i, pkgname=a)
    for i in get_cached_urls(pkgname=a):
        assert not is_url_in_cache(i)
        assert is_url_in_cache(i, pkgname=a)

    # FIXME: need to break a cache to test whether we check the right one
    check_download_cache()
    check_download_cache(pkgname=a)

    # FIXME: check that cache='update' works

    u = get_cached_urls()[0]
    with pytest.raises(KeyError):
        download_file(u, cache=True, sources=[], pkgname=a)
    clear_download_cache(u, pkgname=a)
    assert len(get_cached_urls()) == FEW + 1, "wrong pkgname should do nothing"
    assert len(get_cached_urls(pkgname=a)) == FEW, "wrong pkgname should do nothing"

    f = download_file(u, sources=[], cache=True)
    with pytest.raises(RuntimeError):
        clear_download_cache(f, pkgname=a)

    ua = get_cached_urls(pkgname=a)[0]
    with pytest.raises(KeyError):
        download_file(ua, cache=True, sources=[])

    fa = download_file(ua, sources=[], cache=True, pkgname=a)
    with pytest.raises(RuntimeError):
        clear_download_cache(fa)

    clear_download_cache(ua, pkgname=a)
    assert len(get_cached_urls()) == FEW + 1
    assert len(get_cached_urls(pkgname=a)) == FEW - 1

    clear_download_cache(u)
    assert len(get_cached_urls()) == FEW
    assert len(get_cached_urls(pkgname=a)) == FEW - 1

    clear_download_cache(pkgname=a)
    assert len(get_cached_urls()) == FEW
    assert not get_cached_urls(pkgname=a)

    clear_download_cache()
    assert not get_cached_urls()
    assert not get_cached_urls(pkgname=a)


def test_transport_cache_via_zip(temp_cache, valid_urls):
    a = "bogus_cache_name"

    assert not get_cached_urls()
    assert not get_cached_urls(pkgname=a)

    for u, _ in islice(valid_urls, FEW):
        download_file(u, cache=True)

    with io.BytesIO() as f:
        export_download_cache(f)
        b = f.getvalue()
    with io.BytesIO(b) as f:
        import_download_cache(f, pkgname=a)

    check_download_cache()
    check_download_cache(pkgname=a)

    assert set(get_cached_urls()) == set(get_cached_urls(pkgname=a))
    cca = cache_contents(pkgname=a)
    for k, v in cache_contents().items():
        assert v != cca[k]
        assert get_file_contents(v) == get_file_contents(cca[k])
    clear_download_cache()

    with io.BytesIO() as f:
        export_download_cache(f, pkgname=a)
        b = f.getvalue()
    with io.BytesIO(b) as f:
        import_download_cache(f)

    assert set(get_cached_urls()) == set(get_cached_urls(pkgname=a))


def test_download_parallel_respects_pkgname(temp_cache, valid_urls):
    a = "bogus_cache_name"

    assert not get_cached_urls()
    assert not get_cached_urls(pkgname=a)

    download_files_in_parallel([u for (u, c) in islice(valid_urls, FEW)], pkgname=a)
    assert not get_cached_urls()
    assert len(get_cached_urls(pkgname=a)) == FEW


@pytest.mark.skipif(
    not CAN_RENAME_DIRECTORY_IN_USE,
    reason="This platform is unable to rename directories that are in use.",
)
def test_removal_of_open_files(temp_cache, valid_urls):
    u, c = next(valid_urls)
    with open(download_file(u, cache=True)):
        clear_download_cache(u)
        assert not is_url_in_cache(u)
        check_download_cache()


@pytest.mark.skipif(
    not CAN_RENAME_DIRECTORY_IN_USE,
    reason="This platform is unable to rename directories that are in use.",
)
def test_update_of_open_files(temp_cache, valid_urls):
    u, c = next(valid_urls)
    with open(download_file(u, cache=True)):
        u2, c2 = next(valid_urls)
        f = download_file(u, cache="update", sources=[u2])
        check_download_cache()
        assert is_url_in_cache(u)
        assert get_file_contents(f) == c2
    assert is_url_in_cache(u)


def test_removal_of_open_files_windows(temp_cache, valid_urls, monkeypatch):
    def no_rmtree(*args, **kwargs):
        warnings.warn(CacheMissingWarning("in use"))
        raise PermissionError

    if CAN_RENAME_DIRECTORY_IN_USE:
        # This platform is able to remove files while in use.
        monkeypatch.setattr(astropy.utils.data, "_rmtree", no_rmtree)

    if PYTEST_LT_8_0:
        ctx = nullcontext()
    else:
        ctx = pytest.warns(CacheMissingWarning, match=".*PermissionError.*")

    u, c = next(valid_urls)
    with open(download_file(u, cache=True)):
        with pytest.warns(CacheMissingWarning, match=".*in use.*"), ctx:
            clear_download_cache(u)


def test_update_of_open_files_windows(temp_cache, valid_urls, monkeypatch):
    def no_rmtree(*args, **kwargs):
        warnings.warn(CacheMissingWarning("in use"))
        raise PermissionError

    if CAN_RENAME_DIRECTORY_IN_USE:
        # This platform is able to remove files while in use.
        monkeypatch.setattr(astropy.utils.data, "_rmtree", no_rmtree)

    if PYTEST_LT_8_0:
        ctx = nullcontext()
    else:
        ctx = pytest.warns(CacheMissingWarning, match=".*read-only.*")

    u, c = next(valid_urls)
    with open(download_file(u, cache=True)):
        u2, c2 = next(valid_urls)
        with pytest.warns(CacheMissingWarning, match=".*in use.*"), ctx:
            f = download_file(u, cache="update", sources=[u2])
        check_download_cache()
        assert is_url_in_cache(u)
        assert get_file_contents(f) == c2
    assert get_file_contents(download_file(u, cache=True, sources=[])) == c


def test_no_allow_internet(temp_cache, valid_urls):
    u, c = next(valid_urls)
    with conf.set_temp("allow_internet", False):
        with pytest.raises(urllib.error.URLError):
            download_file(u)
        assert not is_url_in_cache(u)
        with pytest.raises(urllib.error.URLError):
            # This will trigger the remote data error if it's allowed to touch the internet
            download_file(TESTURL)


def test_clear_download_cache_not_too_aggressive(temp_cache, valid_urls):
    u, c = next(valid_urls)
    download_file(u, cache=True)
    dldir = _get_download_cache_loc()

    bad_filename = os.path.join(dldir, "contents")
    assert is_url_in_cache(u)
    clear_download_cache(bad_filename)
    assert is_url_in_cache(u)


def test_clear_download_cache_variants(temp_cache, valid_urls):
    # deletion by contents filename
    u, c = next(valid_urls)
    f = download_file(u, cache=True)
    clear_download_cache(f)
    assert not is_url_in_cache(u)

    # deletion by url filename
    u, c = next(valid_urls)
    f = download_file(u, cache=True)
    clear_download_cache(os.path.join(os.path.dirname(f), "url"))
    assert not is_url_in_cache(u)

    # deletion by hash directory name
    u, c = next(valid_urls)
    f = download_file(u, cache=True)
    clear_download_cache(os.path.dirname(f))
    assert not is_url_in_cache(u)

    # deletion by directory name with trailing slash
    u, c = next(valid_urls)
    f = download_file(u, cache=True)
    clear_download_cache(os.path.dirname(f) + "/")
    assert not is_url_in_cache(u)

    # deletion by hash of file contents
    u, c = next(valid_urls)
    f = download_file(u, cache=True)
    h = compute_hash(f)
    clear_download_cache(h)
    assert not is_url_in_cache(u)


def test_clear_download_cache_invalid_cross_device_link(
    temp_cache, valid_urls, monkeypatch
):
    def no_rename(path, mode=None):
        raise OSError(errno.EXDEV, "os.rename monkeypatched out")

    u, c = next(valid_urls)
    download_file(u, cache=True)

    monkeypatch.setattr(os, "rename", no_rename)

    assert is_url_in_cache(u)
    with pytest.warns(AstropyWarning, match="os.rename monkeypatched out"):
        clear_download_cache(u)
    assert not is_url_in_cache(u)


def test_clear_download_cache_raises_os_error(temp_cache, valid_urls, monkeypatch):
    def no_rename(path, mode=None):
        raise OSError(errno.EBUSY, "os.rename monkeypatched out")

    u, c = next(valid_urls)
    download_file(u, cache=True)

    monkeypatch.setattr(os, "rename", no_rename)

    assert is_url_in_cache(u)
    with pytest.warns(CacheMissingWarning, match="os.rename monkeypatched out"):
        clear_download_cache(u)


@pytest.mark.skipif(
    CI and not IS_CRON,
    reason="Flaky/too much external traffic for regular CI",
)
@pytest.mark.remote_data
def test_ftp_tls_auto(temp_cache):
    """Test that download automatically enables TLS/SSL when required"""

    url = "ftp://anonymous:mail%40astropy.org@gdc.cddis.eosdis.nasa.gov/pub/products/iers/finals2000A.daily"
    download_file(url)


@pytest.mark.parametrize("base", ["http://example.com", "https://example.com"])
def test_url_trailing_slash(temp_cache, valid_urls, base):
    slash = base + "/"
    no_slash = base

    u, c = next(valid_urls)

    download_file(slash, cache=True, sources=[u])

    assert is_url_in_cache(no_slash)
    download_file(no_slash, cache=True, sources=[])
    clear_download_cache(no_slash)
    assert not is_url_in_cache(no_slash)
    assert not is_url_in_cache(slash)

    download_file(no_slash, cache=True, sources=[u])
    # see if implicit check_download_cache squawks


def test_empty_url(temp_cache, valid_urls):
    u, c = next(valid_urls)
    download_file("file://", cache=True, sources=[u])
    assert not is_url_in_cache("file:///")


@pytest.mark.remote_data
def test_download_ftp_file_properly_handles_socket_error():
    faulty_url = "ftp://anonymous:mail%40astropy.org@nonexisting/pub/products/iers/finals2000A.all"
    with pytest.raises(urllib.error.URLError) as excinfo:
        download_file(faulty_url)
    errmsg = excinfo.exconly()
    found_msg = False
    possible_msgs = [
        "Name or service not known",
        "nodename nor servname provided, or not known",
        "getaddrinfo failed",
        "Temporary failure in name resolution",
        "No address associated with hostname",
    ]
    for cur_msg in possible_msgs:
        if cur_msg in errmsg:
            found_msg = True
            break
    assert found_msg, f'Got {errmsg}, expected one of these: {",".join(possible_msgs)}'


@pytest.mark.parametrize(
    ("s", "ans"),
    [
        ("http://googlecom", True),
        ("https://google.com", True),
        ("ftp://google.com", True),
        ("sftp://google.com", True),
        ("ssh://google.com", True),
        ("file:///c:/path/to/the%20file.txt", True),
        ("google.com", False),
        ("C:\\\\path\\\\file.docx", False),
        ("data://file", False),
    ],
)
def test_string_is_url_check(s, ans):
    assert is_url(s) is ans
