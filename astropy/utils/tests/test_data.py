# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
import dbm
import sys
import stat
import base64
import random
import shutil
import hashlib
import pathlib
import platform
import tempfile
import importlib
import itertools
import contextlib
import urllib.error
import urllib.parse
import urllib.request
from itertools import islice
from concurrent.futures import ThreadPoolExecutor
from tempfile import NamedTemporaryFile, TemporaryDirectory

import py.path
import pytest

from astropy.utils import data
from astropy.config import paths
from astropy.utils.compat.context import nullcontext
from astropy.utils.data import (
    CacheMissingWarning,
    CacheDamaged,
    conf,
    _cache,
    compute_hash,
    download_file,
    cache_contents,
    _tempfilestodel,
    get_cached_urls,
    is_url_in_cache,
    cache_total_size,
    get_file_contents,
    check_download_cache,
    clear_download_cache,
    get_pkg_data_fileobj,
    get_readable_fileobj,
    export_download_cache,
    get_pkg_data_contents,
    get_pkg_data_filename,
    import_download_cache,
    get_free_space_in_dir,
    check_free_space_in_dir,
    _get_download_cache_locs,
    download_files_in_parallel,
)
from astropy.tests.helper import catch_warnings

TESTURL = "http://www.astropy.org"
TESTURL2 = "http://www.astropy.org/about.html"
TESTLOCAL = get_pkg_data_filename(os.path.join("data", "local.dat"))

try:
    import bz2  # noqa
except ImportError:
    HAS_BZ2 = False
else:
    HAS_BZ2 = True

try:
    import lzma  # noqa
except ImportError:
    HAS_XZ = False
else:
    HAS_XZ = True

# For when we need "some" test URLs
FEW = 5

# For stress testing the locking system using multiprocessing
N_PARALLEL_HAMMER = 10  # as high as 500 to replicate a bug

# For stress testing the locking system using threads
# (cheaper, works with coverage)
N_THREAD_HAMMER = 20  # as high as 1000 to replicate a bug


def url_to(path):
    return pathlib.Path(path).resolve().as_uri()


@pytest.fixture
def valid_urls(tmpdir):
    def _valid_urls(tmpdir):
        for i in itertools.count():
            c = os.urandom(16).hex()
            fn = os.path.join(tmpdir, "valid_" + str(i))
            with open(fn, "w") as f:
                f.write(c)
            u = url_to(fn)
            yield u, c

    return _valid_urls(tmpdir)


@pytest.fixture
def invalid_urls(tmpdir):
    def _invalid_urls(tmpdir):
        for i in itertools.count():
            fn = os.path.join(tmpdir, "invalid_" + str(i))
            if not os.path.exists(fn):
                yield url_to(fn)

    return _invalid_urls(tmpdir)


@pytest.fixture
def temp_cache(tmpdir):
    with paths.set_temp_cache(tmpdir):
        yield None
        check_download_cache(check_hashes=True)


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
def readonly_cache(tmpdir, valid_urls):
    with TemporaryDirectory(dir=tmpdir) as d:
        # other fixtures use the same tmpdir so we need a subdirectory
        # to make into the cache
        d = pathlib.Path(d)
        with paths.set_temp_cache(d):
            us = set(u for u, c in islice(valid_urls, FEW))
            urls = {u: download_file(u, cache=True) for u in us}
            files = set(d.iterdir())
            with readonly_dir(d):
                if not is_dir_readonly(d):
                    pytest.skip("Unable to make directory readonly")
                yield urls
            assert set(d.iterdir()) == files
            check_download_cache(check_hashes=True)


@pytest.fixture
def fake_readonly_cache(tmpdir, valid_urls, monkeypatch):
    def no_mkdir(p):
        raise PermissionError("os.mkdir monkeypatched out")

    with TemporaryDirectory(dir=tmpdir) as d:
        # other fixtures use the same tmpdir so we need a subdirectory
        # to make into the cache
        d = pathlib.Path(d)
        with paths.set_temp_cache(d):
            us = set(u for u, c in islice(valid_urls, FEW))
            urls = {u: download_file(u, cache=True) for u in us}
            files = set(d.iterdir())
            monkeypatch.setattr(os, "mkdir", no_mkdir)
            yield urls
            assert set(d.iterdir()) == files
            check_download_cache(check_hashes=True)


_shelve_possible_backends = ["dbm.dumb", "dbm.ndbm", "dbm.gnu"]


@contextlib.contextmanager
def shelve_backend(name):
    """Ensure that shelve has access only to one backend."""
    defaultmod = dbm._defaultmod
    modules = dbm._modules
    try:
        try:
            dbm._defaultmod = __import__(name, fromlist=['open'])
        except ImportError:
            pytest.skip(f"module {name} not available")
        dbm._modules = {name: dbm._defaultmod}
        yield
    finally:
        dbm._defaultmod = defaultmod
        dbm._modules = modules


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


@pytest.mark.parametrize("parallel", [False, True])
def test_download_with_sources_and_bogus_original(
        valid_urls, invalid_urls, temp_cache, parallel):

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
        sources[um] = []
        # For many of them the sources list starts with invalid URLs
        for iu in islice(invalid_urls, i):
            sources[um].append(iu)
        u, c = next(valid_urls)
        sources[um].append(u)
        urls.append((um, c, c_bad))

    # Now fetch them all
    if parallel:
        rs = download_files_in_parallel([u for (u, c, c_bad) in urls],
                                        cache=True,
                                        sources=sources)
    else:
        rs = [
            download_file(u, cache=True, sources=sources.get(u, None))
            for (u, c, c_bad) in urls
        ]
    assert len(rs) == len(urls)
    for r, (u, c, c_bad) in zip(rs, urls):
        assert get_file_contents(r) == c
        assert get_file_contents(r) != c_bad
        assert is_url_in_cache(u)


@pytest.mark.parametrize("b", _shelve_possible_backends)
def test_download_file_threaded_many(b, temp_cache, valid_urls):
    """Hammer download_file with multiple threaded requests.

    The goal is to stress-test the locking system. Normal parallel downloading
    also does this but coverage tools lose track of which paths are explored.

    """
    with shelve_backend(b):
        urls = list(islice(valid_urls, N_THREAD_HAMMER))
        with ThreadPoolExecutor(max_workers=len(urls)) as P:
            r = list(P.map(lambda u: download_file(u, cache=True),
                           [u for (u, c) in urls]))
        check_download_cache()
        assert len(r) == len(urls)
        for r, (u, c) in zip(r, urls):
            assert get_file_contents(r) == c


def test_download_file_threaded_many_partial_success(
        temp_cache, valid_urls, invalid_urls):
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
    for r, u in zip(r, urls):
        if u in contents:
            assert get_file_contents(r) == contents[u]
        else:
            assert r is None


def test_clear_download_multiple_references_doesnt_corrupt_storage(temp_cache, tmpdir):
    """Check that files with the same hash don't confuse the storage."""
    content = "Test data; doesn't matter much.\n"

    def make_url():
        with NamedTemporaryFile("w", dir=str(tmpdir), delete=False) as f:
            f.write(content)
        url = url_to(f.name)
        clear_download_cache(url)
        hash = download_file(url, cache=True)
        return url, hash

    a_url, a_hash = make_url()
    clear_download_cache(a_hash)
    assert not is_url_in_cache(a_url)

    f_url, f_hash = make_url()
    g_url, g_hash = make_url()

    assert f_url != g_url
    assert f_hash == g_hash
    assert is_url_in_cache(f_url)
    assert is_url_in_cache(g_url)

    clear_download_cache(f_url)
    assert not is_url_in_cache(f_url)
    assert is_url_in_cache(g_url)
    assert os.path.exists(
        g_hash
    ), "Contents should not be deleted while a reference exists"

    clear_download_cache(g_url)
    assert not os.path.exists(
        g_hash
    ), "No reference exists any more, file should be deleted"


def test_download_file_basic(valid_urls):
    primary, contents = next(valid_urls)
    f = download_file(primary, cache=False)
    assert get_file_contents(f) == contents


@pytest.mark.parametrize("use_cache", [False, True])
def test_download_file_local_cache_survives(tmpdir, temp_cache, use_cache):
    """Confirm that downloading a local file does not delete it.

    When implemented with urlretrieve (rather than urlopen) local files are
    not copied to create temporaries, so importing them to the cache deleted
    the original from wherever it was in the filesystem. I lost some built-in
    astropy data.

    """
    fn = tmpdir / "file"
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


def test_update_url(tmpdir, temp_cache):
    with TemporaryDirectory(dir=tmpdir) as d:
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
    assert get_file_contents(download_file(f_url, cache=True)) == "new", \
        "Cached version should still exist"
    with pytest.raises(urllib.error.URLError):
        # cannot download new version to check for updates
        download_file(f_url, cache="update")
    assert get_file_contents(download_file(f_url, cache=True)) == "new", \
        "Failed update should not remove the current version"


@pytest.mark.remote_data(source="astropy")
def test_download_noprogress():
    fnout = download_file(TESTURL, cache=False, show_progress=False)
    assert os.path.isfile(fnout)


@pytest.mark.remote_data(source="astropy")
def test_download_cache():

    download_dir = _get_download_cache_locs()[0]

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


def test_download_cache_after_clear(tmpdir, temp_cache, valid_urls):
    testurl, contents = next(valid_urls)
    # Test issues raised in #4427 with clear_download_cache() without a URL,
    # followed by subsequent download.
    download_dir = _get_download_cache_locs()[0]

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
        sources[urls[-1]] = [urls[-1], mirror_url+s]
    fnout = download_files_in_parallel(urls, sources=sources)
    assert all([os.path.isfile(f) for f in fnout]), fnout


@pytest.mark.parametrize("method", [None, "spawn"])
def test_download_parallel_fills_cache(tmpdir, valid_urls, method):
    urls = []
    # tmpdir is shared between many tests, and that can cause weird
    # interactions if we set the temporary cache too directly
    with paths.set_temp_cache(tmpdir):
        for um, c in islice(valid_urls, FEW):
            assert not is_url_in_cache(um)
            urls.append((um, c))
        rs = download_files_in_parallel(
            [u for (u, c) in urls], multiprocessing_start_method=method
        )
        assert len(rs) == len(urls)
        url_set = set(u for (u, c) in urls)
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
        sources[um] = []
        for iu in islice(invalid_urls, i):
            sources[um].append(iu)
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
    for r, (u, c) in zip(r, td):
        assert get_file_contents(r) == c


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


def test_download_parallel_partial_success_lock_safe(temp_cache, valid_urls, invalid_urls):
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


def test_download_parallel_update(temp_cache, tmpdir):
    td = []
    for i in range(N_PARALLEL_HAMMER):
        c = f"{i:04d}"
        fn = os.path.join(tmpdir, c)
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
    for (fn, u, c) in td:
        c_plus = c + " updated"
        fn = os.path.join(tmpdir, c)
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


@pytest.mark.remote_data(source="astropy")
def test_url_nocache():
    with get_readable_fileobj(TESTURL, cache=False, encoding="utf-8") as page:
        assert page.read().find("Astropy") > -1


@pytest.mark.remote_data(source="astropy")
def test_find_by_hash(valid_urls, temp_cache):
    testurl, contents = next(valid_urls)
    p = download_file(testurl, cache=True)
    hash = compute_hash(p)

    hashstr = "hash/" + hash

    fnout = get_pkg_data_filename(hashstr)
    assert os.path.isfile(fnout)
    clear_download_cache(hashstr[5:])
    assert not os.path.isfile(fnout)

    lockdir = os.path.join(_get_download_cache_locs()[0], "lock")
    assert not os.path.isdir(lockdir), "Cache dir lock was not released!"


@pytest.mark.remote_data(source="astropy")
def test_find_invalid():
    # this is of course not a real data file and not on any remote server, but
    # it should *try* to go to the remote server
    with pytest.raises(urllib.error.URLError):
        get_pkg_data_filename(
            "kjfrhgjklahgiulrhgiuraehgiurhgiuhreglhurieghruelighiuerahiulruli"
        )


# Package data functions
@pytest.mark.parametrize(
    ("filename"), ["local.dat", "local.dat.gz", "local.dat.bz2", "local.dat.xz"]
)
def test_local_data_obj(filename):
    if (not HAS_BZ2 and "bz2" in filename) or (not HAS_XZ and "xz" in filename):
        with pytest.raises(ValueError) as e:
            with get_pkg_data_fileobj(
                os.path.join("data", filename), encoding="binary"
            ) as f:
                f.readline()
                # assert f.read().rstrip() == b'CONTENT'
        assert " format files are not supported" in str(e.value)
    else:
        with get_pkg_data_fileobj(
            os.path.join("data", filename), encoding="binary"
        ) as f:
            f.readline()
            assert f.read().rstrip() == b"CONTENT"


@pytest.fixture(params=["invalid.dat.bz2", "invalid.dat.gz"])
def bad_compressed(request, tmpdir):
    # These contents have valid headers for their respective file formats, but
    # are otherwise malformed and invalid.
    bz_content = b"BZhinvalid"
    gz_content = b"\x1f\x8b\x08invalid"

    datafile = tmpdir.join(request.param)
    filename = datafile.strpath

    if filename.endswith(".bz2"):
        contents = bz_content
    elif filename.endswith(".gz"):
        contents = gz_content
    else:
        contents = "invalid"

    datafile.write(contents, mode="wb")

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
    if (not HAS_BZ2 and is_bz2) or (not HAS_XZ and is_xz):
        with pytest.raises(ValueError) as e:
            with get_readable_fileobj(bad_compressed, encoding="binary") as f:
                f.read()
        assert " format files are not supported" in str(e.value)
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
        assert filename == os.path.join(data_dir, "test_package", "data", "foo.txt")
    finally:
        sys.path.pop(0)


def test_local_data_nonlocalfail():
    # this would go *outside* the astropy tree
    with pytest.raises(RuntimeError):
        get_pkg_data_filename("../../../data/README.rst")


def test_compute_hash(tmpdir):

    rands = b"1234567890abcdefghijklmnopqrstuvwxyz"

    filename = tmpdir.join("tmp.dat").strpath

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

    # needed for testing the *real* lock at the end
    lockdir = os.path.join(_get_download_cache_locs('astropy')[0], 'lock')

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
        raise OSError
    monkeypatch.setattr(paths, '_find_or_create_root_dir', osraiser)

    with pytest.raises(OSError):
        # make sure the config dir search fails
        paths.get_cache_dir(rootname='astropy')

    # first try with cache
    with catch_warnings(CacheMissingWarning) as w:
        fnout = data.download_file(TESTURL, cache=True)

    assert os.path.isfile(fnout)

    has_inaccessible_warning = False
    has_temporary_warning = False

    for wi in w:
        if wi.category == CacheMissingWarning:
            if "Remote data cache could not be accessed" in wi.message.args[0]:
                has_inaccessible_warning = True
            if "temporary" in wi.message.args[0]:
                assert fnout == wi.message.args[1]
                has_temporary_warning = True
    assert has_inaccessible_warning
    assert has_temporary_warning

    # clearing the cache should be a no-up that doesn't affect fnout
    with catch_warnings(CacheMissingWarning) as w:
        data.clear_download_cache(TESTURL)
    assert os.path.isfile(fnout)

    # now remove it so tests don't clutter up the temp dir this should get
    # called at exit, anyway, but we do it here just to make sure it's working
    # correctly
    data._deltemps()
    assert not os.path.isfile(fnout)

    has_noclear_warning = False
    for wi in w:
        if wi.category == CacheMissingWarning:
            if "Not clearing data cache - cache inaccessible" in str(wi.message):
                has_noclear_warning = True
    assert has_noclear_warning

    # now try with no cache
    with catch_warnings(CacheMissingWarning) as w:
        fnnocache = data.download_file(TESTURL, cache=False)
    with open(fnnocache, "rb") as page:
        assert page.read().decode("utf-8").find("Astropy") > -1

    # no warnings should be raise in fileobj because cache is unnecessary
    assert len(w) == 0

    # lockdir determined above as the *real* lockdir, not the temp one
    assert not os.path.isdir(lockdir), "Cache dir lock was not released!"


@pytest.mark.parametrize(
    ("filename"),
    [
        "unicode.txt",
        "unicode.txt.gz",
        pytest.param(
            "unicode.txt.bz2",
            marks=pytest.mark.xfail(not HAS_BZ2, reason="no bz2 support"),
        ),
        pytest.param(
            "unicode.txt.xz",
            marks=pytest.mark.xfail(not HAS_XZ, reason="no lzma support"),
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
    assert x == (
        b"\xff\xd7\x94\xd7\x90\xd7\xa1\xd7\x98\xd7\xa8\xd7\x95\xd7\xa0"
        b"\xd7\x95\xd7\x9e\xd7\x99 \xd7\xa4\xd7\x99\xd7\x99\xd7\xaa\xd7\x95\xd7\x9f"[1:]
    )


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


def test_check_download_cache(tmpdir, temp_cache, valid_urls, invalid_urls):
    testurl, testurl_contents = next(valid_urls)
    testurl2, testurl2_contents = next(valid_urls)

    zip_file_name = os.path.join(tmpdir, "the.zip")
    clear_download_cache()
    check_download_cache()

    download_file(testurl, cache=True)
    # normal files probably corresponding to the urlmap
    normal = check_download_cache()
    download_file(testurl2, cache=True)
    assert check_download_cache() == normal

    export_download_cache(zip_file_name, [testurl, testurl2])
    assert check_download_cache(check_hashes=True) == normal

    clear_download_cache(testurl2)
    assert check_download_cache() == normal

    import_download_cache(zip_file_name, [testurl])
    assert check_download_cache(check_hashes=True) == normal


def test_export_import_roundtrip_one(tmpdir, temp_cache, valid_urls):
    testurl, contents = next(valid_urls)
    f = download_file(testurl, cache=True, show_progress=False)
    assert get_file_contents(f) == contents
    normal = check_download_cache()
    initial_urls_in_cache = set(get_cached_urls())
    zip_file_name = os.path.join(tmpdir, "the.zip")

    export_download_cache(zip_file_name, [testurl])
    clear_download_cache(testurl)
    import_download_cache(zip_file_name)
    assert is_url_in_cache(testurl)
    assert set(get_cached_urls()) == initial_urls_in_cache
    assert (
        get_file_contents(download_file(testurl, cache=True, show_progress=False))
        == contents
    )
    assert check_download_cache(check_hashes=True) == normal


def test_export_url_not_present(temp_cache, valid_urls):
    testurl, contents = next(valid_urls)
    with NamedTemporaryFile("wb") as zip_file:
        assert not is_url_in_cache(testurl)
        with pytest.raises(KeyError):
            export_download_cache(zip_file, [testurl])


def test_import_one(tmpdir, temp_cache, valid_urls):
    testurl, testurl_contents = next(valid_urls)
    testurl2, testurl2_contents = next(valid_urls)
    zip_file_name = os.path.join(tmpdir, "the.zip")

    download_file(testurl, cache=True)
    download_file(testurl2, cache=True)
    assert is_url_in_cache(testurl2)
    export_download_cache(zip_file_name, [testurl, testurl2])
    clear_download_cache(testurl)
    clear_download_cache(testurl2)
    import_download_cache(zip_file_name, [testurl])
    assert is_url_in_cache(testurl)
    assert not is_url_in_cache(testurl2)


def test_export_import_roundtrip(tmpdir, temp_cache, valid_urls):
    zip_file_name = os.path.join(tmpdir, "the.zip")
    for u, c in islice(valid_urls, FEW):
        download_file(u, cache=True)
    normal = check_download_cache()
    initial_urls_in_cache = set(get_cached_urls())

    export_download_cache(zip_file_name)
    clear_download_cache()
    import_download_cache(zip_file_name)

    assert set(get_cached_urls()) == initial_urls_in_cache
    assert check_download_cache(check_hashes=True) == normal


def test_export_import_roundtrip_stream(temp_cache, valid_urls):
    for u, c in islice(valid_urls, FEW):
        download_file(u, cache=True)
    normal = check_download_cache()
    initial_urls_in_cache = set(get_cached_urls())

    with io.BytesIO() as f:
        export_download_cache(f)
        b = f.getvalue()
    clear_download_cache()
    with io.BytesIO(b) as f:
        import_download_cache(f)

    assert set(get_cached_urls()) == initial_urls_in_cache
    assert check_download_cache(check_hashes=True) == normal


def test_export_overwrite_flag_works(temp_cache, valid_urls, tmpdir):
    fn = tmpdir / "f.zip"
    c = b"Some contents\nto check later"
    with open(fn, "wb") as f:
        f.write(c)
    for u, _ in islice(valid_urls, FEW):
        download_file(u, cache=True)

    with pytest.raises(FileExistsError):
        export_download_cache(fn)
    assert get_file_contents(fn, encoding='binary') == c

    export_download_cache(fn, overwrite=True)
    assert get_file_contents(fn, encoding='binary') != c


def test_export_import_roundtrip_different_location(tmpdir, valid_urls):
    original_cache = tmpdir / "original"
    os.mkdir(original_cache)
    zip_file_name = tmpdir / "the.zip"

    urls = list(islice(valid_urls, FEW))
    initial_urls_in_cache = set(u for (u, c) in urls)
    with paths.set_temp_cache(original_cache):
        for u, c in urls:
            download_file(u, cache=True)
        assert set(get_cached_urls()) == initial_urls_in_cache
        export_download_cache(zip_file_name)

    new_cache = tmpdir / "new"
    os.mkdir(new_cache)
    with paths.set_temp_cache(new_cache):
        import_download_cache(zip_file_name)
        check_download_cache(check_hashes=True)
        assert set(get_cached_urls()) == initial_urls_in_cache
        for (u, c) in urls:
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
    assert cache_total_size() == s_i + len(c)
    clear_download_cache(u)
    assert cache_total_size() == s_i


def test_cache_contents_agrees_with_get_urls(temp_cache, valid_urls):
    r = []
    for a, a_c in islice(valid_urls, FEW):
        a_f = download_file(a, cache=True)
        r.append((a, a_c, a_f))
    assert set(cache_contents().keys()) == set(get_cached_urls())
    for (u, c, h) in r:
        assert cache_contents()[u] == h


def test_free_space_checker_huge(tmpdir):
    with pytest.raises(OSError):
        check_free_space_in_dir(tmpdir, 1_000_000_000_000_000_000)


def test_get_free_space_file_directory(tmpdir):
    fn = tmpdir / "file"
    with open(fn, "w"):
        pass
    with pytest.raises(OSError):
        get_free_space_in_dir(fn)
    assert get_free_space_in_dir(tmpdir) > 0


def test_download_file_bogus_settings(invalid_urls, temp_cache):
    u = next(invalid_urls)
    with pytest.raises(KeyError):
        download_file(u, sources=[])


def test_download_file_local_directory(tmpdir):
    """Make sure we get a URLError rather than OSError even if it's a
    local directory."""
    with pytest.raises(urllib.request.URLError):
        download_file(url_to(tmpdir))


def test_download_file_schedules_deletion(valid_urls):
    u, c = next(valid_urls)
    f = download_file(u)
    assert f in _tempfilestodel
    # how to test deletion actually occurs?


def test_clear_download_cache_refuses_to_delete_outside_the_cache(tmpdir):
    fn = os.path.abspath(os.path.join(tmpdir, "file"))
    with open(fn, "w") as f:
        f.write("content")
    assert os.path.exists(fn)
    with pytest.raises(RuntimeError):
        clear_download_cache(fn)
    assert os.path.exists(fn)


def test_check_download_cache_finds_unreferenced_files(temp_cache, valid_urls):
    u, c = next(valid_urls)
    download_file(u, cache=True)
    with _cache(pkgname='astropy', write=True) as (dldir, urlmap):
        del urlmap[u]
    with pytest.raises(ValueError):
        check_download_cache()
    clear_download_cache()


def test_check_download_cache_finds_missing_files(temp_cache, valid_urls):
    u, c = next(valid_urls)
    os.remove(download_file(u, cache=True))
    with pytest.raises(ValueError):
        check_download_cache()
    clear_download_cache()


def test_check_download_cache_finds_bogus_entries(temp_cache, valid_urls):
    u, c = next(valid_urls)
    download_file(u, cache=True)
    with _cache(pkgname='astropy', write=True) as (dldir, urlmap):
        bd = os.path.join(dldir, "bogus")
        os.mkdir(bd)
        bf = os.path.join(bd, "file")
        with open(bf, "wt") as f:
            f.write("bogus file that exists")
        urlmap[u] = bf
    with pytest.raises(ValueError):
        check_download_cache()
    clear_download_cache()


def test_check_download_cache_finds_bogus_hashes(temp_cache, valid_urls):
    u, c = next(valid_urls)
    fn = download_file(u, cache=True)
    with open(fn, "w") as f:
        f.write("bogus contents")
    with pytest.raises(ValueError):
        check_download_cache(check_hashes=True)
    clear_download_cache()


def test_download_cache_update_doesnt_damage_cache(temp_cache, valid_urls):
    u, _ = next(valid_urls)
    download_file(u, cache=True)
    download_file(u, cache="update")


def test_cache_dir_is_actually_a_file(tmpdir, valid_urls):
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

    # set_temp_cache acts weird if it is pointed at a file (see below)
    # but we want to see what happens when the cache is pointed
    # at a file instead of a directory, so make a directory we can
    # replace later.
    fn = str(tmpdir / "file")
    ct = "contents\n"
    os.mkdir(fn)
    with paths.set_temp_cache(fn):
        shutil.rmtree(fn)
        with open(fn, "w") as f:
            f.write(ct)
        with pytest.raises(OSError):
            paths.get_cache_dir()
        check_quietly_ignores_bogus_cache()
    assert get_file_contents(fn) == ct, "File should not be harmed."

    # See what happens when set_temp_cache is pointed at a file
    with pytest.raises(OSError):
        with paths.set_temp_cache(fn):
            pass
    assert get_file_contents(str(fn)) == ct

    # Now the cache directory is normal but the subdirectory it wants
    # to make is a file
    cd = str(tmpdir / "astropy")
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmpdir):
        check_quietly_ignores_bogus_cache()
    assert get_file_contents(cd) == ct
    os.remove(cd)

    # Ditto one level deeper
    os.makedirs(cd)
    cd = str(tmpdir / "astropy" / "download")
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmpdir):
        check_quietly_ignores_bogus_cache()
    assert get_file_contents(cd) == ct
    os.remove(cd)

    # Ditto another level deeper
    os.makedirs(cd)
    py_version = "py" + str(sys.version_info.major)
    cd = str(tmpdir / "astropy" / "download" / py_version)
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmpdir):
        check_quietly_ignores_bogus_cache()
    assert get_file_contents(cd) == ct
    os.remove(cd)

    # Now interfere with creating the shelve object; this might actually
    # be okay if the shelve object has a funny naming convention
    # (the relation between the string you hand shelve and the names of
    # any files it may create is explicitly system-dependent)
    cd = str(tmpdir / "astropy" / "download" / py_version / "urlmap")
    os.makedirs(cd)
    with paths.set_temp_cache(tmpdir):
        check_quietly_ignores_bogus_cache()


def test_get_fileobj_str(a_file):
    fn, c = a_file
    with get_readable_fileobj(str(fn)) as rf:
        assert rf.read() == c


def test_get_fileobj_localpath(a_file):
    fn, c = a_file
    with get_readable_fileobj(py.path.local(fn)) as rf:
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
    with open(fn, "r") as f:
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


def test_cache_read_not_writable(temp_cache, valid_urls):
    with _cache(pkgname='astropy') as (dldir, urlmap):
        with pytest.raises(TypeError):
            urlmap["foo"] = 7
    u, _ = next(valid_urls)
    download_file(u, cache=True)
    with _cache(pkgname='astropy') as (dldir, urlmap):
        assert u in urlmap
        with pytest.raises(TypeError):
            urlmap["foo"] = 7


def test_cache_not_relocatable(tmpdir, valid_urls):
    u, c = next(valid_urls)
    d1 = tmpdir / "1"
    d2 = tmpdir / "2"
    os.mkdir(d1)
    with paths.set_temp_cache(d1):
        p1 = download_file(u, cache=True)
        assert is_url_in_cache(u)
        assert get_file_contents(p1) == c
        shutil.copytree(d1, d2)
        clear_download_cache()
    # this will not work! The filenames listed in the shelve are absolute
    # and so point back to the first cache
    with paths.set_temp_cache(d2):
        assert is_url_in_cache(u)
        p2 = download_file(u, cache=True)
        assert p1 == p2
        assert not os.path.exists(p2)
        with pytest.raises(RuntimeError):
            clear_download_cache(p2)
        with pytest.raises(CacheDamaged):
            check_download_cache()


def test_get_readable_fileobj_cleans_up_temporary_files(tmpdir, monkeypatch):
    """checks that get_readable_fileobj leaves no temporary files behind"""
    # Create a 'file://' URL pointing to a path on the local filesystem
    url = url_to(TESTLOCAL)

    # Save temporary files to a known location
    monkeypatch.setattr(tempfile, "tempdir", str(tmpdir))

    # Call get_readable_fileobj() as a context manager
    with get_readable_fileobj(url) as f:
        f.read()

    # Get listing of files in temporary directory
    tempdir_listing = tmpdir.listdir()

    # Assert that the temporary file was empty after get_readable_fileobj()
    # context manager finished running
    assert len(tempdir_listing) == 0


def test_path_objects_get_readable_fileobj():
    fpath = pathlib.Path(TESTLOCAL)
    with get_readable_fileobj(fpath) as f:
        assert f.read().rstrip() == (
            "This file is used in the test_local_data_* testing functions\nCONTENT"
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
    import contextlib
    from astropy.utils.data import download_file

    @contextlib.contextmanager
    def mockurl(remote_url, timeout=None):
        yield MockURL()

    class MockURL:
        def __init__(self):
            self.reader = io.BytesIO(b"a" * real_length)

        def info(self):
            return {"Content-Length": str(report_length)}

        def read(self, length=None):
            return self.reader.read(length)

    monkeypatch.setattr(urllib.request, "urlopen", mockurl)

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


def test_can_make_directories_readonly(tmpdir):
    try:
        with readonly_dir(tmpdir):
            assert is_dir_readonly(tmpdir)
    except AssertionError:
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            pytest.skip(
                "We are root, we can't make a directory un-writable with chmod."
            )
        elif platform.system() == "Windows":
            pytest.skip(
                "It seems we can't make a driectory un-writable under Windows "
                "with chmod, in spite of the documentation."
            )
        else:
            raise


def test_can_make_files_readonly(tmpdir):
    fn = tmpdir / "test"
    c = "contents\n"
    with open(fn, "w") as f:
        f.write(c)
    with readonly_dir(tmpdir):
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
        assert compute_hash(f) == os.path.basename(readonly_cache[u])


def test_check_download_cache_works_if_readonly(readonly_cache):
    check_download_cache(check_hashes=True)


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


def test_download_file_cache_fake_readonly_cache_miss(fake_readonly_cache, valid_urls):
    u, c = next(valid_urls)
    with pytest.warns(CacheMissingWarning):
        f = download_file(u, cache=True)
    assert get_file_contents(f) == c
    assert not is_url_in_cache(u)


def test_download_file_cache_fake_readonly_update(fake_readonly_cache):
    for u in fake_readonly_cache:
        with pytest.warns(CacheMissingWarning):
            f = download_file(u, cache="update")
        assert f != fake_readonly_cache[u]
        assert compute_hash(f) == os.path.basename(fake_readonly_cache[u])


def test_check_download_cache_works_if_fake_readonly(fake_readonly_cache):
    check_download_cache(check_hashes=True)


def test_pkgname_isolation(temp_cache, valid_urls):
    a = "bogus_cache_name"

    assert not get_cached_urls()
    assert not get_cached_urls(pkgname=a)

    for u, _ in islice(valid_urls, FEW):
        download_file(u, cache=True, pkgname=a)
    assert not get_cached_urls()
    assert len(get_cached_urls(pkgname=a)) == FEW
    assert cache_total_size() < cache_total_size(pkgname=a)

    for u, _ in islice(valid_urls, FEW+1):
        download_file(u, cache=True)
    assert len(get_cached_urls()) == FEW+1
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
    assert len(get_cached_urls()) == FEW+1, "wrong pkgname should do nothing"
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
    assert len(get_cached_urls()) == FEW+1
    assert len(get_cached_urls(pkgname=a)) == FEW-1

    clear_download_cache(u)
    assert len(get_cached_urls()) == FEW
    assert len(get_cached_urls(pkgname=a)) == FEW-1

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

    download_files_in_parallel([u for (u, c) in islice(valid_urls, FEW)],
                               pkgname=a)
    assert not get_cached_urls()
    assert len(get_cached_urls(pkgname=a)) == FEW


@pytest.mark.parametrize("b", _shelve_possible_backends)
def test_cache_with_different_shelve_backends(b, temp_cache, valid_urls):
    """Without special handling this emits a warning for dbm.dumb."""
    with shelve_backend(b):
        clear_download_cache()

        uc = list(islice(valid_urls, FEW))
        for u, c in uc:
            download_file(u, cache=True)

        check_download_cache()

        for u, c in uc:
            assert get_file_contents(
                download_file(u, cache=True, sources=[])) == c

        clear_download_cache()


@pytest.mark.parametrize("b", _shelve_possible_backends)
def test_lock_behaviour_if_directory_disappears(b, temp_cache):
    dldir, urlmapfn = _get_download_cache_locs()
    try:
        dbm = importlib.import_module(b)
    except ImportError:
        pytest.skip(f"module {b} not available")
    if b == "dbm.dumb":
        c = pytest.raises(FileNotFoundError)
    else:
        c = nullcontext()
    with dbm.open(urlmapfn, "c"):
        pass
    with c:
        with _cache("astropy", write=True) as (dldir, url2hash):
            url2hash["1"] = 2
            shutil.rmtree(dldir)


@pytest.mark.parametrize("b1b2", [(b1, b2)
                                  for b1 in _shelve_possible_backends
                                  for b2 in _shelve_possible_backends
                                  if b1 != b2
                                  ])
def test_wrong_backend_reports_useful_error(b1b2, temp_cache, valid_urls):
    b1, b2 = b1b2
    with shelve_backend(b1):
        for u, c in islice(valid_urls, FEW):
            download_file(u, cache=True)
        with shelve_backend(b2):
            for u, c in islice(valid_urls, FEW):
                with pytest.raises(dbm.error) as e:
                    download_file(u, cache=True)
                assert "module" in str(e.value)
                assert b1 in str(e.value)
            with pytest.raises(CacheDamaged):
                check_download_cache()
