# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
import sys
import base64
import shutil
import hashlib
import pathlib
import tempfile
import warnings
import itertools
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from tempfile import NamedTemporaryFile, TemporaryDirectory

import py.path
import pytest

from astropy.utils import data
from astropy.config import paths
from astropy.utils.data import (
    CacheMissingWarning,
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
from astropy.tests.helper import raises, catch_warnings

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

n_parallel_hammer = 10
n_thread_hammer = 20


def url_to(path):
    return pathlib.Path(path).resolve().as_uri()


@pytest.fixture
def valid_urls(tmpdir):
    def _valid_urls(tmpdir):
        for i in itertools.count():
            c = os.urandom(16).hex()
            fn = os.path.join(tmpdir, "valid_"+str(i))
            with open(fn, "w") as f:
                f.write(c)
            u = url_to(fn)
            yield u, c
    return _valid_urls(tmpdir)


@pytest.fixture
def invalid_urls(tmpdir):
    def _invalid_urls(tmpdir):
        for i in itertools.count():
            fn = os.path.join(tmpdir, "invalid_"+str(i))
            assert not os.path.exists(fn)
            yield url_to(fn)
    return _invalid_urls(tmpdir)


@pytest.fixture
def temp_cache(tmpdir):
    with paths.set_temp_cache(tmpdir):
        yield None
        check_download_cache(check_hashes=True)


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


def test_download_with_sources_and_bogus_original(valid_urls,
                                                  invalid_urls,
                                                  temp_cache):
    u, c = next(valid_urls)
    urls = [(u, c, None)]
    sources = {}
    for i in range(5):
        um, c_bad = next(valid_urls)
        assert not is_url_in_cache(um)
        sources[um] = []
        for j in range(i):
            sources[um].append(next(invalid_urls))
        u, c = next(valid_urls)
        sources[um].append(u)
        urls.append((um, c, c_bad))
    rs = [download_file(u, cache=True, sources=sources.get(u, None))
          for (u, c, c_bad) in urls]
    assert len(rs) == len(urls)
    for r, (u, c, c_bad) in zip(rs, urls):
        assert get_file_contents(r) == c
        assert get_file_contents(r) != c_bad
        assert is_url_in_cache(u)


def test_download_file_threaded_many(temp_cache, valid_urls):
    """Hammer download_file with multiple threaded requests.

    The goal is to stress-test the locking system. Normal parallel downloading
    also does this but coverage tools lose track of which paths are explored.
    The hope is that if the parallelism happens in threads, and from here, the
    code coverage tools should be able to cope, maybe?

    """
    urls = [next(valid_urls) for i in range(n_thread_hammer)]
    with ThreadPoolExecutor(max_workers=len(urls)) as P:
        r = list(P.map(lambda u: download_file(u, cache=True),
                       [u for (u, c) in urls]))
    check_download_cache()
    assert len(r) == len(urls)
    for r, (u, c) in zip(r, urls):
        assert get_file_contents(r) == c


def test_clear_download_multiple_references_doesnt_corrupt_storage(
        temp_cache,
        tmpdir
):
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
    assert os.path.exists(g_hash), \
        "Contents should not be deleted while a reference exists"

    clear_download_cache(g_url)
    assert not os.path.exists(g_hash), \
        "No reference exists any more, file should be deleted"


def test_download_file_basic(valid_urls):
    primary, contents = next(valid_urls)
    f = download_file(primary, cache=False)
    assert get_file_contents(f) == contents


def test_download_file_local_survives(tmpdir):
    fn = tmpdir / "file"
    contents = "some text"
    with open(fn, "w") as f:
        f.write(contents)
    u = url_to(fn)
    f = download_file(u, cache=False)
    assert fn not in _tempfilestodel
    assert os.path.isfile(fn)
    assert get_file_contents(f) == contents


def test_download_file_local_cache_survives(tmpdir, temp_cache):
    fn = tmpdir / "file"
    contents = "some other text"
    with open(fn, "w") as f:
        f.write(contents)
    u = url_to(fn)
    f = download_file(u, cache=True)
    assert fn not in _tempfilestodel
    assert os.path.isfile(fn)
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
    f = download_file(primary,
                      cache=True,
                      sources=[primary, fallback1, fallback2])
    assert get_file_contents(f) == contents
    assert is_url_in_cache(primary)
    assert not is_url_in_cache(fallback1)
    assert not is_url_in_cache(fallback2)


def test_sources_multiple_missing(temp_cache, valid_urls, invalid_urls):
    primary = next(invalid_urls)
    fallback1 = next(invalid_urls)
    fallback2 = next(invalid_urls)
    with pytest.raises(urllib.error.URLError):
        download_file(primary,
                      cache=True,
                      sources=[primary, fallback1, fallback2])
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
        assert get_file_contents(download_file(f_url,
                                               cache=True,
                                               update_cache=True)) == "new"
    assert not os.path.exists(f_name)
    assert get_file_contents(download_file(f_url, cache=True)) == "new"
    with pytest.raises(urllib.error.URLError):
        download_file(f_url, cache=True, update_cache=True)
    assert get_file_contents(download_file(f_url, cache=True)) == "new"


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
def test_download_parallel_from_internet_works():
    main_url = conf.dataurl
    mirror_url = conf.dataurl_mirror
    fileloc = "intersphinx/README"
    try:
        fnout = download_files_in_parallel([main_url, main_url + fileloc])
    except urllib.error.URLError:  # Use mirror if timed out
        fnout = download_files_in_parallel([mirror_url, mirror_url + fileloc])
    assert all([os.path.isfile(f) for f in fnout]), fnout


@pytest.mark.parametrize("method", [None, "spawn"])
def test_download_parallel_fills_cache(tmpdir, valid_urls, method):
    urls = []
    # tmpdir is shared between many tests, and that can cause weird
    # interactions if we set the temporary cache too directly
    with paths.set_temp_cache(tmpdir):
        for i in range(5):
            um, c = next(valid_urls)
            assert not is_url_in_cache(um)
            urls.append((um, c))
        rs = download_files_in_parallel([u for (u, c) in urls],
                                        multiprocessing_start_method=method)
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
    for i in range(5):
        um, c = next(valid_urls)
        assert not is_url_in_cache(um)
        urls.append((um, c))
    rs = download_files_in_parallel([u for (u, c) in urls], sources=sources)
    assert len(rs) == len(urls)
    # u = set(u for (u, c) in urls)
    # assert u <= set(get_cached_urls())
    check_download_cache()
    for r, (u, c) in zip(rs, urls):
        assert get_file_contents(r) == c


def test_download_parallel_with_sources_and_bogus_original(valid_urls,
                                                           invalid_urls,
                                                           temp_cache):
    u, c = next(valid_urls)
    urls = [(u, c, None)]
    sources = {}
    for i in range(5):
        um, c_bad = next(valid_urls)
        assert not is_url_in_cache(um)
        sources[um] = []
        for j in range(i):
            sources[um].append(next(invalid_urls))
        u, c = next(valid_urls)
        sources[um].append(u)
        urls.append((um, c, c_bad))
    rs = download_files_in_parallel([u for (u, c, c_bad) in urls],
                                    sources=sources)
    assert len(rs) == len(urls)
    # u = set(u for (u, c, c_bad) in urls)
    # assert u <= set(get_cached_urls())
    for r, (u, c, c_bad) in zip(rs, urls):
        assert get_file_contents(r) == c
        assert get_file_contents(r) != c_bad


def test_download_parallel_many(temp_cache, valid_urls):
    td = []
    for i in range(n_parallel_hammer):
        u, c = next(valid_urls)
        clear_download_cache(u)
        td.append((u, c))

    r = download_files_in_parallel([u for (u, c) in td])
    assert len(r) == len(td)
    for r, (u, c) in zip(r, td):
        assert get_file_contents(r) == c


def test_download_parallel_partial_success(temp_cache,
                                           valid_urls,
                                           invalid_urls):
    td = []
    for i in range(n_parallel_hammer):
        u, c = next(valid_urls)
        clear_download_cache(u)
        td.append((u, c))

    u_bad = next(invalid_urls)

    with pytest.raises(urllib.request.URLError):
        download_files_in_parallel([u_bad] + [u for (u, c) in td])
    # Actually some files may get downloaded, others not.
    # Is this good? Should we stubbornly keep trying?
    # assert not any([is_url_in_cache(u) for (u, c) in td])


def test_download_parallel_update(temp_cache, tmpdir):
    td = []
    for i in range(n_parallel_hammer):
        c = "%04d" % i
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

    r2 = download_files_in_parallel([u for (fn, u, c) in td],
                                    update_cache=False)
    assert len(r2) == len(td)
    for r_2, (fn, u, c, c_plus) in zip(r2, td2):
        assert get_file_contents(r_2) == c
        assert c != c_plus
    r3 = download_files_in_parallel([u for (fn, u, c) in td],
                                    update_cache=True)

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
            with get_pkg_data_fileobj(os.path.join("data", filename), encoding="binary") as f:
                f.readline()
                # assert f.read().rstrip() == b'CONTENT'
        assert " format files are not supported" in str(e.value)
    else:
        with get_pkg_data_fileobj(os.path.join("data", filename), encoding="binary") as f:
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


@raises(RuntimeError)
def test_local_data_nonlocalfail():
    # this would go *outside* the astropy tree
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
        warnings.simplefilter("ignore", ResourceWarning)
        fnout = data.download_file(TESTURL, cache=True)

    assert os.path.isfile(fnout)

    assert len(w) > 1

    w1 = w.pop(0)
    w2 = w.pop(0)

    assert w1.category == CacheMissingWarning
    assert "Remote data cache could not be accessed" in w1.message.args[0]
    assert w2.category == CacheMissingWarning
    assert "File downloaded to temporary location" in w2.message.args[0]
    assert fnout == w2.message.args[1]

    # clearing the cache should be a no-up that doesn't affect fnout
    with catch_warnings(CacheMissingWarning) as w:
        data.clear_download_cache(TESTURL)
    assert os.path.isfile(fnout)

    # now remove it so tests don't clutter up the temp dir this should get
    # called at exit, anyway, but we do it here just to make sure it's working
    # correctly
    data._deltemps()
    assert not os.path.isfile(fnout)

    assert len(w) > 0
    w3 = w.pop()

    assert w3.category == data.CacheMissingWarning
    assert "Not clearing data cache - cache inacessible" in str(w3.message)

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
    assert get_file_contents(download_file(testurl,
                                           cache=True,
                                           show_progress=False)) == contents
    assert check_download_cache(check_hashes=True) == normal


def test_export_with_download(temp_cache, valid_urls):
    testurl, contents = next(valid_urls)
    with NamedTemporaryFile("wb") as zip_file:
        clear_download_cache(testurl)
        export_download_cache(zip_file, [testurl])
        assert is_url_in_cache(testurl)


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
    for u, c in itertools.islice(valid_urls, 7):
        download_file(u, cache=True)
    normal = check_download_cache()
    initial_urls_in_cache = set(get_cached_urls())

    export_download_cache(zip_file_name)
    clear_download_cache()
    import_download_cache(zip_file_name)

    assert set(get_cached_urls()) == initial_urls_in_cache
    assert check_download_cache(check_hashes=True) == normal


def test_export_import_roundtrip_different_location(tmpdir, valid_urls):
    original_cache = tmpdir / "original"
    os.mkdir(original_cache)
    zip_file_name = tmpdir / "the.zip"

    urls = list(itertools.islice(valid_urls, 7))
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
        temp_cache,
        valid_urls):
    u, c = next(valid_urls)
    clear_download_cache(u)
    s_i = cache_total_size()
    download_file(u, cache=True)
    assert cache_total_size() == s_i + len(c)
    clear_download_cache(u)
    assert cache_total_size() == s_i


def test_cache_contents_agrees_with_get_urls(temp_cache, valid_urls):
    r = []
    for i in range(7):
        a, a_c = next(valid_urls)
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
    with pytest.raises(ValueError):
        download_file(u, cache=False, update_cache=True)
    with pytest.raises(ValueError):
        download_file(u, sources=[])


def test_download_file_local_directory(tmpdir):
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
    with _cache(write=True) as (dldir, urlmap):
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
    with _cache(write=True) as (dldir, urlmap):
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
    download_file(u, cache=True, update_cache=True)


def test_cache_dir_is_actually_a_file(tmpdir, valid_urls):
    def check_quietly_ignores_bogus_cache():
        with catch_warnings(CacheMissingWarning) as w:
            assert not get_cached_urls()
        assert any(wi.category == CacheMissingWarning for wi in w)
        with catch_warnings(CacheMissingWarning) as w:
            assert not is_url_in_cache("http://www.example.com/")
        assert any(wi.category == CacheMissingWarning for wi in w)
        with catch_warnings(CacheMissingWarning) as w:
            assert not cache_contents()
        assert any(wi.category == CacheMissingWarning for wi in w)
        with catch_warnings(CacheMissingWarning) as w:
            u, c = next(valid_urls)
            r = download_file(u, cache=True)
            assert get_file_contents(r) == c
        assert any(wi.category == CacheMissingWarning for wi in w)
        with catch_warnings(CacheMissingWarning) as w:
            assert not is_url_in_cache(u)
        assert any(wi.category == CacheMissingWarning for wi in w)
        with pytest.raises(OSError):
            check_download_cache()

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
    assert get_file_contents(fn) == ct

    with pytest.raises(OSError):
        with paths.set_temp_cache(fn):
            pass
    assert get_file_contents(str(fn)) == ct

    cd = str(tmpdir / "astropy")
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmpdir):
        check_quietly_ignores_bogus_cache()
    assert get_file_contents(cd) == ct
    os.remove(cd)

    os.makedirs(cd)
    cd = str(tmpdir / "astropy" / "download")
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmpdir):
        check_quietly_ignores_bogus_cache()
    assert get_file_contents(cd) == ct
    os.remove(cd)

    os.makedirs(cd)
    py_version = 'py' + str(sys.version_info.major)
    cd = str(tmpdir / "astropy" / "download" / py_version)
    with open(cd, "w") as f:
        f.write(ct)
    with paths.set_temp_cache(tmpdir):
        check_quietly_ignores_bogus_cache()
    assert get_file_contents(cd) == ct
    os.remove(cd)

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
    with _cache() as (dldir, urlmap):
        with pytest.raises(TypeError):
            urlmap["foo"] = 7
    u, _ = next(valid_urls)
    download_file(u, cache=True)
    with _cache() as (dldir, urlmap):
        assert u in urlmap
        with pytest.raises(TypeError):
            urlmap["foo"] = 7


def test_cache_not_relocatable(tmpdir, valid_urls):
    u, c = next(valid_urls)
    d1 = tmpdir / "1"
    os.mkdir(d1)
    with paths.set_temp_cache(d1):
        p1 = download_file(u, cache=True)
        assert is_url_in_cache(u)
        assert get_file_contents(p1) == c
    d2 = tmpdir / "2"
    # this will not work!
    shutil.copytree(d1, d2)
    with paths.set_temp_cache(d2):
        assert is_url_in_cache(u)
        p2 = download_file(u, cache=True)
        assert p1 == p2


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
            "This file is used in the test_local_data_* "
            "testing functions\nCONTENT"
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
            self.reader = io.BytesIO(b"a"*real_length)

        def info(self):
            return {"Content-Length": str(report_length)}

        def read(self, length=None):
            return self.reader.read(length)

    monkeypatch.setattr(urllib.request, "urlopen", mockurl)

    with pytest.raises(urllib.error.URLError):
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
        assert f.read() == b"a"*real_length

    report_length = None
    real_length = 1023
    fn = download_file(TESTURL, cache=False)
    with open(fn, "rb") as f:
        assert f.read() == b"a"*real_length
