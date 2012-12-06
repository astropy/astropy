# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ...tests.helper import remote_data, raises

import os
import io

import pytest
from ..data import _get_download_cache_locs

TESTURL = 'http://www.google.com/index.html'

# General file object function


@remote_data
def test_download_nocache():

    from ..data import download_file

    fnout = download_file(TESTURL)
    assert os.path.isfile(fnout)


@remote_data
def test_download_cache():

    from ..data import download_file, clear_download_cache

    fnout = download_file(TESTURL, cache=True)
    assert os.path.isfile(fnout)
    clear_download_cache(TESTURL)
    assert not os.path.isfile(fnout)

    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')
    assert not os.path.isdir(lockdir), 'Cache dir lock was not released!'


@remote_data
def test_url_nocache():

    from ..data import get_readable_fileobj

    with get_readable_fileobj(TESTURL, cache=False) as googlepage:
        assert googlepage.read().find('oogle</title>') > -1

@remote_data
def test_find_by_hash():

    from ..data import get_readable_fileobj, get_pkg_data_filename, clear_download_cache

    import hashlib

    with get_readable_fileobj(TESTURL, encoding="binary", cache=True) as googlepage:
        hash = hashlib.md5(googlepage.read())

    hashstr = 'hash/' + hash.hexdigest()

    fnout = get_pkg_data_filename(hashstr)
    assert os.path.isfile(fnout)
    clear_download_cache(hashstr[5:])
    assert not os.path.isfile(fnout)

    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')
    assert not os.path.isdir(lockdir), 'Cache dir lock was not released!'


# Package data functions
@pytest.mark.parametrize(('filename'), ['local.dat', 'local.dat.gz', 'local.dat.bz2'])
def test_local_data_obj(filename):
    from ..data import get_pkg_data_fileobj

    with get_pkg_data_fileobj(os.path.join('data', filename), encoding='binary') as f:
        f.readline()
        assert f.read().rstrip() == b'CONTENT'


@pytest.mark.parametrize(('filename'), ['invalid.dat.gz', 'invalid.dat.bz2'])
def test_local_data_obj_invalid(filename):
    from ..data import get_pkg_data_fileobj

    with get_pkg_data_fileobj(os.path.join('data', filename), encoding='binary') as f:
        assert f.read().rstrip().endswith(b'invalid')


def test_local_data_name():
    from ..data import get_pkg_data_filename

    fnout = get_pkg_data_filename('data/local.dat')
    assert os.path.isfile(fnout) and fnout.endswith('local.dat')

    #TODO: if in the future, the root data/ directory is added in, the below
    #test should be uncommented and the README.rst should be replaced with
    #whatever file is there

    #get something in the astropy root
    #fnout2 = get_pkg_data_filename('../../data/README.rst')
    #assert os.path.isfile(fnout2) and fnout2.endswith('README.rst')


@raises(AssertionError)
def test_local_data_nonlocalfail():
    from ..data import get_pkg_data_filename

    #this would go *outside* the atropy tree
    get_pkg_data_filename('../../../data/README.rst')


def test_compute_hash(tmpdir):
    import hashlib
    from ..data import compute_hash

    rands = b'1234567890abcdefghijklmnopqrstuvwxyz'

    filename = tmpdir.join('tmp.dat').strpath

    with io.open(filename, 'wb') as ntf:
        ntf.write(rands)
        ntf.flush()

    chhash = compute_hash(filename)
    shash = hashlib.md5(rands).hexdigest()

    assert chhash == shash


def test_get_pkg_data_contents():
    from ..data import get_pkg_data_fileobj, get_pkg_data_contents

    with get_pkg_data_fileobj('data/local.dat') as f:
        contents1 = f.read()

    contents2 = get_pkg_data_contents('data/local.dat')

    assert contents1 == contents2


@remote_data
def test_data_noastropy_fallback(monkeypatch, recwarn):
    """
    Tests to make sure the default behavior when the cache directory can't
    be located is correct
    """
    from pytest import raises
    from .. import data
    from ...config import paths

    # needed for testing the *real* lock at the end
    lockdir = os.path.join(_get_download_cache_locs()[0], 'lock')

    #better yet, set the configuration to make sure the temp files are deleted
    data.DELETE_TEMPORARY_DOWNLOADS_AT_EXIT.set(True)

    #make sure the config and cache directories are not searched
    monkeypatch.setenv('XDG_CONFIG_HOME', 'foo')
    monkeypatch.delenv('XDG_CONFIG_HOME')
    monkeypatch.setenv('XDG_CACHE_HOME', 'bar')
    monkeypatch.delenv('XDG_CACHE_HOME')

    # make sure the _find_or_create_astropy_dir function fails as though the
    # astropy dir could not be accessed
    def osraiser(dirnm, linkto):
        raise OSError
    monkeypatch.setattr(paths, '_find_or_create_astropy_dir', osraiser)

    with raises(OSError):
        #make sure the config dir search fails
        paths.get_cache_dir()

    #first try with cache
    fnout = data.download_file(TESTURL, cache=True)
    assert os.path.isfile(fnout)

    assert len(recwarn.list) > 1
    w1 = recwarn.pop()
    w2 = recwarn.pop()

    assert w1.category == data.CacheMissingWarning
    assert 'Remote data cache could not be accessed' in w1.message.args[0]
    assert w2.category == data.CacheMissingWarning
    assert 'File downloaded to temporary location' in w2.message.args[0]
    assert fnout == w2.message.args[1]

    #clearing the cache should be a no-up that doesn't affect fnout
    data.clear_download_cache(TESTURL)
    assert os.path.isfile(fnout)

    #now remove it so tests don't clutter up the temp dir
    #this should get called at exit, anyway, but we do it here just to make
    #sure it's working correctly
    data._deltemps()
    assert not os.path.isfile(fnout)

    assert len(recwarn.list) > 0
    w3 = recwarn.pop()

    assert w3.category == data.CacheMissingWarning
    assert 'Not clearing data cache - cache inacessable' in str(w3.message)

    #now try with no cache
    fnnocache = data.download_file(TESTURL, cache=False)
    with open(fnnocache, 'rb') as googlepage:
        assert googlepage.read().decode().find('oogle</title>') > -1

    #no warnings should be raise in fileobj because cache is unnecessary
    assert len(recwarn.list) == 0

    # lockdir determined above as the *real* lockdir, not the temp one
    assert not os.path.isdir(lockdir), 'Cache dir lock was not released!'


def test_read_unicode():
    from ..data import get_pkg_data_contents

    contents = get_pkg_data_contents('data/unicode.txt', encoding='utf-8')
    assert isinstance(contents, unicode)
    contents = contents.splitlines()[1]
    assert contents == u"\u05d4\u05d0\u05e1\u05d8\u05e8\u05d5\u05e0\u05d5\u05de\u05d9 \u05e4\u05d9\u05d9\u05ea\u05d5\u05df"

    contents = get_pkg_data_contents('data/unicode.txt', encoding='binary')
    assert isinstance(contents, bytes)
    contents = contents.splitlines()[1]
    assert contents == b"\xd7\x94\xd7\x90\xd7\xa1\xd7\x98\xd7\xa8\xd7\x95\xd7\xa0\xd7\x95\xd7\x9e\xd7\x99 \xd7\xa4\xd7\x99\xd7\x99\xd7\xaa\xd7\x95\xd7\x9f"


def test_compressed_stream():
    import base64
    from ..data import get_readable_fileobj

    gzipped_data = b"H4sICIxwG1AAA2xvY2FsLmRhdAALycgsVkjLzElVANKlxakpCpl5CiUZqQolqcUl8Tn5yYk58SmJJYnxWmCRzLx0hbTSvOSSzPy8Yi5nf78QV78QLgAlLytnRQAAAA=="
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
                self.data = b''
            else:
                result = self.data[:nbytes]
                self.data = self.data[nbytes:]
            return result

    stream = FakeStream(gzipped_data)
    with get_readable_fileobj(stream, encoding='binary') as f:
        f.readline()
        assert f.read().rstrip() == b'CONTENT'
