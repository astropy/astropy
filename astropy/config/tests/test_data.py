# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ...tests.helper import remote_data,raises

import io

TESTURL = 'http://www.google.com/index.html'

@remote_data
def test_url_cache():
    from ..data import get_data_filename,clear_data_cache
    from os.path import isfile

    fnout = get_data_filename(TESTURL)
    assert isfile(fnout)
    clear_data_cache(TESTURL)
    assert not isfile(fnout)

@remote_data
def test_url_nocache():
    from ..data import get_data_fileobj

    with get_data_fileobj(TESTURL,cache=False) as googlepage:
        assert googlepage.read().decode().find('oogle</title>')>-1

@remote_data
def test_find_by_hash():
    from ..data import get_data_fileobj,get_data_filename,clear_data_cache
    from os.path import isfile
    import hashlib

    with get_data_fileobj(TESTURL) as googlepage:
        hash = hashlib.md5(googlepage.read())

    hashstr = 'hash/'+hash.hexdigest()

    fnout = get_data_filename(hashstr)
    assert isfile(fnout)
    clear_data_cache(hashstr[5:])
    assert not isfile(fnout)

def test_local_data_obj():
    from ..data import get_data_fileobj

    with get_data_fileobj('data/local.dat') as f:
        f.readline()
        assert f.read().rstrip() == b'CONTENT'

def test_local_data_name():
    from ..data import get_data_filename
    from os.path import isfile

    fnout = get_data_filename('data/local.dat')
    assert isfile(fnout) and fnout.endswith('local.dat')

    #get something in the astropy root
    fnout2 = get_data_filename('../../data/README.rst')
    assert isfile(fnout2) and fnout2.endswith('README.rst')

@raises(AssertionError)
def test_local_data_nonlocalfail():
    from ..data import get_data_filename

    #this would go *outside* the atropy tree
    fn = get_data_filename('../../../data/README.rst')

def test_compute_hash(tmpdir):
    import string
    import tempfile
    import hashlib
    from ..data import compute_hash

    rands = b'1234567890abcdefghijklmnopqrstuvwxyz'

    filename = tmpdir.join('tmp.dat').strpath

    with io.open(filename, 'wb') as ntf:
        ntf.write(rands)
        ntf.flush()

    chhash = compute_hash(filename)
    shash = hashlib.md5(rands).hexdigest()

    assert chhash==shash


def test_get_data_contents():
    from ..data import get_data_fileobj, get_data_contents

    with get_data_fileobj('data/local.dat') as f:
        contents1 = f.read()

    contents2 = get_data_contents('data/local.dat')

    assert contents1 == contents2

@remote_data
def test_data_noastropy_fallback(monkeypatch, recwarn):
    """
    Tests to make sure configuration items fall back to their defaults when
    there's a problem accessing the astropy directory
    """
    from os import path, remove
    from pytest import raises
    from .. import paths, data

    #better yet, set the configuration to make sure the temp files are deleted
    data.DELETE_TEMPORARY_DOWNLOADS_AT_EXIT.set(True)

    #make sure the config directory is not searched
    monkeypatch.setenv('XDG_CONFIG_HOME', 'foo')
    monkeypatch.delenv('XDG_CONFIG_HOME')

    # make sure the _find_or_create_astropy_dir function fails as though the
    # astropy dir could not be accessed
    def osraiser(dirnm, linkto):
        raise OSError
    monkeypatch.setattr(paths, '_find_or_create_astropy_dir', osraiser)


    with raises(OSError):
        #make sure the config dir search fails
        paths.get_cache_dir()

    #first try with cache
    fnout = data.get_data_filename(TESTURL)
    assert path.isfile(fnout)

    assert len(recwarn.list) > 1
    w1 = recwarn.pop()
    w2 = recwarn.pop()

    assert w1.category == data.CacheMissingWarning
    assert 'Remote data cache could not be accessed' in w1.message.args[0]
    assert w2.category == data.CacheMissingWarning
    assert 'File downloaded to temp file' in w2.message.args[0]
    assert fnout == w2.message.args[1]

    #clearing the cache should be a no-up that doesn't affect fnout
    data.clear_data_cache(TESTURL)
    assert path.isfile(fnout)

    #now remove it so tests don't clutter up the temp dir
    #this should get called at exit, anyway, but we do it here just to make
    #sure it's working correctly
    data._deltemps()
    assert not path.isfile(fnout)

    assert len(recwarn.list) > 0
    w3 = recwarn.pop()

    assert w3.category == data.CacheMissingWarning
    assert 'Not clearing data cache - cache inacessable' in str(w3.message)

    #now try with no cache
    with data.get_data_fileobj(TESTURL, cache=False) as googlepage:
        assert googlepage.read().decode().find('oogle</title>') > -1

    #no warnings should be raise in fileobj because cache is unnecessary
    assert len(recwarn.list) == 0
