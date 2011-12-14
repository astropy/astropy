# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ...tests.helper import remote_data,raises

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
        assert f.read()==b'CONTENT\n'

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

def test_compute_hash():
    import string
    import random
    import tempfile
    import hashlib
    from ..data import compute_hash

    #generate a random string of 25 characters
    rands = ''.join(random.choice(string.ascii_letters) for x in range(25))

    with tempfile.NamedTemporaryFile('w+') as ntf:
        ntf.write(rands)
        ntf.flush()

        chhash = compute_hash(ntf.name)
        # the encode() call is necessary for py3.x compatibility
        shash = hashlib.md5(rands.encode()).hexdigest()

        assert chhash==shash


def test_get_data_file_contents():
    from ..data import get_data_fileobj, get_data_file_contents

    with get_data_fileobj('data/local.dat') as f:
        contents1 = f.read()

    contents2 = get_data_file_contents('data/local.dat')

    assert contents1 == contents2
