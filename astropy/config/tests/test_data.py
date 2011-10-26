from ...tests.helper import remote_data

TESTURL = 'http://www.google.com/index.html'

@remote_data
def test_url_cache():
    from ..data import get_data_filename,clear_data_cache
    from os.path import isfile
    
    fnout = get_data_filename(TESTURL)
    assert isfile(fnout)
    clear_data_cache(fnout)
    assert not isfile(fnout)

@remote_data
def test_url_cache_custom_fn():
    from ..data import get_data_filename,clear_data_cache
    from os.path import isfile
    
    fnout = get_data_filename(TESTURL,cachename='tester_custom_filename.html')
    assert fnout[-27:] == 'tester_custom_filename.html'
    clear_data_cache(fnout)
    assert not isfile(fnout)    

@remote_data
def test_url_nocache():
    from ..data import get_data_fileobj
    
    googledata = get_data_fileobj(TESTURL,cache=False)
    assert googledata.read().find('oogle</title>')>-1
    
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
        shash = hashlib.md5(rands).hexdigest()
        
        assert chhash==shash