TESTURL = 'http://www.google.com'

def test_url_cache():
    from ..data import get_data_filename,clear_data_cache
    from os.path import isfile
    
    fnout = get_data_filename(TESTURL)
    assert isfile(fnout)
    clear_data_cache(fnout)
    assert not isfile(fnout)
    
def test_url_cache_custom_fn():
    from ..data import get_data_filename,clear_data_cache
    from os.path import isfile
    
    fnout = get_data_filename(TESTURL,cachename='tester_custom_filename.html')
    assert fnout.endswith('tester_custom_filename.html')
    clear_data_cache(fnout)
    assert not isfile(fnout)    
    
def test_url_nocache():
    from ..data import get_data_fileobj
    
    googledata = get_data_fileobj(TESTURL,cache=False)