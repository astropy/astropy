# test helper.run_tests function

from astropy.tests import helper

def test_module_not_found():
    with helper.pytest.raises(ValueError):
        helper.run_tests('fake.module')
        
def test_pastebin_keyword():
    with helper.pytest.raises(ValueError):
        helper.run_tests(pastebin='not_an_option')
