# test helper.run_tests function

from .. import helper

# run_tests should raise ValueError when asked to run on a module it can't find
def test_module_not_found():
    with helper.pytest.raises(ValueError):
        helper.run_tests('fake.module')


# run_tests should raise ValueError when passed an invalid pastebin= option
def test_pastebin_keyword():
    with helper.pytest.raises(ValueError):
        helper.run_tests(pastebin='not_an_option')
