from ... import ascii as asciitable

from .common import (raises,
                     assert_equal, assert_almost_equal, assert_true,
                     setup_function, teardown_function, has_isnan)

def test_besancon():
    B = asciitable.read('t/besancon_test.txt',Reader=asciitable.besancon.BesanconFixed,guess=False)
    assert_equal(len(B),12)
