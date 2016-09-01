# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from ... import constants as const
from ...tests.helper import pytest, pickle_protocol, check_pickling_recovery
from ...extern.six.moves import zip

originals = [const.Constant('h_fake', 'Not Planck',
                            0.0, 'J s', 0.0, 'fakeref',
                            system='si'),
             const.h,
             const.e]
xfails = [True, True, True]

@pytest.mark.parametrize(("original","xfail"), zip(originals, xfails))
def test_new_constant(pickle_protocol, original, xfail):
    if xfail:
        pytest.xfail()
    check_pickling_recovery(original, pickle_protocol)
