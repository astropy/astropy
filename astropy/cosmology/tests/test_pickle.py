# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pickle

import pytest

from astropy import cosmology
from astropy.cosmology import FLRW
from astropy.tests.helper import check_pickling_recovery, pickle_protocol

originals = [FLRW]
xfails = [False]


@pytest.mark.parametrize(("original", "xfail"), zip(originals, xfails))
def test_flrw(pickle_protocol, original, xfail):
    if xfail:
        pytest.xfail()
    check_pickling_recovery(original, pickle_protocol)


@pytest.mark.parametrize("name", cosmology.parameters.available)
@pytest.mark.parametrize("protocol", [0, 1, 2, 3, 4])  # add [5] when drop 3.7
def test_builtin_realizations(name, protocol):
    """
    Test in-built realizations can pickle and unpickle.
    Also a regression test for #12008.
    """
    # get class instance
    original = getattr(cosmology, name)

    # pickle and unpickle
    f = pickle.dumps(original, protocol=protocol)
    unpickled = pickle.loads(f)

    # test equality
    assert unpickled == original
    assert unpickled.meta == original.meta
