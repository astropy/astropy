# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .. import utils
from ...tests.helper import pytest
from numpy.testing.utils import assert_allclose

def test_InputParameterError():
    message = 'no no no'
    err = utils.InputParameterError(message)
    assert str(err) == message

pars = [
        (1, (0, 10), (20, 40), 22),
        ([1.1, 2.2], (1.1, 2.2), (11, 22), [11, 22])
        ]

@pytest.mark.parametrize(('oldx', 'domain', 'window', 'newx'), pars)
def test_pmapdomain(oldx, domain, window, newx):
    assert_allclose(utils.pmapdomain(oldx, domain, window), newx)


pars = [
        (0, 0, 1),
        (3, 2, 3),
        (-1, 0, 0),
        (0, -1, 0),
        (2, 3, 0),
        ]

@pytest.mark.parametrize(('N', 'k', 'comb'), pars)
def test_comb(N, k, comb):
    assert utils.comb(N, k) == comb

