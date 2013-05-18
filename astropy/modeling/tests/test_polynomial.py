# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from numpy.testing.utils import assert_equal
from ...tests.helper import pytest
from .. import polynomial

# TODO: these tests fail at the moment.
# I'll file an issue with changes to make using these polynomials more uniform

pars = [
        (polynomial.Poly1DModel, 1),
        (polynomial.Poly2DModel, 2),
        (polynomial.Legendre1DModel, 1),
        (polynomial.Legendre2DModel, 2),
        (polynomial.Chebyshev1DModel, 1),
        (polynomial.Chebyshev2DModel, 2),
        ]

@pytest.mark.xfail
@pytest.mark.parametrize(('poly', 'ndim'), pars)
def test_constructor(poly, ndim):
    poly(1)
    poly(16)

    with pytest.raises(ValueError):
        poly(0)
        
    with pytest.raises(ValueError):
        poly(17)


@pytest.mark.xfail
@pytest.mark.parametrize(('poly', 'ndim'), pars)
def test_set_domain(poly, ndim):
    # TODO: Ask @nden if these results are OK.
    if ndim == 1:
        p = poly(3)
        xold = [1, 2]
        xnew = p.set_domain(xold)
        assert_equal(xnew, [1, 2])
    else: # must be ndim=2
        p = poly(3, 4)
        xold, yold = [1, 2], [3, 4]
        xnew, ynew = p.set_domain(xold, yold)
        assert_equal(xnew, [1, 2])
        assert_equal(ynew, [1, 2])

@pytest.mark.xfail
@pytest.mark.parametrize(('poly', 'ndim'), pars)
def test_eval(poly, ndim):
    # Try different kinds of input:
    for xx in [42, [42], [42, 43],
               np.array(42), np.array([42]), np.array([42, 43])]:
        if ndim == 1:
            p = poly(3)
            out = p(xx)
            # assert_equal(out, TODO)
        else: # must be ndim=2
            p = poly(3, 4)
            out = p(xx, xx)
            # assert_equal(out, TODO)

    