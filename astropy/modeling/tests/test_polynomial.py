# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for polynomial models."""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import numpy as np

from numpy.testing.utils import assert_allclose

from .. import fitting
from ...tests.helper import pytest
from ... import wcs
from ...io import fits
from ..polynomial import (Chebyshev1D, Legendre1D, Polynomial1D,
                          Chebyshev2D, Legendre2D, Polynomial2D, SIP,
                          OrthoPolynomialBase)
from ..functional_models import Linear1D
from ...utils.data import get_pkg_data_filename

try:
    from scipy import optimize  # pylint: disable=W0611
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


linear1d = {
    Chebyshev1D: {'parameters': [3],
                  'kwargs': {'c0': 1.2, 'c1': 2, 'c2': 2.3, 'c3': 0.2,
                             'domain': [1, 10]}},
    Legendre1D: {'parameters': [3],
                 'kwargs': {'c0': 1.2, 'c1': 2, 'c2': 2.3, 'c3': 0.2,
                            'domain': [1, 10]}},
    Polynomial1D: {'parameters': [3],
                   'kwargs': {'c0': 1.2, 'c1': 2, 'c2': 2.3, 'c3': 0.2}},
    Linear1D: {'parameters': [1.2, 23.1],
               'kwargs': {}}
}


linear2d = {
    Chebyshev2D: {'parameters': [1, 1],
                  'kwargs': {'c0_0': 1.2, 'c1_0': 2, 'c0_1': 2.3, 'c1_1': 0.2,
                             'x_domain': [0, 99], 'y_domain': [0, 82]}},
    Legendre2D: {'parameters': [1, 1],
                 'kwargs': {'c0_0': 1.2, 'c1_0': 2, 'c0_1': 2.3, 'c1_1': 0.2,
                            'x_domain': [0, 99], 'y_domain': [0, 82]}},
    Polynomial2D: {'parameters': [1],
                   'kwargs': {'c0_0': 1.2, 'c1_0': 2, 'c0_1': 2.3}},
}


@pytest.mark.skipif('not HAS_SCIPY')
class TestFitting(object):
    """Test linear fitter with polynomial models."""

    def setup_class(self):
        self.N = 100
        self.M = 100
        self.x1 = np.linspace(1, 10, 100)
        self.y2, self.x2 = np.mgrid[:100, :83]
        rsn = np.random.RandomState(0)
        self.n1 = rsn.randn(self.x1.size) * .1
        self.n2 = rsn.randn(self.x2.size)
        self.n2.shape = self.x2.shape
        self.linear_fitter = fitting.LinearLSQFitter()
        self.non_linear_fitter = fitting.LevMarLSQFitter()

    @pytest.mark.parametrize(('model_class'), linear1d.keys())
    def test_linear_fitter_1D(self, model_class):
        """Test fitting with LinearLSQFitter"""

        parameters = linear1d[model_class]['parameters']
        kwargs = linear1d[model_class]['kwargs']
        model = model_class(*parameters, **kwargs)
        y1 = model(self.x1)
        model_lin = self.linear_fitter(model, self.x1, y1 + self.n1)
        assert_allclose(model_lin.parameters, model.parameters,
                        atol=0.2)

    @pytest.mark.parametrize(('model_class'), linear1d.keys())
    def test_non_linear_fitter_1D(self, model_class):
        """Test fitting with non-linear LevMarLSQFitter"""

        parameters = linear1d[model_class]['parameters']
        kwargs = linear1d[model_class]['kwargs']
        model = model_class(*parameters, **kwargs)
        y1 = model(self.x1)
        model_nlin = self.non_linear_fitter(model, self.x1, y1 + self.n1)
        assert_allclose(model_nlin.parameters, model.parameters,
                        atol=0.2)

    @pytest.mark.parametrize(('model_class'), linear2d.keys())
    def test_linear_fitter_2D(self, model_class):
        """Test fitting with LinearLSQFitter"""

        parameters = linear2d[model_class]['parameters']
        kwargs = linear2d[model_class]['kwargs']
        model = model_class(*parameters, **kwargs)
        z = model(self.x2, self.y2)
        model_lin = self.linear_fitter(model, self.x2, self.y2, z + self.n2)
        assert_allclose(model_lin.parameters, model.parameters,
                        atol=0.2)

    @pytest.mark.parametrize(('model_class'), linear2d.keys())
    def test_non_linear_fitter_2D(self, model_class):
        """Test fitting with non-linear LevMarLSQFitter"""

        parameters = linear2d[model_class]['parameters']
        kwargs = linear2d[model_class]['kwargs']
        model = model_class(*parameters, **kwargs)
        z = model(self.x2, self.y2)
        model_nlin = self.non_linear_fitter(model, self.x2, self.y2,
                                            z + self.n2)
        assert_allclose(model_nlin.parameters, model.parameters,
                        atol=0.2)


def test_sip_hst():
    """Test SIP against astropy.wcs"""

    test_file = get_pkg_data_filename(os.path.join('data', 'hst_sip.hdr'))
    hdr = fits.Header.fromtextfile(test_file)
    crpix1 = hdr['CRPIX1']
    crpix2 = hdr['CRPIX2']
    wobj = wcs.WCS(hdr)
    a_pars = dict(**hdr['A_*'])
    b_pars = dict(**hdr['B_*'])
    a_order = a_pars.pop('A_ORDER')
    b_order = b_pars.pop('B_ORDER')
    sip = SIP([crpix1, crpix2], a_order, b_order, a_pars, b_pars)
    coords = [1, 1]
    rel_coords = [1 - crpix1, 1 - crpix2]
    astwcs_result = wobj.sip_pix2foc([coords], 1)[0] - rel_coords
    assert_allclose(sip(1, 1), astwcs_result)


def test_sip_irac():
    """Test forward and inverse SIP againts astropy.wcs"""

    test_file = get_pkg_data_filename(os.path.join('data', 'irac_sip.hdr'))
    hdr = fits.Header.fromtextfile(test_file)
    crpix1 = hdr['CRPIX1']
    crpix2 = hdr['CRPIX2']
    wobj = wcs.WCS(hdr)
    a_pars = dict(**hdr['A_*'])
    b_pars = dict(**hdr['B_*'])
    ap_pars = dict(**hdr['AP_*'])
    bp_pars = dict(**hdr['BP_*'])
    a_order = a_pars.pop('A_ORDER')
    b_order = b_pars.pop('B_ORDER')
    ap_order = ap_pars.pop('AP_ORDER')
    bp_order = bp_pars.pop('BP_ORDER')
    del a_pars['A_DMAX']
    del b_pars['B_DMAX']
    pix = [200, 200]
    rel_pix = [200 - crpix1, 200 - crpix2]
    sip = SIP([crpix1, crpix2], a_order, b_order, a_pars, b_pars,
              ap_order=ap_order, ap_coeff=ap_pars, bp_order=bp_order,
              bp_coeff=bp_pars)

    foc = wobj.sip_pix2foc([pix], 1)
    newpix = wobj.sip_foc2pix(foc, 1)[0]
    assert_allclose(sip(*pix), foc[0] - rel_pix)
    assert_allclose(sip.inverse(*foc[0]) +
                    foc[0] - rel_pix, newpix - pix)


def test_sip_no_coeff():
    sip = SIP([10,12], 2, 2)
    assert_allclose(sip.sip1d_a.parameters, [0., 0., 0])
    assert_allclose(sip.sip1d_b.parameters, [0., 0., 0])
    with pytest.raises(NotImplementedError):
        sip.inverse


@pytest.mark.parametrize('cls', (Polynomial1D, Chebyshev1D, Legendre1D,
                                 Polynomial2D, Chebyshev2D, Legendre2D))
def test_zero_degree_polynomial(cls):
    """
    A few tests that degree=0 polynomials are correctly evaluated and
    fitted.

    Regression test for https://github.com/astropy/astropy/pull/3589
    """

    if cls.n_inputs == 1:  # Test 1D polynomials
        p1 = cls(degree=0, c0=1)
        assert p1(0) == 1
        assert np.all(p1(np.zeros(5)) == np.ones(5))

        x = np.linspace(0, 1, 100)
        # Add a little noise along a straight line
        y = 1 + np.random.uniform(0, 0.1, len(x))

        p1_init = cls(degree=0)
        fitter = fitting.LinearLSQFitter()
        p1_fit = fitter(p1_init, x, y)

        # The fit won't be exact of course, but it should get close to within
        # 1%
        assert_allclose(p1_fit.c0, 1, atol=0.10)
    elif cls.n_inputs == 2:  # Test 2D polynomials
        if issubclass(cls, OrthoPolynomialBase):
            p2 = cls(x_degree=0, y_degree=0, c0_0=1)
        else:
            p2 = cls(degree=0, c0_0=1)
        assert p2(0, 0) == 1
        assert np.all(p2(np.zeros(5), np.zeros(5)) == np.ones(5))

        y, x = np.mgrid[0:1:100j,0:1:100j]
        z = (1 + np.random.uniform(0, 0.1, x.size)).reshape(100, 100)

        if issubclass(cls, OrthoPolynomialBase):
            p2_init = cls(x_degree=0, y_degree=0)
        else:
            p2_init = cls(degree=0)
        fitter = fitting.LinearLSQFitter()
        p2_fit = fitter(p2_init, x, y, z)

        assert_allclose(p2_fit.c0_0, 1, atol=0.10)
