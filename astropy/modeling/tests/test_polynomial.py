# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for polynomial models."""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import numpy as np

from numpy.testing import utils
from .. import fitting
from ...tests.helper import pytest
from ... import wcs
from ...io import fits
from ..polynomial import (Chebyshev1D, Legendre1D, Polynomial1D,
                          Chebyshev2D, Legendre2D, Polynomial2D, SIP)
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
        utils.assert_allclose(model_lin.parameters, model.parameters,
                              atol=0.2)

    @pytest.mark.parametrize(('model_class'), linear1d.keys())
    def test_non_linear_fitter_1D(self, model_class):
        """Test fitting with non-linear LevMarLSQFitter"""

        parameters = linear1d[model_class]['parameters']
        kwargs = linear1d[model_class]['kwargs']
        model = model_class(*parameters, **kwargs)
        y1 = model(self.x1)
        model_nlin = self.non_linear_fitter(model, self.x1, y1 + self.n1)
        utils.assert_allclose(model_nlin.parameters, model.parameters,
                              atol=0.2)

    @pytest.mark.parametrize(('model_class'), linear2d.keys())
    def test_linear_fitter_2D(self, model_class):
        """Test fitting with LinearLSQFitter"""

        parameters = linear2d[model_class]['parameters']
        kwargs = linear2d[model_class]['kwargs']
        model = model_class(*parameters, **kwargs)
        z = model(self.x2, self.y2)
        model_lin = self.linear_fitter(model, self.x2, self.y2, z + self.n2)
        utils.assert_allclose(model_lin.parameters, model.parameters,
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
        utils.assert_allclose(model_nlin.parameters, model.parameters,
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
    utils.assert_allclose(sip(1, 1), astwcs_result)


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
    utils.assert_allclose(sip(*pix), foc[0] - rel_pix)
    utils.assert_allclose(sip.inverse(*foc[0]) +
                          foc[0] - rel_pix, newpix - pix)


def test_sip_no_coeff():
    sip = SIP([10,12], 2, 2)
    utils.assert_allclose(sip.sip1d_a.parameters, [0., 0., 0])
    utils.assert_allclose(sip.sip1d_b.parameters, [0., 0., 0])
    with pytest.raises(NotImplementedError):
        sip.inverse
