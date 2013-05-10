# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test sky projections defined in WCS Paper II
"""
from __future__ import division
import os.path
import numpy as np
from numpy.testing import utils
from .. import projections
from ...io import fits
from ... import wcs
from ...utils.data import get_pkg_data_filename
from ...tests.helper import pytest

def test_Projection_properties():
    projection = projections.Sky2Pix_CAR()
    assert projection.n_inputs == 2
    assert projection.n_outputs == 2
    assert projection.pdim == 1

PIX_COORDINATES = [-10, 30]
#PIX_COORDINATES = [-0.1, 0.3]

pars = [
        ('TAN', projections.Sky2Pix_TAN, {}),
        ('STG', projections.Sky2Pix_STG, {}),
        ('SIN', projections.Sky2Pix_SIN, {}),
#        ('CYP', projections.Sky2Pix_CYP, dict(mu=1, lam=1)),
#        ('CEA', projections.Sky2Pix_CEA, {}),
#        ('CAR', projections.Sky2Pix_CAR, {}),
#        ('MER', projections.Sky2Pix_MER, {}),
#        ('COP', projections.Sky2Pix_COP, {}),
        ]

@pytest.mark.parametrize(('ID', 'model', 'args'), pars)
def test_Sky2Pix(ID, model, args):
    """Check astropy model eval against wcslib eval"""
    w = wcs.WCS()
    w.wcs.ctype = ['RA---{0}'.format(ID),
                   'DEC--{0}'.format(ID)]
    wcslibout = w.wcs.p2s([PIX_COORDINATES],1)
    wcs_pix = w.wcs.s2p(wcslibout['world'], 1)['pixcrd']
    tinv = model(**args)
    x, y = tinv(wcslibout['phi'], wcslibout['theta'])
    utils.assert_almost_equal(np.asarray(x), wcs_pix[:,0])
    utils.assert_almost_equal(np.asarray(y), wcs_pix[:,1])

pars = [
        ('TAN', projections.Pix2Sky_TAN, {}),
        ('STG', projections.Pix2Sky_STG, {}),
        ('SIN', projections.Pix2Sky_SIN, {}),
#        ('CYP', projections.Pix2Sky_CYP, dict(mu=1, lam=1)),
#        ('CEA', projections.Pix2Sky_CEA, {}),
#        ('CAR', projections.Pix2Sky_CAR, {}),
#        ('MER', projections.Pix2Sky_MER, {}),
#        ('COP', projections.Pix2Sky_COP, {}),
        ]

@pytest.mark.parametrize(('ID', 'model', 'args'), pars)
def test_Pix2Sky(ID, model, args):
    """Check astropy model eval against wcslib eval"""
    w = wcs.WCS()
    w.wcs.ctype = ['RA---{0}'.format(ID),
                   'DEC--{0}'.format(ID)]
    wcslibout = w.wcs.p2s([PIX_COORDINATES],1)
    wcs_phi = wcslibout['phi']
    wcs_theta = wcslibout['theta']
    tanprj = model(**args)
    phi, theta = tanprj(*PIX_COORDINATES)
    utils.assert_almost_equal(np.asarray(phi), wcs_phi)
    utils.assert_almost_equal(np.asarray(theta), wcs_theta)

class TestAZP(object):
    """
    Test AZP projection
    """
    def setup_class(self):
        self.w = wcs.WCS()
        test_file = get_pkg_data_filename(os.path.join('data', '1904-66_AZP.fits'))
        hdr = fits.getheader(test_file)
        self.wazp = wcs.WCS(hdr)
        self.wazp.wcs.crpix = np.array([0., 0.])
        self.wazp.wcs.crval = np.array([0., 0.])
        self.wazp.wcs.cdelt = np.array([1., 1.])
        
    def test_AZP_p2s(self):
        azp = projections.Pix2Sky_AZP(mu=2, gamma=30)
        wcslibout = self.wazp.wcs.p2s([[-10, 30]],1)
        wcs_phi = wcslibout['phi']
        wcs_theta = wcslibout['theta']
        phi, theta = azp(-10, 30)
        utils.assert_almost_equal(np.asarray(phi), wcs_phi)
        utils.assert_almost_equal(np.asarray(theta), wcs_theta)

    def test_AZP_s2p(self):
        wcslibout = self.wazp.wcs.p2s([[-10, 30]],1)
        wcs_pix = self.wazp.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        azpinv = projections.Sky2Pix_AZP(mu=2, gamma=30)
        x, y = azpinv(wcslibout['phi'], wcslibout['theta'])
        utils.assert_almost_equal(np.asarray(x), wcs_pix[:,0])
        utils.assert_almost_equal(np.asarray(y), wcs_pix[:,1])
