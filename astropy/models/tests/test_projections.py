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

def b(s):
    return s.encode('ascii')

class TestProjections(object):
    """
    Test composite models evaluation in series
    """
    def setup_class(self):
        self.w = wcs.WCS()
        test_file = get_pkg_data_filename(os.path.join('data', '1904-66_AZP.fits'))
        hdr = fits.getheader(test_file)
        self.wazp = wcs.WCS(hdr)
        self.wazp.wcs.crpix = np.array([0., 0.])
        self.wazp.wcs.crval = np.array([0., 0.])
        self.wazp.wcs.cdelt = np.array([1., 1.])
        
    def test_TAN_p2s(self):
        self.w.wcs.ctype = [b('RA---TAN'), b('DEC--TAN')]
        wcslibout = self.w.wcs.p2s([[-10, 30]],1)
        wcs_phi = wcslibout['phi']
        wcs_theta = wcslibout['theta']
        tanprj = projections.Pix2Sky_TAN()
        phi, theta = tanprj(-10, 30)
        utils.assert_almost_equal(np.asarray(phi), wcs_phi)
        utils.assert_almost_equal(np.asarray(theta), wcs_theta)
        
    def test_TAN_s2p(self):
        self.w.wcs.ctype = [b('RA---TAN'), b('DEC--TAN')]
        wcslibout = self.w.wcs.p2s([[-10, 30]],1)
        wcs_pix = self.w.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        tinv = projections.Sky2Pix_TAN()
        x, y = tinv(wcslibout['phi'], wcslibout['theta'])
        utils.assert_almost_equal(np.asarray(x), wcs_pix[:,0])
        utils.assert_almost_equal(np.asarray(y), wcs_pix[:,1])
        
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
        
    def test_STG_p2s(self):
        self.w.wcs.ctype = [b('RA---STG'), b('DEC--STG')]
        wcslibout = self.w.wcs.p2s([[-10, 30]],1)
        wcs_phi = wcslibout['phi']
        wcs_theta = wcslibout['theta']
        tanprj = projections.Pix2Sky_STG()
        phi, theta = tanprj(-10, 30)
        utils.assert_almost_equal(np.asarray(phi), wcs_phi)
        utils.assert_almost_equal(np.asarray(theta), wcs_theta)
        
    def test_STG_s2p(self):
        self.w.wcs.ctype = [b('RA---STG'), b('DEC--STG')]
        wcslibout = self.w.wcs.p2s([[-10, 30]],1)
        wcs_pix = self.w.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        tinv = projections.Sky2Pix_STG()
        x, y = tinv(wcslibout['phi'], wcslibout['theta'])
        utils.assert_almost_equal(np.asarray(x), wcs_pix[:,0])
        utils.assert_almost_equal(np.asarray(y), wcs_pix[:,1])
