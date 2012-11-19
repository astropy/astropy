"""
Test sky projections defined in WCS Paper II
"""
from __future__ import division
import os.path
import pytest
import unittest
import numpy as np
from numpy.testing import utils
import pywcs, pyfits
from .. import projections
from data import dpath

class TestProjections(unittest.TestCase):
    """
    Test composite models evaluation in series
    """
    def setUp(self):
        self.w = pywcs.WCS()
        hdr = pyfits.getheader(os.path.join(dpath,'1904-66_AZP.fits'))
        self.wazp = pywcs.WCS(hdr)
        self.wazp.wcs.crpix = np.array([0., 0.])
        self.wazp.wcs.crval = np.array([0., 0.])
        self.wazp.wcs.cdelt = np.array([1., 1.])
        
    def test_TAN_p2s(self):
        self.w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        wcslibout = self.w.wcs.p2s([[-10, 30]],1)
        wcs_phi = wcslibout['phi']
        wcs_theta = wcslibout['theta']
        tanprj = projections.Pix2Sky_TAN()
        phi, theta = tanprj(-10, 30)
        utils.assert_almost_equal(np.asarray(phi), wcs_phi)
        utils.assert_almost_equal(np.asarray(theta), wcs_theta)
        
    def test_TAN_s2p(self):
        self.w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
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
        azp = projections.Pix2Sky_AZP(mu=2, gamma=30)
        wcslibout = self.wazp.wcs.p2s([[-10, 30]],1)
        wcs_pix = self.wazp.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        azpinv = projections.Sky2Pix_AZP(mu=2, gamma=30)
        x, y = azpinv(wcslibout['phi'], wcslibout['theta'])
        utils.assert_almost_equal(np.asarray(x), wcs_pix[:,0])
        utils.assert_almost_equal(np.asarray(y), wcs_pix[:,1])
        
    def test_STG_p2s(self):
        self.w.wcs.ctype = ['RA---STG', 'DEC--STG']
        wcslibout = self.w.wcs.p2s([[-10, 30]],1)
        wcs_phi = wcslibout['phi']
        wcs_theta = wcslibout['theta']
        tanprj = projections.Pix2Sky_STG()
        phi, theta = tanprj(-10, 30)
        utils.assert_almost_equal(np.asarray(phi), wcs_phi)
        utils.assert_almost_equal(np.asarray(theta), wcs_theta)
        
    def test_STG_s2p(self):
        self.w.wcs.ctype = ['RA---STG', 'DEC--STG']
        wcslibout = self.w.wcs.p2s([[-10, 30]],1)
        wcs_pix = self.w.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        tinv = projections.Sky2Pix_STG()
        x, y = tinv(wcslibout['phi'], wcslibout['theta'])
        utils.assert_almost_equal(np.asarray(x), wcs_pix[:,0])
        utils.assert_almost_equal(np.asarray(y), wcs_pix[:,1])