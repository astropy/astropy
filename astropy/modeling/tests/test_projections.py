# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test sky projections defined in WCS Paper II"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import os
import numpy as np
from numpy.testing import utils
from .. import projections
from ..parameters import InputParameterError

from ...io import fits
from ... import wcs
from ...utils.data import get_pkg_data_filename
from ...tests.helper import pytest


def test_Projection_properties():
    projection = projections.Sky2Pix_PlateCarree()
    assert projection.n_inputs == 2
    assert projection.n_outputs == 2


PIX_COORDINATES = [-10, 30]


pars = [
    ('TAN', projections.Sky2Pix_Gnomonic, {}),
    ('STG', projections.Sky2Pix_Stereographic, {}),
    ('SIN', projections.Sky2Pix_SlantOrthographic, {}),
    ('CEA', projections.Sky2Pix_CylindricalEqualArea, {}),
    ('CAR', projections.Sky2Pix_PlateCarree, {}),
    ('MER', projections.Sky2Pix_Mercator, {})
]


@pytest.mark.parametrize(('ID', 'model', 'args'), pars)
def test_Sky2Pix(ID, model, args):
    """Check astropy model eval against wcslib eval"""

    wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                           "1904-66_{0}.hdr".format(ID))
    test_file = get_pkg_data_filename(wcs_map)
    header = fits.Header.fromfile(test_file, endcard=False, padding=False)
    w = wcs.WCS(header)
    w.wcs.crval = [0., 0.]
    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = [1, 1]
    wcslibout = w.wcs.p2s([PIX_COORDINATES], 1)
    wcs_pix = w.wcs.s2p(wcslibout['world'], 1)['pixcrd']
    tinv = model(**args)
    x, y = tinv(wcslibout['phi'], wcslibout['theta'])
    utils.assert_almost_equal(np.asarray(x), wcs_pix[:, 0])
    utils.assert_almost_equal(np.asarray(y), wcs_pix[:, 1])

pars = [
    ('TAN', projections.Pix2Sky_Gnomonic, {}),
    ('STG', projections.Pix2Sky_Stereographic, {}),
    ('SIN', projections.Pix2Sky_SlantOrthographic, {}),
    ('CEA', projections.Pix2Sky_CylindricalEqualArea, {}),
    ('CAR', projections.Pix2Sky_PlateCarree, {}),
    ('MER', projections.Pix2Sky_Mercator, {}),
]


@pytest.mark.parametrize(('ID', 'model', 'args'), pars)
def test_Pix2Sky(ID, model, args):
    """Check astropy model eval against wcslib eval"""

    wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                           "1904-66_{0}.hdr".format(ID))
    test_file = get_pkg_data_filename(wcs_map)
    header = fits.Header.fromfile(test_file, endcard=False, padding=False)
    w = wcs.WCS(header)
    w.wcs.crval = [0., 0.]
    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = [1, 1]
    wcslibout = w.wcs.p2s([PIX_COORDINATES], 1)
    wcs_phi = wcslibout['phi']
    wcs_theta = wcslibout['theta']
    tanprj = model(**args)
    phi, theta = tanprj(*PIX_COORDINATES)
    utils.assert_almost_equal(np.asarray(phi), wcs_phi)
    utils.assert_almost_equal(np.asarray(theta), wcs_theta)


class TestZenithalPerspective(object):
    """Test Zenithal Perspective projection"""

    def setup_class(self):
        ID = 'AZP'
        wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                               "1904-66_{0}.hdr".format(ID))
        test_file = get_pkg_data_filename(wcs_map)
        header = fits.Header.fromfile(test_file, endcard=False, padding=False)
        self.wazp = wcs.WCS(header)
        self.wazp.wcs.crpix = np.array([0., 0.])
        self.wazp.wcs.crval = np.array([0., 0.])
        self.wazp.wcs.cdelt = np.array([1., 1.])
        self.pv_kw = [kw[2] for kw in self.wazp.wcs.get_pv()]
        self.azp = projections.Pix2Sky_ZenithalPerspective(*self.pv_kw)

    def test_AZP_p2s(self):
        wcslibout = self.wazp.wcs.p2s([[-10, 30]], 1)
        wcs_phi = wcslibout['phi']
        wcs_theta = wcslibout['theta']
        phi, theta = self.azp(-10, 30)
        utils.assert_almost_equal(np.asarray(phi), wcs_phi)
        utils.assert_almost_equal(np.asarray(theta), wcs_theta)

    def test_AZP_s2p(self):
        wcslibout = self.wazp.wcs.p2s([[-10, 30]], 1)
        wcs_pix = self.wazp.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        x, y = self.azp.inverse(wcslibout['phi'], wcslibout['theta'])
        utils.assert_almost_equal(np.asarray(x), wcs_pix[:, 0])
        utils.assert_almost_equal(np.asarray(y), wcs_pix[:, 1])


class TestCylindricalPerspective(object):
    """Test cylindrical perspective projection"""

    def setup_class(self):
        ID = "CYP"
        wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                               "1904-66_{0}.hdr".format(ID))
        test_file = get_pkg_data_filename(wcs_map)
        header = fits.Header.fromfile(test_file, endcard=False, padding=False)
        self.wazp = wcs.WCS(header)
        self.wazp.wcs.crpix = np.array([0., 0.])
        self.wazp.wcs.crval = np.array([0., 0.])
        self.wazp.wcs.cdelt = np.array([1., 1.])
        self.pv_kw = [kw[2] for kw in self.wazp.wcs.get_pv()]
        self.azp = projections.Pix2Sky_CylindricalPerspective(*self.pv_kw)

    def test_CYP_p2s(self):
        wcslibout = self.wazp.wcs.p2s([[-10, 30]], 1)
        wcs_phi = wcslibout['phi']
        wcs_theta = wcslibout['theta']
        phi, theta = self.azp(-10, 30)
        utils.assert_almost_equal(np.asarray(phi), wcs_phi)
        utils.assert_almost_equal(np.asarray(theta), wcs_theta)

    def test_CYP_s2p(self):
        wcslibout = self.wazp.wcs.p2s([[-10, 30]], 1)
        wcs_pix = self.wazp.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        x, y = self.azp.inverse(wcslibout['phi'], wcslibout['theta'])
        utils.assert_almost_equal(np.asarray(x), wcs_pix[:, 0])
        utils.assert_almost_equal(np.asarray(y), wcs_pix[:, 1])


def test_AffineTransformation2D():
    # Simple test with a scale and translation
    model = projections.AffineTransformation2D(
        matrix=[[2, 0], [0, 2]], translation=[1, 1])

    # Coordinates for vertices of a rectangle
    rect = [[0, 0], [1, 0], [0, 3], [1, 3]]

    x, y = zip(*rect)

    new_rect = np.vstack(model(x, y)).T

    assert np.all(new_rect == [[1, 1], [3, 1], [1, 7], [3, 7]])


def test_AffineTransformation2D_inverse():
    # Test non-invertible model
    model1 = projections.AffineTransformation2D(
        matrix=[[1, 1], [1, 1]])

    with pytest.raises(InputParameterError):
        model1.inverse

    model2 = projections.AffineTransformation2D(
        matrix=[[1.2, 3.4], [5.6, 7.8]], translation=[9.1, 10.11])

    # Coordinates for vertices of a rectangle
    rect = [[0, 0], [1, 0], [0, 3], [1, 3]]

    x, y = zip(*rect)

    x_new, y_new = model2.inverse(*model2(x, y))

    utils.assert_allclose([x, y], [x_new, y_new], atol=1e-10)


def test_c_projection_striding():
    # This is just a simple test to make sure that the striding is
    # handled correctly in the projection C extension
    coords = np.arange(10).reshape((5, 2))

    model = projections.Sky2Pix_ZenithalPerspective(2, 30)

    phi, theta = model(coords[:, 0], coords[:, 1])

    utils.assert_allclose(
        phi,
        [0., 2.31188456, 4.62283148, 6.9319002, 9.23814431])
    utils.assert_allclose(
        theta,
        [-76.48169176, -76.44546157, -76.33675621, -76.15553145, -75.90171376])
