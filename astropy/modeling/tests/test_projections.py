# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test sky projections defined in WCS Paper II"""

import os

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal

from .. import projections
from ..parameters import InputParameterError

from ... import units as u
from ...io import fits
from ... import wcs
from ...utils.data import get_pkg_data_filename
from ...tests.helper import assert_quantity_allclose


def test_Projection_properties():
    projection = projections.Sky2Pix_PlateCarree()
    assert projection.n_inputs == 2
    assert projection.n_outputs == 2


PIX_COORDINATES = [-10, 30]


pars = [(x,) for x in projections.projcodes]
# There is no groundtruth file for the XPH projection available here:
#   http://www.atnf.csiro.au/people/mcalabre/WCS/example_data.html
pars.remove(('XPH',))


@pytest.mark.parametrize(('code',), pars)
def test_Sky2Pix(code):
    """Check astropy model eval against wcslib eval"""

    wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                           "1904-66_{0}.hdr".format(code))
    test_file = get_pkg_data_filename(wcs_map)
    header = fits.Header.fromfile(test_file, endcard=False, padding=False)

    params = []
    for i in range(3):
        key = 'PV2_{0}'.format(i + 1)
        if key in header:
            params.append(header[key])

    w = wcs.WCS(header)
    w.wcs.crval = [0., 0.]
    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = [1, 1]
    wcslibout = w.wcs.p2s([PIX_COORDINATES], 1)
    wcs_pix = w.wcs.s2p(wcslibout['world'], 1)['pixcrd']
    model = getattr(projections, 'Sky2Pix_' + code)
    tinv = model(*params)
    x, y = tinv(wcslibout['phi'], wcslibout['theta'])
    assert_almost_equal(np.asarray(x), wcs_pix[:, 0])
    assert_almost_equal(np.asarray(y), wcs_pix[:, 1])


@pytest.mark.parametrize(('code',), pars)
def test_Pix2Sky(code):
    """Check astropy model eval against wcslib eval"""

    wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                           "1904-66_{0}.hdr".format(code))
    test_file = get_pkg_data_filename(wcs_map)
    header = fits.Header.fromfile(test_file, endcard=False, padding=False)

    params = []
    for i in range(3):
        key = 'PV2_{0}'.format(i + 1)
        if key in header:
            params.append(header[key])

    w = wcs.WCS(header)
    w.wcs.crval = [0., 0.]
    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = [1, 1]
    wcslibout = w.wcs.p2s([PIX_COORDINATES], 1)
    wcs_phi = wcslibout['phi']
    wcs_theta = wcslibout['theta']
    model = getattr(projections, 'Pix2Sky_' + code)
    tanprj = model(*params)
    phi, theta = tanprj(*PIX_COORDINATES)
    assert_almost_equal(np.asarray(phi), wcs_phi)
    assert_almost_equal(np.asarray(theta), wcs_theta)


@pytest.mark.parametrize(('code',), pars)
def test_Sky2Pix_unit(code):
    """Check astropy model eval against wcslib eval"""

    wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                           "1904-66_{0}.hdr".format(code))
    test_file = get_pkg_data_filename(wcs_map)
    header = fits.Header.fromfile(test_file, endcard=False, padding=False)

    params = []
    for i in range(3):
        key = 'PV2_{0}'.format(i + 1)
        if key in header:
            params.append(header[key])

    w = wcs.WCS(header)
    w.wcs.crval = [0., 0.]
    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = [1, 1]
    wcslibout = w.wcs.p2s([PIX_COORDINATES], 1)
    wcs_pix = w.wcs.s2p(wcslibout['world'], 1)['pixcrd']
    model = getattr(projections, 'Sky2Pix_' + code)
    tinv = model(*params)
    x, y = tinv(wcslibout['phi'] * u.deg, wcslibout['theta'] * u.deg)
    assert_quantity_allclose(x, wcs_pix[:, 0] * u.deg)
    assert_quantity_allclose(y, wcs_pix[:, 1] * u.deg)


@pytest.mark.parametrize(('code',), pars)
def test_Pix2Sky_unit(code):
    """Check astropy model eval against wcslib eval"""

    wcs_map = os.path.join(os.pardir, os.pardir, "wcs", "tests", "maps",
                           "1904-66_{0}.hdr".format(code))
    test_file = get_pkg_data_filename(wcs_map)
    header = fits.Header.fromfile(test_file, endcard=False, padding=False)

    params = []
    for i in range(3):
        key = 'PV2_{0}'.format(i + 1)
        if key in header:
            params.append(header[key])

    w = wcs.WCS(header)
    w.wcs.crval = [0., 0.]
    w.wcs.crpix = [0, 0]
    w.wcs.cdelt = [1, 1]
    wcslibout = w.wcs.p2s([PIX_COORDINATES], 1)
    wcs_phi = wcslibout['phi']
    wcs_theta = wcslibout['theta']
    model = getattr(projections, 'Pix2Sky_' + code)
    tanprj = model(*params)
    phi, theta = tanprj(*PIX_COORDINATES * u.deg)
    assert_quantity_allclose(phi, wcs_phi * u.deg)
    assert_quantity_allclose(theta, wcs_theta * u.deg)
    phi, theta = tanprj(*(PIX_COORDINATES * u.deg).to(u.rad))
    assert_quantity_allclose(phi, wcs_phi * u.deg)
    assert_quantity_allclose(theta, wcs_theta * u.deg)
    phi, theta = tanprj(*(PIX_COORDINATES * u.deg).to(u.arcmin))
    assert_quantity_allclose(phi, wcs_phi * u.deg)
    assert_quantity_allclose(theta, wcs_theta * u.deg)


@pytest.mark.parametrize(('code',), pars)
def test_projection_default(code):
    """Check astropy model eval with default parameters"""
    # Just makes sure that the default parameter values are reasonable
    # and accepted by wcslib.

    model = getattr(projections, 'Sky2Pix_' + code)
    tinv = model()
    x, y = tinv(45, 45)

    model = getattr(projections, 'Pix2Sky_' + code)
    tinv = model()
    x, y = tinv(0, 0)


class TestZenithalPerspective:
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
        assert_almost_equal(np.asarray(phi), wcs_phi)
        assert_almost_equal(np.asarray(theta), wcs_theta)

    def test_AZP_s2p(self):
        wcslibout = self.wazp.wcs.p2s([[-10, 30]], 1)
        wcs_pix = self.wazp.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        x, y = self.azp.inverse(wcslibout['phi'], wcslibout['theta'])
        assert_almost_equal(np.asarray(x), wcs_pix[:, 0])
        assert_almost_equal(np.asarray(y), wcs_pix[:, 1])


class TestCylindricalPerspective:
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
        assert_almost_equal(np.asarray(phi), wcs_phi)
        assert_almost_equal(np.asarray(theta), wcs_theta)

    def test_CYP_s2p(self):
        wcslibout = self.wazp.wcs.p2s([[-10, 30]], 1)
        wcs_pix = self.wazp.wcs.s2p(wcslibout['world'], 1)['pixcrd']
        x, y = self.azp.inverse(wcslibout['phi'], wcslibout['theta'])
        assert_almost_equal(np.asarray(x), wcs_pix[:, 0])
        assert_almost_equal(np.asarray(y), wcs_pix[:, 1])


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

    assert_allclose([x, y], [x_new, y_new], atol=1e-10)


def test_c_projection_striding():
    # This is just a simple test to make sure that the striding is
    # handled correctly in the projection C extension
    coords = np.arange(10).reshape((5, 2))

    model = projections.Sky2Pix_ZenithalPerspective(2, 30)

    phi, theta = model(coords[:, 0], coords[:, 1])

    assert_almost_equal(
        phi,
        [0., 2.2790416, 4.4889294, 6.6250643, 8.68301])

    assert_almost_equal(
        theta,
        [-76.4816918, -75.3594654, -74.1256332, -72.784558, -71.3406629])


def test_c_projections_shaped():
    nx, ny = (5, 2)
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x, y)

    model = projections.Pix2Sky_TAN()

    phi, theta = model(xv, yv)

    assert_allclose(
        phi,
        [[0., 90., 90., 90., 90.],
         [180., 165.96375653, 153.43494882, 143.13010235, 135.]])

    assert_allclose(
        theta,
        [[90., 89.75000159, 89.50001269, 89.25004283, 89.00010152],
         [89.00010152, 88.96933478, 88.88210788, 88.75019826, 88.58607353]])


def test_affine_with_quantities():
    x = 1
    y = 2
    xdeg = (x * u.pix).to(u.deg, equivalencies=u.pixel_scale(2.5 * u.deg / u.pix))
    ydeg = (y * u.pix).to(u.deg, equivalencies=u.pixel_scale(2.5 * u.deg / u.pix))
    xpix = x * u.pix
    ypix = y * u.pix

    # test affine with matrix only
    qaff = projections.AffineTransformation2D(matrix=[[1, 2], [2, 1]] * u.deg)
    with pytest.raises(ValueError):
        qx1, qy1 = qaff(xpix, ypix, equivalencies={
            'x': u.pixel_scale(2.5 * u.deg / u.pix),
            'y': u.pixel_scale(2.5 * u.deg / u.pix)})

    # test affine with matrix and translation
    qaff = projections.AffineTransformation2D(matrix=[[1, 2], [2, 1]] * u.deg,
                                              translation=[1, 2] * u.deg)
    qx1, qy1 = qaff(xpix, ypix, equivalencies={
        'x': u.pixel_scale(2.5 * u.deg / u.pix),
        'y': u.pixel_scale(2.5 * u.deg / u.pix)})
    aff = projections.AffineTransformation2D(matrix=[[1, 2], [2, 1]], translation=[1, 2])
    x1, y1 = aff(xdeg.value, ydeg.value)
    assert_quantity_allclose(qx1, x1 * u.deg)
    assert_quantity_allclose(qy1, y1 * u.deg)

    # test the case of WCS PC and CDELT transformations
    pc = np.array([[0.86585778922708, 0.50029020461607],
                   [-0.50029020461607, 0.86585778922708]])
    cdelt = np.array([[1, 3.0683055555556E-05], [3.0966944444444E-05, 1]])
    matrix = cdelt * pc
    qaff = projections.AffineTransformation2D(matrix=matrix * u.deg,
                                              translation=[0, 0] * u.deg)

    inv_matrix = np.linalg.inv(matrix)
    inv_qaff = projections.AffineTransformation2D(matrix=inv_matrix * u.pix,
                                                  translation=[0, 0] * u.pix)
    qaff.inverse = inv_qaff
    qx1, qy1 = qaff(xpix, ypix, equivalencies={
            'x': u.pixel_scale(1 * u.deg / u.pix),
            'y': u.pixel_scale(1 * u.deg / u.pix)})
    x1, y1 = qaff.inverse(qx1, qy1, equivalencies={
        'x': u.pixel_scale(1 * u.deg / u.pix),
        'y': u.pixel_scale(1 * u.deg / u.pix)})
    assert_quantity_allclose(x1, xpix)
    assert_quantity_allclose(y1, ypix)
