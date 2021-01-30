# Licensed under a 3-clause BSD style license - see LICENSE.rst

# pylint: disable=invalid-name
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal, assert_array_less

from astropy.modeling import models, InputParameterError
from astropy.coordinates import Angle
from astropy.modeling import fitting
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning

try:
    from scipy.integrate import quad  # pylint: disable=W0611  # noqa
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def test_sigma_constant():
    """
    Test that the GAUSSIAN_SIGMA_TO_FWHM constant matches the
    gaussian_sigma_to_fwhm constant in astropy.stats. We define
    it manually in astropy.modeling to avoid importing from
    astropy.stats.
    """
    from astropy.stats.funcs import gaussian_sigma_to_fwhm
    from astropy.modeling.functional_models import GAUSSIAN_SIGMA_TO_FWHM
    assert gaussian_sigma_to_fwhm == GAUSSIAN_SIGMA_TO_FWHM


def test_Trapezoid1D():
    """Regression test for https://github.com/astropy/astropy/issues/1721"""

    model = models.Trapezoid1D(amplitude=4.2, x_0=2.0, width=1.0, slope=3)
    xx = np.linspace(0, 4, 8)
    yy = model(xx)
    yy_ref = [0., 1.41428571, 3.12857143, 4.2, 4.2, 3.12857143, 1.41428571, 0.]
    assert_allclose(yy, yy_ref, rtol=0, atol=1e-6)


def test_Gaussian2D():
    """
    Test rotated elliptical Gaussian2D model.
    https://github.com/astropy/astropy/pull/2038
    """

    model = models.Gaussian2D(4.2, 1.7, 3.1, x_stddev=5.1, y_stddev=3.3,
                              theta=np.pi/6.)
    y, x = np.mgrid[0:5, 0:5]
    g = model(x, y)
    g_ref = [[3.01907812, 2.99051889, 2.81271552, 2.5119566, 2.13012709],
             [3.55982239, 3.6086023, 3.4734158, 3.17454575, 2.75494838],
             [3.88059142, 4.0257528, 3.96554926, 3.70908389, 3.29410187],
             [3.91095768, 4.15212857, 4.18567526, 4.00652015, 3.64146544],
             [3.6440466, 3.95922417, 4.08454159, 4.00113878, 3.72161094]]
    assert_allclose(g, g_ref, rtol=0, atol=1e-6)
    assert_allclose([model.x_fwhm, model.y_fwhm],
                    [12.009582229657841, 7.7709061486021325])


def test_Gaussian2DCovariance():
    """
    Test rotated elliptical Gaussian2D model when cov_matrix is input.
    https://github.com/astropy/astropy/pull/2199
    """

    cov_matrix = [[49., -16.], [-16., 9.]]
    model = models.Gaussian2D(17., 2.0, 2.5, cov_matrix=cov_matrix)
    y, x = np.mgrid[0:5, 0:5]
    g = model(x, y)
    g_ref = [[4.3744505, 5.8413977, 7.42988694, 9.00160175, 10.38794269],
             [8.83290201, 10.81772851, 12.61946384, 14.02225593, 14.84113227],
             [13.68528889, 15.37184621, 16.44637743, 16.76048705, 16.26953638],
             [16.26953638, 16.76048705, 16.44637743, 15.37184621, 13.68528889],
             [14.84113227, 14.02225593, 12.61946384, 10.81772851, 8.83290201]]
    assert_allclose(g, g_ref, rtol=0, atol=1e-6)


def test_Gaussian2DRotation():
    amplitude = 42
    x_mean, y_mean = 0, 0
    x_stddev, y_stddev = 2, 3
    theta = Angle(10, 'deg')
    pars = dict(amplitude=amplitude, x_mean=x_mean, y_mean=y_mean,
                x_stddev=x_stddev, y_stddev=y_stddev)
    rotation = models.Rotation2D(angle=theta.degree)
    point1 = (x_mean + 2 * x_stddev, y_mean + 2 * y_stddev)
    point2 = rotation(*point1)
    g1 = models.Gaussian2D(theta=0, **pars)
    g2 = models.Gaussian2D(theta=theta.radian, **pars)
    value1 = g1(*point1)
    value2 = g2(*point2)
    assert_allclose(value1, value2)


def test_Gaussian2D_invalid_inputs():
    x_stddev = 5.1
    y_stddev = 3.3
    theta = 10
    cov_matrix = [[49., -16.], [-16., 9.]]

    # first make sure the valid ones are OK
    models.Gaussian2D()
    models.Gaussian2D(x_stddev=x_stddev, y_stddev=y_stddev, theta=theta)
    models.Gaussian2D(x_stddev=None, y_stddev=y_stddev, theta=theta)
    models.Gaussian2D(x_stddev=x_stddev, y_stddev=None, theta=theta)
    models.Gaussian2D(x_stddev=x_stddev, y_stddev=y_stddev, theta=None)
    models.Gaussian2D(cov_matrix=cov_matrix)

    with pytest.raises(InputParameterError):
        models.Gaussian2D(x_stddev=0, cov_matrix=cov_matrix)
    with pytest.raises(InputParameterError):
        models.Gaussian2D(y_stddev=0, cov_matrix=cov_matrix)
    with pytest.raises(InputParameterError):
        models.Gaussian2D(theta=0, cov_matrix=cov_matrix)


@pytest.mark.parametrize('gamma', (10, -10))
def test_moffat_fwhm(gamma):
    ans = 34.641016151377542
    kwargs = {'gamma': gamma, 'alpha': 0.5}
    m1 = models.Moffat1D(**kwargs)
    m2 = models.Moffat2D(**kwargs)
    assert_allclose([m1.fwhm, m2.fwhm], ans)
    assert_array_less(0, [m1.fwhm, m2.fwhm])


def test_RedshiftScaleFactor():
    """Like ``test_ScaleModel()``."""

    # Scale by a scalar
    m = models.RedshiftScaleFactor(0.4)
    assert m(0) == 0
    assert_array_equal(m([1, 2]), [1.4, 2.8])

    assert_allclose(m.inverse(m([1, 2])), [1, 2])

    # Scale by a list
    m = models.RedshiftScaleFactor([-0.5, 0, 0.5], n_models=3)
    assert_array_equal(m(0), 0)
    assert_array_equal(m([1, 2], model_set_axis=False),
                       [[0.5, 1], [1, 2], [1.5, 3]])

    assert_allclose(m.inverse(m([1, 2], model_set_axis=False)),
                    [[1, 2], [1, 2], [1, 2]])


def test_Ellipse2D():
    """Test Ellipse2D model."""
    amplitude = 7.5
    x0, y0 = 15, 15
    theta = Angle(45, 'deg')
    em = models.Ellipse2D(amplitude, x0, y0, 7, 3, theta.radian)
    y, x = np.mgrid[0:30, 0:30]
    e = em(x, y)
    assert np.all(e[e > 0] == amplitude)
    assert e[y0, x0] == amplitude

    rotation = models.Rotation2D(angle=theta.degree)
    point1 = [2, 0]      # Rotation2D center is (0, 0)
    point2 = rotation(*point1)
    point1 = np.array(point1) + [x0, y0]
    point2 = np.array(point2) + [x0, y0]
    e1 = models.Ellipse2D(amplitude, x0, y0, 7, 3, theta=0.)
    e2 = models.Ellipse2D(amplitude, x0, y0, 7, 3, theta=theta.radian)
    assert e1(*point1) == e2(*point2)


def test_Ellipse2D_circular():
    """Test that circular Ellipse2D agrees with Disk2D [3736]."""
    amplitude = 7.5
    radius = 10
    size = (radius * 2) + 1
    y, x = np.mgrid[0:size, 0:size]
    ellipse = models.Ellipse2D(amplitude, radius, radius, radius, radius,
                               theta=0)(x, y)
    disk = models.Disk2D(amplitude, radius, radius, radius)(x, y)
    assert np.all(ellipse == disk)


def test_Scale_inverse():
    m = models.Scale(1.2345)
    assert_allclose(m.inverse(m(6.789)), 6.789)


def test_Multiply_inverse():
    m = models.Multiply(1.2345)
    assert_allclose(m.inverse(m(6.789)), 6.789)


def test_Shift_inverse():
    m = models.Shift(1.2345)
    assert_allclose(m.inverse(m(6.789)), 6.789)


@pytest.mark.skipif('not HAS_SCIPY')
def test_Shift_model_levmar_fit():
    """Test fitting Shift model with LevMarLSQFitter (issue #6103)."""

    init_model = models.Shift()

    x = np.arange(10)
    y = x+0.1

    fitter = fitting.LevMarLSQFitter()
    with pytest.warns(AstropyUserWarning,
                      match='Model is linear in parameters'):
        fitted_model = fitter(init_model, x, y)

    assert_allclose(fitted_model.parameters, [0.1], atol=1e-15)


def test_Shift_model_set_linear_fit():
    """Test linear fitting of Shift model (issue #6103)."""

    init_model = models.Shift(offset=[0, 0], n_models=2)

    x = np.arange(10)
    yy = np.array([x+0.1, x-0.2])

    fitter = fitting.LinearLSQFitter()
    fitted_model = fitter(init_model, x, yy)

    assert_allclose(fitted_model.parameters, [0.1, -0.2], atol=1e-15)


@pytest.mark.parametrize('Model', (models.Scale, models.Multiply))
def test_Scale_model_set_linear_fit(Model):
    """Test linear fitting of Scale model (#6103)."""

    init_model = Model(factor=[0, 0], n_models=2)

    x = np.arange(-3, 7)
    yy = np.array([1.15*x, 0.96*x])

    fitter = fitting.LinearLSQFitter()
    fitted_model = fitter(init_model, x, yy)

    assert_allclose(fitted_model.parameters, [1.15, 0.96], atol=1e-15)


# https://github.com/astropy/astropy/issues/6178
def test_Ring2D_rout():
    m = models.Ring2D(amplitude=1, x_0=1, y_0=1, r_in=2, r_out=5)
    assert m.width.value == 3


@pytest.mark.skipif("not HAS_SCIPY")
def test_Voigt1D():
    voi = models.Voigt1D(amplitude_L=-0.5, x_0=1.0, fwhm_L=5.0, fwhm_G=5.0)
    xarr = np.linspace(-5.0, 5.0, num=40)
    yarr = voi(xarr)
    voi_init = models.Voigt1D(amplitude_L=-1.0, x_0=1.0, fwhm_L=5.0, fwhm_G=5.0)
    fitter = fitting.LevMarLSQFitter()
    voi_fit = fitter(voi_init, xarr, yarr)
    assert_allclose(voi_fit.param_sets, voi.param_sets)


@pytest.mark.skipif("not HAS_SCIPY")
@pytest.mark.parametrize('algorithm', ('humlicek2', 'wofz'))
def test_Voigt1D_norm(algorithm):
    """Test integral of normalized Voigt profile."""
    voi = models.Voigt1D(amplitude_L=1.0/np.pi, x_0=0.0, fwhm_L=2.0, fwhm_G=1.5, method=algorithm)
    if algorithm == 'wofz':
        atol = 1e-14
    else:
        atol = 1e-8
    assert_allclose(quad(voi, -np.inf, np.inf)[0], 1.0, atol=atol)


@pytest.mark.skipif("not HAS_SCIPY")
@pytest.mark.parametrize('doppler', (1.e-3, 1.e-2, 0.1, 0.5, 1.0, 2.5, 5.0, 10))
def test_Voigt1D_hum2(doppler):
    """Verify accuracy of Voigt profile in Humlicek approximation to Faddeeva.cc (SciPy)."""
    x = np.linspace(-20, 20, 400001)

    voi_w = models.Voigt1D(amplitude_L=2.0/np.pi, fwhm_L=1.0, fwhm_G=doppler, method='wofz')
    vf_w = voi_w(x)
    dvda_w = voi_w.fit_deriv(x, x_0=0, amplitude_L=2.0/np.pi, fwhm_L=1.0, fwhm_G=doppler)

    voi_h = models.Voigt1D(amplitude_L=2.0/np.pi, fwhm_L=1.0, fwhm_G=doppler, method='humlicek2')
    vf_h = voi_h(x)
    dvda_h = voi_h.fit_deriv(x, x_0=0, amplitude_L=2.0/np.pi, fwhm_L=1.0, fwhm_G=doppler)

    assert_allclose(vf_h, vf_w, rtol=1e-7 * (2 + 1 / np.sqrt(doppler)))
    assert_allclose(dvda_h, dvda_w, rtol=1e-9, atol=1e-7 * (1 + 30 / doppler))


@pytest.mark.skipif("not HAS_SCIPY")
def test_KingProjectedAnalytic1D_fit():
    km = models.KingProjectedAnalytic1D(amplitude=1, r_core=1, r_tide=2)
    xarr = np.linspace(0.1, 2, 10)
    yarr = km(xarr)
    km_init = models.KingProjectedAnalytic1D(amplitude=1, r_core=1, r_tide=1)
    fitter = fitting.LevMarLSQFitter()
    km_fit = fitter(km_init, xarr, yarr)
    assert_allclose(km_fit.param_sets, km.param_sets)


def test_ExponentialAndLogarithmic1D_fit():
    xarr = np.linspace(0.1, 10., 200)
    em_model = models.Exponential1D(amplitude=1, tau=1)
    log_model = models.Logarithmic1D(amplitude=1, tau=1)
    assert_allclose(xarr, em_model.inverse(em_model(xarr)))
    assert_allclose(xarr, log_model.inverse(log_model(xarr)))


def test_deprecated_hat_kernel():

    # 'MexicanHat' was deprecated as a name for the models which are now
    # 'RickerWavelet'. This test ensures that the models are correctly
    # deprecated.

    with pytest.warns(AstropyDeprecationWarning):
        models.MexicanHat1D()

    with pytest.warns(AstropyDeprecationWarning):
        models.MexicanHat2D()
