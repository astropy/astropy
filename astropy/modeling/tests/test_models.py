# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for model evaluation.
Compare the results of some models with other programs.
"""
# pylint: disable=invalid-name, no-member
import pytest
import numpy as np

from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.modeling import fitting, models
from astropy.modeling.models import Gaussian2D
from astropy.modeling.core import FittableModel
from astropy.modeling.parameters import Parameter
from astropy.modeling.polynomial import PolynomialBase
from astropy.utils import minversion
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils import NumpyRNGContext
from .example_models import models_1D, models_2D

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


@pytest.mark.skipif('not HAS_SCIPY')
def test_custom_model(amplitude=4, frequency=1):

    def sine_model(x, amplitude=4, frequency=1):
        """
        Model function
        """
        return amplitude * np.sin(2 * np.pi * frequency * x)

    def sine_deriv(x, amplitude=4, frequency=1):
        """
        Jacobian of model function, e.g. derivative of the function with
        respect to the *parameters*
        """
        da = np.sin(2 * np.pi * frequency * x)
        df = 2 * np.pi * x * amplitude * np.cos(2 * np.pi * frequency * x)
        return np.vstack((da, df))

    SineModel = models.custom_model(sine_model, fit_deriv=sine_deriv)

    x = np.linspace(0, 4, 50)
    sin_model = SineModel()

    sin_model.evaluate(x, 5., 2.)
    sin_model.fit_deriv(x, 5., 2.)

    np.random.seed(0)
    data = sin_model(x) + np.random.rand(len(x)) - 0.5
    fitter = fitting.LevMarLSQFitter()
    model = fitter(sin_model, x, data)
    assert np.all((np.array([model.amplitude.value, model.frequency.value]) -
                   np.array([amplitude, frequency])) < 0.001)


def test_custom_model_init():
    @models.custom_model
    def SineModel(x, amplitude=4, frequency=1):
        """Model function"""

        return amplitude * np.sin(2 * np.pi * frequency * x)

    sin_model = SineModel(amplitude=2., frequency=0.5)
    assert sin_model.amplitude == 2.
    assert sin_model.frequency == 0.5


def test_custom_model_defaults():
    @models.custom_model
    def SineModel(x, amplitude=4, frequency=1):
        """Model function"""

        return amplitude * np.sin(2 * np.pi * frequency * x)

    sin_model = SineModel()
    assert SineModel.amplitude.default == 4
    assert SineModel.frequency.default == 1

    assert sin_model.amplitude == 4
    assert sin_model.frequency == 1


def test_inconsistent_input_shapes():
    g = Gaussian2D()
    x = np.arange(-1., 1, .2)
    y = x.copy()
    # check scalar input broadcasting works
    assert np.abs(g(x, 0) - g(x, 0 * x)).sum() == 0
    # but not array broadcasting
    x.shape = (10, 1)
    y.shape = (1, 10)
    with pytest.raises(ValueError):
        g(x, y)


def test_custom_model_bounding_box():
    """Test bounding box evaluation for a 3D model"""

    def ellipsoid(x, y, z, x0=13, y0=10, z0=8, a=4, b=3, c=2, amp=1):
        rsq = ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 + ((z - z0) / c) ** 2
        val = (rsq < 1) * amp
        return val

    class Ellipsoid3D(models.custom_model(ellipsoid)):
        @property
        def bounding_box(self):
            return ((self.z0 - self.c, self.z0 + self.c),
                    (self.y0 - self.b, self.y0 + self.b),
                    (self.x0 - self.a, self.x0 + self.a))

    model = Ellipsoid3D()
    bbox = model.bounding_box

    zlim, ylim, xlim = bbox
    dz, dy, dx = np.diff(bbox) / 2
    z1, y1, x1 = np.mgrid[slice(zlim[0], zlim[1] + 1),
                          slice(ylim[0], ylim[1] + 1),
                          slice(xlim[0], xlim[1] + 1)]
    z2, y2, x2 = np.mgrid[slice(zlim[0] - dz, zlim[1] + dz + 1),
                          slice(ylim[0] - dy, ylim[1] + dy + 1),
                          slice(xlim[0] - dx, xlim[1] + dx + 1)]

    arr = model(x2, y2, z2)
    sub_arr = model(x1, y1, z1)

    # check for flux agreement
    assert abs(arr.sum() - sub_arr.sum()) < arr.sum() * 1e-7


class Fittable2DModelTester:
    """
    Test class for all two dimensional parametric models.

    Test values have to be defined in example_models.py. It currently test the
    model with different input types, evaluates the model at different
    positions and assures that it gives the correct values. And tests if the
    model works with non-linear fitters.

    This can be used as a base class for user defined model testing.
    """

    def setup_class(self):
        self.N = 100
        self.M = 100
        self.eval_error = 0.0001
        self.fit_error = 0.1
        self.x = 5.3
        self.y = 6.7
        self.x1 = np.arange(1, 10, .1)
        self.y1 = np.arange(1, 10, .1)
        self.y2, self.x2 = np.mgrid[:10, :8]

    def test_input2D(self, model_class, test_parameters):
        """Test model with different input types."""

        model = create_model(model_class, test_parameters)
        model(self.x, self.y)
        model(self.x1, self.y1)
        model(self.x2, self.y2)

    def test_eval2D(self, model_class, test_parameters):
        """Test model values add certain given points"""

        model = create_model(model_class, test_parameters)
        x = test_parameters['x_values']
        y = test_parameters['y_values']
        z = test_parameters['z_values']
        assert np.all(np.abs(model(x, y) - z) < self.eval_error)

    def test_bounding_box2D(self, model_class, test_parameters):
        """Test bounding box evaluation"""

        model = create_model(model_class, test_parameters)

        # testing setter
        model.bounding_box = ((-5, 5), (-5, 5))
        assert model.bounding_box == ((-5, 5), (-5, 5))

        model.bounding_box = None
        with pytest.raises(NotImplementedError):
            model.bounding_box

        # test the exception of dimensions don't match
        with pytest.raises(ValueError):
            model.bounding_box = (-5, 5)

        del model.bounding_box

        try:
            bbox = model.bounding_box
        except NotImplementedError:
            pytest.skip("Bounding_box is not defined for model.")

        ylim, xlim = bbox
        dy, dx = np.diff(bbox)/2
        y1, x1 = np.mgrid[slice(ylim[0], ylim[1] + 1),
                          slice(xlim[0], xlim[1] + 1)]
        y2, x2 = np.mgrid[slice(ylim[0] - dy, ylim[1] + dy + 1),
                          slice(xlim[0] - dx, xlim[1] + dx + 1)]

        arr = model(x2, y2)
        sub_arr = model(x1, y1)

        # check for flux agreement
        assert abs(arr.sum() - sub_arr.sum()) < arr.sum() * 1e-7

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_fitter2D(self, model_class, test_parameters):
        """Test if the parametric model works with the fitter."""

        x_lim = test_parameters['x_lim']
        y_lim = test_parameters['y_lim']

        parameters = test_parameters['parameters']
        model = create_model(model_class, test_parameters)

        if isinstance(parameters, dict):
            parameters = [parameters[name] for name in model.param_names]

        if "log_fit" in test_parameters:
            if test_parameters['log_fit']:
                x = np.logspace(x_lim[0], x_lim[1], self.N)
                y = np.logspace(y_lim[0], y_lim[1], self.N)
        else:
            x = np.linspace(x_lim[0], x_lim[1], self.N)
            y = np.linspace(y_lim[0], y_lim[1], self.N)
        xv, yv = np.meshgrid(x, y)

        np.random.seed(0)
        # add 10% noise to the amplitude
        noise = np.random.rand(self.N, self.N) - 0.5
        data = model(xv, yv) + 0.1 * parameters[0] * noise
        fitter = fitting.LevMarLSQFitter()
        new_model = fitter(model, xv, yv, data)

        params = [getattr(new_model, name) for name in new_model.param_names]
        fixed = [param.fixed for param in params]
        expected = np.array([val for val, fixed in zip(parameters, fixed)
                             if not fixed])
        fitted = np.array([param.value for param in params
                           if not param.fixed])
        assert_allclose(fitted, expected,
                        atol=self.fit_error)

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_deriv_2D(self, model_class, test_parameters):
        """
        Test the derivative of a model by fitting with an estimated and
        analytical derivative.
        """

        x_lim = test_parameters['x_lim']
        y_lim = test_parameters['y_lim']

        if model_class.fit_deriv is None:
            pytest.skip("Derivative function is not defined for model.")
        if issubclass(model_class, PolynomialBase):
            pytest.skip("Skip testing derivative of polynomials.")

        if "log_fit" in test_parameters:
            if test_parameters['log_fit']:
                x = np.logspace(x_lim[0], x_lim[1], self.N)
                y = np.logspace(y_lim[0], y_lim[1], self.M)
        else:
            x = np.linspace(x_lim[0], x_lim[1], self.N)
            y = np.linspace(y_lim[0], y_lim[1], self.M)
        xv, yv = np.meshgrid(x, y)

        try:
            model_with_deriv = create_model(model_class, test_parameters,
                                            use_constraints=False,
                                            parameter_key='deriv_initial')
            model_no_deriv = create_model(model_class, test_parameters,
                                          use_constraints=False,
                                          parameter_key='deriv_initial')
            model = create_model(model_class, test_parameters,
                                 use_constraints=False,
                                 parameter_key='deriv_initial')
        except KeyError:
            model_with_deriv = create_model(model_class, test_parameters,
                                            use_constraints=False)
            model_no_deriv = create_model(model_class, test_parameters,
                                          use_constraints=False)
            model = create_model(model_class, test_parameters,
                                 use_constraints=False)

        # add 10% noise to the amplitude
        rsn = np.random.RandomState(1234567890)
        amplitude = test_parameters['parameters'][0]
        n = 0.1 * amplitude * (rsn.rand(self.M, self.N) - 0.5)

        data = model(xv, yv) + n
        fitter_with_deriv = fitting.LevMarLSQFitter()
        new_model_with_deriv = fitter_with_deriv(model_with_deriv, xv, yv,
                                                 data)
        fitter_no_deriv = fitting.LevMarLSQFitter()
        new_model_no_deriv = fitter_no_deriv(model_no_deriv, xv, yv, data,
                                             estimate_jacobian=True)
        assert_allclose(new_model_with_deriv.parameters,
                        new_model_no_deriv.parameters,
                        rtol=0.1)


class Fittable1DModelTester:
    """
    Test class for all one dimensional parametric models.

    Test values have to be defined in example_models.py. It currently test the
    model with different input types, evaluates the model at different
    positions and assures that it gives the correct values. And tests if the
    model works with non-linear fitters.

    This can be used as a base class for user defined model testing.
    """

    def setup_class(self):
        self.N = 100
        self.M = 100
        self.eval_error = 0.0001
        self.fit_error = 0.1
        self.x = 5.3
        self.y = 6.7
        self.x1 = np.arange(1, 10, .1)
        self.y1 = np.arange(1, 10, .1)
        self.y2, self.x2 = np.mgrid[:10, :8]

    @pytest.mark.filterwarnings(r'ignore:.*:RuntimeWarning')
    def test_input1D(self, model_class, test_parameters):
        """Test model with different input types."""

        model = create_model(model_class, test_parameters)
        model(self.x)
        model(self.x1)
        model(self.x2)

    def test_eval1D(self, model_class, test_parameters):
        """
        Test model values at certain given points
        """
        model = create_model(model_class, test_parameters)
        x = test_parameters['x_values']
        y = test_parameters['y_values']
        assert_allclose(model(x), y, atol=self.eval_error)

    def test_bounding_box1D(self, model_class, test_parameters):
        """Test bounding box evaluation"""

        model = create_model(model_class, test_parameters)

        # testing setter
        model.bounding_box = (-5, 5)
        model.bounding_box = None

        with pytest.raises(NotImplementedError):
            model.bounding_box

        del model.bounding_box

        # test exception if dimensions don't match
        with pytest.raises(ValueError):
            model.bounding_box = 5

        try:
            bbox = model.bounding_box
        except NotImplementedError:
            pytest.skip("Bounding_box is not defined for model.")

        if isinstance(model, models.Lorentz1D) or isinstance(model, models.Drude1D):
            rtol = 0.01  # 1% agreement is enough due to very extended wings
            ddx = 0.1  # Finer sampling to "integrate" flux for narrow peak
        else:
            rtol = 1e-7
            ddx = 1

        dx = np.diff(bbox) / 2
        x1 = np.mgrid[slice(bbox[0], bbox[1] + 1, ddx)]
        x2 = np.mgrid[slice(bbox[0] - dx, bbox[1] + dx + 1, ddx)]
        arr = model(x2)
        sub_arr = model(x1)

        # check for flux agreement
        assert abs(arr.sum() - sub_arr.sum()) < arr.sum() * rtol

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_fitter1D(self, model_class, test_parameters):
        """
        Test if the parametric model works with the fitter.
        """
        x_lim = test_parameters['x_lim']
        parameters = test_parameters['parameters']
        model = create_model(model_class, test_parameters)

        if isinstance(parameters, dict):
            parameters = [parameters[name] for name in model.param_names]

        if "log_fit" in test_parameters:
            if test_parameters['log_fit']:
                x = np.logspace(x_lim[0], x_lim[1], self.N)
        else:
            x = np.linspace(x_lim[0], x_lim[1], self.N)

        np.random.seed(0)
        # add 10% noise to the amplitude
        relative_noise_amplitude = 0.01
        data = ((1 + relative_noise_amplitude * np.random.randn(len(x))) *
                model(x))
        fitter = fitting.LevMarLSQFitter()
        new_model = fitter(model, x, data)

        # Only check parameters that were free in the fit
        params = [getattr(new_model, name) for name in new_model.param_names]
        fixed = [param.fixed for param in params]
        expected = np.array([val for val, fixed in zip(parameters, fixed)
                             if not fixed])
        fitted = np.array([param.value for param in params
                           if not param.fixed])
        assert_allclose(fitted, expected, atol=self.fit_error)

    @pytest.mark.skipif('not HAS_SCIPY')
    @pytest.mark.filterwarnings(r'ignore:.*:RuntimeWarning')
    def test_deriv_1D(self, model_class, test_parameters):
        """
        Test the derivative of a model by comparing results with an estimated
        derivative.
        """

        x_lim = test_parameters['x_lim']

        if model_class.fit_deriv is None:
            pytest.skip("Derivative function is not defined for model.")
        if issubclass(model_class, PolynomialBase):
            pytest.skip("Skip testing derivative of polynomials.")

        if "log_fit" in test_parameters:
            if test_parameters['log_fit']:
                x = np.logspace(x_lim[0], x_lim[1], self.N)
        else:
            x = np.linspace(x_lim[0], x_lim[1], self.N)

        parameters = test_parameters['parameters']
        model_with_deriv = create_model(model_class, test_parameters,
                                        use_constraints=False)
        model_no_deriv = create_model(model_class, test_parameters,
                                      use_constraints=False)

        # add 10% noise to the amplitude
        rsn = np.random.RandomState(1234567890)
        n = 0.1 * parameters[0] * (rsn.rand(self.N) - 0.5)

        data = model_with_deriv(x) + n
        fitter_with_deriv = fitting.LevMarLSQFitter()
        new_model_with_deriv = fitter_with_deriv(model_with_deriv, x, data)
        fitter_no_deriv = fitting.LevMarLSQFitter()
        new_model_no_deriv = fitter_no_deriv(model_no_deriv, x, data,
                                             estimate_jacobian=True)
        assert_allclose(new_model_with_deriv.parameters,
                        new_model_no_deriv.parameters, atol=0.15)


def create_model(model_class, test_parameters, use_constraints=True,
                 parameter_key='parameters'):
    """Create instance of model class."""

    constraints = {}
    if issubclass(model_class, PolynomialBase):
        return model_class(**test_parameters[parameter_key])
    elif issubclass(model_class, FittableModel):
        if "requires_scipy" in test_parameters and not HAS_SCIPY:
            pytest.skip("SciPy not found")
        if use_constraints:
            if 'constraints' in test_parameters:
                constraints = test_parameters['constraints']
        return model_class(*test_parameters[parameter_key], **constraints)


@pytest.mark.filterwarnings(r'ignore:Model is linear in parameters.*')
@pytest.mark.filterwarnings(r'ignore:The fit may be unsuccessful.*')
@pytest.mark.parametrize(('model_class', 'test_parameters'),
                         sorted(models_1D.items(), key=lambda x: str(x[0])))
class TestFittable1DModels(Fittable1DModelTester):
    pass


@pytest.mark.filterwarnings(r'ignore:Model is linear in parameters.*')
@pytest.mark.parametrize(('model_class', 'test_parameters'),
                         sorted(models_2D.items(), key=lambda x: str(x[0])))
class TestFittable2DModels(Fittable2DModelTester):
    pass


def test_ShiftModel():
    # Shift by a scalar
    m = models.Shift(42)
    assert m(0) == 42
    assert_equal(m([1, 2]), [43, 44])

    # Shift by a list
    m = models.Shift([42, 43], n_models=2)
    assert_equal(m(0), [42, 43])
    assert_equal(m([1, 2], model_set_axis=False),
                 [[43, 44], [44, 45]])


def test_ScaleModel():
    # Scale by a scalar
    m = models.Scale(42)
    assert m(0) == 0
    assert_equal(m([1, 2]), [42, 84])

    # Scale by a list
    m = models.Scale([42, 43], n_models=2)
    assert_equal(m(0), [0, 0])
    assert_equal(m([1, 2], model_set_axis=False),
                 [[42, 84], [43, 86]])


def test_voigt_model():
    """
    Currently just tests that the model peaks at its origin.
    Regression test for https://github.com/astropy/astropy/issues/3942
    """

    m = models.Voigt1D(x_0=5, amplitude_L=10, fwhm_L=0.5, fwhm_G=0.9)
    x = np.arange(0, 10, 0.01)
    y = m(x)
    assert y[500] == y.max()  # y[500] is right at the center


def test_model_instance_repr():
    m = models.Gaussian1D(1.5, 2.5, 3.5)
    assert repr(m) == '<Gaussian1D(amplitude=1.5, mean=2.5, stddev=3.5)>'


@pytest.mark.skipif("not HAS_SCIPY")
def test_tabular_interp_1d():
    """
    Test Tabular1D model.
    """
    points = np.arange(0, 5)
    values = [1., 10, 2, 45, -3]
    LookupTable = models.tabular_model(1)
    model = LookupTable(points=points, lookup_table=values)
    xnew = [0., .7, 1.4, 2.1, 3.9]
    ans1 = [1., 7.3, 6.8, 6.3, 1.8]
    assert_allclose(model(xnew), ans1)
    # Test evaluate without passing `points`.
    model = LookupTable(lookup_table=values)
    assert_allclose(model(xnew), ans1)
    # Test bounds error.
    xextrap = [0., .7, 1.4, 2.1, 3.9, 4.1]
    with pytest.raises(ValueError):
        model(xextrap)
    # test extrapolation and fill value
    model = LookupTable(lookup_table=values, bounds_error=False,
                        fill_value=None)
    assert_allclose(model(xextrap),
                    [1., 7.3, 6.8, 6.3, 1.8, -7.8])

    # Test unit support
    xnew = xnew * u.nm
    ans1 = ans1 * u.nJy
    model = LookupTable(points=points*u.nm, lookup_table=values*u.nJy)
    assert_quantity_allclose(model(xnew), ans1)
    assert_quantity_allclose(model(xnew.to(u.nm)), ans1)
    assert model.bounding_box == (0 * u.nm, 4 * u.nm)

    # Test fill value unit conversion and unitless input on table with unit
    model = LookupTable([1, 2, 3], [10, 20, 30] * u.nJy, bounds_error=False,
                        fill_value=1e-33*(u.W / (u.m * u.m * u.Hz)))
    assert_quantity_allclose(model(np.arange(5)),
                             [100, 10, 20, 30, 100] * u.nJy)


@pytest.mark.skipif("not HAS_SCIPY")
def test_tabular_interp_2d():
    table = np.array([
        [-0.04614432, -0.02512547, -0.00619557, 0.0144165, 0.0297525],
        [-0.04510594, -0.03183369, -0.01118008, 0.01201388, 0.02496205],
        [-0.05464094, -0.02804499, -0.00960086, 0.01134333, 0.02284104],
        [-0.04879338, -0.02539565, -0.00440462, 0.01795145, 0.02122417],
        [-0.03637372, -0.01630025, -0.00157902, 0.01649774, 0.01952131]])

    points = np.arange(0, 5)
    points = (points, points)

    xnew = np.array([0., .7, 1.4, 2.1, 3.9])
    LookupTable = models.tabular_model(2)
    model = LookupTable(points, table)
    znew = model(xnew, xnew)
    result = np.array(
        [-0.04614432, -0.03450009, -0.02241028, -0.0069727, 0.01938675])
    assert_allclose(znew, result, atol=1e-7)

    # test 2D arrays as input
    a = np.arange(12).reshape((3, 4))
    y, x = np.mgrid[:3, :4]
    t = models.Tabular2D(lookup_table=a)
    r = t(y, x)
    assert_allclose(a, r)

    with pytest.raises(ValueError):
        model = LookupTable(points=([1.2, 2.3], [1.2, 6.7], [3, 4]))
    with pytest.raises(ValueError):
        model = LookupTable(lookup_table=[1, 2, 3])
    with pytest.raises(NotImplementedError):
        model = LookupTable(n_models=2)
    with pytest.raises(ValueError):
        model = LookupTable(([1, 2], [3, 4]), [5, 6])
    with pytest.raises(ValueError):
        model = LookupTable(([1, 2] * u.m, [3, 4]), [[5, 6], [7, 8]])
    with pytest.raises(ValueError):
        model = LookupTable(points, table, bounds_error=False,
                            fill_value=1*u.Jy)

    # Test unit support
    points = points[0] * u.nm
    points = (points, points)
    xnew = xnew * u.nm
    model = LookupTable(points, table * u.nJy)
    result = result * u.nJy
    assert_quantity_allclose(model(xnew, xnew), result, atol=1e-7*u.nJy)
    xnew = xnew.to(u.m)
    assert_quantity_allclose(model(xnew, xnew), result, atol=1e-7*u.nJy)
    bbox = (0 * u.nm, 4 * u.nm)
    bbox = (bbox, bbox)
    assert model.bounding_box == bbox


@pytest.mark.skipif("not HAS_SCIPY")
def test_tabular_nd():
    a = np.arange(24).reshape((2, 3, 4))
    x, y, z = np.mgrid[:2, :3, :4]
    tab = models.tabular_model(3)
    t = tab(lookup_table=a)
    result = t(x, y, z)
    assert_allclose(a, result)

    with pytest.raises(ValueError):
        models.tabular_model(0)


def test_with_bounding_box():
    """
    Test the option to evaluate a model respecting
    its bunding_box.
    """
    p = models.Polynomial2D(2) & models.Polynomial2D(2)
    m = models.Mapping((0, 1, 0, 1)) | p
    with NumpyRNGContext(1234567):
        m.parameters = np.random.rand(12)

    m.bounding_box = ((3, 9), (1, 8))
    x, y = np.mgrid[:10, :10]
    a, b = m(x, y)
    aw, bw = m(x, y, with_bounding_box=True)
    ind = (~np.isnan(aw)).nonzero()
    assert_allclose(a[ind], aw[ind])
    assert_allclose(b[ind], bw[ind])

    aw, bw = m(x, y, with_bounding_box=True, fill_value=1000)
    ind = (aw != 1000).nonzero()
    assert_allclose(a[ind], aw[ind])
    assert_allclose(b[ind], bw[ind])

    # test the order of bbox is not reversed for 1D models
    p = models.Polynomial1D(1, c0=12, c1=2.3)
    p.bounding_box = (0, 5)
    assert(p(1) == p(1, with_bounding_box=True))

    t3 = models.Shift(10) & models.Scale(2) & models.Shift(-1)
    t3.bounding_box = ((4.3, 6.9), (6, 15), (-1, 10))
    assert_allclose(t3([1, 1], [7, 7], [3, 5], with_bounding_box=True),
                    [[np.nan, 11], [np.nan, 14], [np.nan, 4]])

    trans3 = models.Shift(10) & models.Scale(2) & models.Shift(-1)
    trans3.bounding_box = ((4.3, 6.9), (6, 15), (-1, 10))
    assert_allclose(trans3(1, 7, 5, with_bounding_box=True), [11, 14, 4])


@pytest.mark.skipif("not HAS_SCIPY")
def test_tabular_with_bounding_box():
    points = np.arange(5)
    values = np.array([1.5, 3.4, 6.7, 7, 32])
    t = models.Tabular1D(points, values)
    result = t(1, with_bounding_box=True)

    assert result == 3.4
    assert t.inverse(result, with_bounding_box=True) == 1.


@pytest.mark.skipif("not HAS_SCIPY")
def test_tabular_bounding_box_with_units():
    points = np.arange(5)*u.pix
    lt = np.arange(5)*u.AA
    t = models.Tabular1D(points, lt)
    result = t(1*u.pix, with_bounding_box=True)

    assert result == 1.*u.AA
    assert t.inverse(result, with_bounding_box=True) == 1*u.pix


@pytest.mark.skipif("not HAS_SCIPY")
def test_tabular1d_inverse():
    """Test that the Tabular1D inverse is defined"""
    points = np.arange(5)
    values = np.array([1.5, 3.4, 6.7, 7, 32])
    t = models.Tabular1D(points, values)
    result = t.inverse((3.4, 6.7))
    assert_allclose(result, np.array((1., 2.)))

    # Check that it works for descending values in lookup_table
    t2 = models.Tabular1D(points, values[::-1])
    assert_allclose(t2.inverse.points[0], t2.lookup_table[::-1])

    result2 = t2.inverse((7, 6.7))
    assert_allclose(result2, np.array((1., 2.)))

    # Check that it errors on double-valued lookup_table
    points = np.arange(5)
    values = np.array([1.5, 3.4, 3.4, 32, 25])
    t = models.Tabular1D(points, values)
    with pytest.raises(NotImplementedError):
        t.inverse((3.4, 7.))

    # Check that Tabular2D.inverse raises an error
    table = np.arange(5*5).reshape(5, 5)
    points = np.arange(0, 5)
    points = (points, points)
    t3 = models.Tabular2D(points=points, lookup_table=table)
    with pytest.raises(NotImplementedError):
        t3.inverse((3, 3))

    # Check that it uses the same kwargs as the original model
    points = np.arange(5)
    values = np.array([1.5, 3.4, 6.7, 7, 32])
    t = models.Tabular1D(points, values)
    with pytest.raises(ValueError):
        t.inverse(100)
    t = models.Tabular1D(points, values, bounds_error=False, fill_value=None)
    result = t.inverse(100)
    assert_allclose(t(result), 100)


@pytest.mark.skipif("not HAS_SCIPY")
def test_tabular_module_name():
    """
    The module name must be set manually because
    these classes are created dynamically.
    """
    for model in [models.Tabular1D, models.Tabular2D]:
        assert model.__module__ == "astropy.modeling.tabular"


class classmodel(FittableModel):
    f = Parameter(default=1)
    x = Parameter(default=0)
    y = Parameter(default=2)

    def __init__(self, f=f.default, x=x.default, y=y.default):
        super().__init__(f, x, y)

    def evaluate(self):
        pass


class subclassmodel(classmodel):
    f = Parameter(default=3, fixed=True)
    x = Parameter(default=10)
    y = Parameter(default=12)
    h = Parameter(default=5)

    def __init__(self, f=f.default, x=x.default, y=y.default, h=h.default):
        super().__init__(f, x, y)

    def evaluate(self):
        pass


def test_parameter_inheritance():
    b = subclassmodel()
    assert b.param_names == ('f', 'x', 'y', 'h')
    assert b.h == 5
    assert b.f == 3
    assert b.f.fixed == True
