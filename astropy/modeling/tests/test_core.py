# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
import os
import subprocess
import sys
import unittest.mock as mk
from inspect import signature

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

import astropy
import astropy.modeling.core as core
import astropy.units as u
from astropy.convolution import convolve_models
from astropy.modeling import models
from astropy.modeling.bounding_box import CompoundBoundingBox, ModelBoundingBox
from astropy.modeling.core import (
    SPECIAL_OPERATORS,
    CompoundModel,
    Model,
    _add_special_operator,
    bind_bounding_box,
    bind_compound_bounding_box,
    custom_model,
    fix_inputs,
)
from astropy.modeling.parameters import Parameter
from astropy.modeling.separable import separability_matrix
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY


class NonFittableModel(Model):
    """An example class directly subclassing Model for testing."""

    a = Parameter()

    def __init__(self, a, model_set_axis=None):
        super().__init__(a, model_set_axis=model_set_axis)

    @staticmethod
    def evaluate():
        pass


def test_Model_instance_repr_and_str():
    m = NonFittableModel(42.5)
    assert repr(m) == "<NonFittableModel(a=42.5)>"
    assert (
        str(m) == "Model: NonFittableModel\n"
        "Inputs: ()\n"
        "Outputs: ()\n"
        "Model set size: 1\n"
        "Parameters:\n"
        "     a  \n"
        "    ----\n"
        "    42.5"
    )

    assert len(m) == 1


def test_Model_array_parameter():
    model = models.Gaussian1D(4, 2, 1)
    assert_allclose(model.param_sets, [[4], [2], [1]])


def test_inputless_model():
    """
    Regression test for
    https://github.com/astropy/astropy/pull/3772#issuecomment-101821641
    """

    class TestModel(Model):
        n_outputs = 1
        a = Parameter()

        @staticmethod
        def evaluate(a):
            return a

    m = TestModel(1)
    assert m.a == 1
    assert m() == 1

    # Test array-like output
    m = TestModel([1, 2, 3], model_set_axis=False)
    assert len(m) == 1
    assert np.all(m() == [1, 2, 3])

    # Test a model set
    m = TestModel(a=[1, 2, 3], model_set_axis=0)
    assert len(m) == 3
    assert np.all(m() == [1, 2, 3])

    # Test a model set
    m = TestModel(a=[[1, 2, 3], [4, 5, 6]], model_set_axis=0)
    assert len(m) == 2
    assert np.all(m() == [[1, 2, 3], [4, 5, 6]])

    # Test a model set
    m = TestModel(a=[[1, 2, 3], [4, 5, 6]], model_set_axis=np.int64(0))
    assert len(m) == 2
    assert np.all(m() == [[1, 2, 3], [4, 5, 6]])


def test_ParametericModel():
    MESSAGE = r"Gaussian1D.__init__.* got an unrecognized parameter 'wrong'"
    with pytest.raises(TypeError, match=MESSAGE):
        models.Gaussian1D(1, 2, 3, wrong=4)


def test_custom_model_signature():
    """
    Tests that the signatures for the __init__ and __call__
    methods of custom models are useful.
    """

    @custom_model
    def model_a(x):
        return x

    assert model_a.param_names == ()
    assert model_a.n_inputs == 1
    sig = signature(model_a.__init__)
    assert list(sig.parameters.keys()) == ["self", "args", "meta", "name", "kwargs"]
    sig = signature(model_a.__call__)
    assert list(sig.parameters.keys()) == [
        "self",
        "inputs",
        "model_set_axis",
        "with_bounding_box",
        "fill_value",
        "equivalencies",
        "inputs_map",
        "new_inputs",
    ]

    @custom_model
    def model_b(x, a=1, b=2):
        return x + a + b

    assert model_b.param_names == ("a", "b")
    assert model_b.n_inputs == 1
    sig = signature(model_b.__init__)
    assert list(sig.parameters.keys()) == ["self", "a", "b", "kwargs"]
    assert [x.default for x in sig.parameters.values()] == [sig.empty, 1, 2, sig.empty]
    sig = signature(model_b.__call__)
    assert list(sig.parameters.keys()) == [
        "self",
        "inputs",
        "model_set_axis",
        "with_bounding_box",
        "fill_value",
        "equivalencies",
        "inputs_map",
        "new_inputs",
    ]

    @custom_model
    def model_c(x, y, a=1, b=2):
        return x + y + a + b

    assert model_c.param_names == ("a", "b")
    assert model_c.n_inputs == 2
    sig = signature(model_c.__init__)
    assert list(sig.parameters.keys()) == ["self", "a", "b", "kwargs"]
    assert [x.default for x in sig.parameters.values()] == [sig.empty, 1, 2, sig.empty]
    sig = signature(model_c.__call__)
    assert list(sig.parameters.keys()) == [
        "self",
        "inputs",
        "model_set_axis",
        "with_bounding_box",
        "fill_value",
        "equivalencies",
        "inputs_map",
        "new_inputs",
    ]


def test_custom_model_subclass():
    """Test that custom models can be subclassed."""

    @custom_model
    def model_a(x, a=1):
        return x * a

    class model_b(model_a):
        # Override the evaluate from model_a
        @classmethod
        def evaluate(cls, x, a):
            return -super().evaluate(x, a)

    b = model_b()
    assert b.param_names == ("a",)
    assert b.a == 1
    assert b(1) == -1

    sig = signature(model_b.__init__)
    assert list(sig.parameters.keys()) == ["self", "a", "kwargs"]
    sig = signature(model_b.__call__)
    assert list(sig.parameters.keys()) == [
        "self",
        "inputs",
        "model_set_axis",
        "with_bounding_box",
        "fill_value",
        "equivalencies",
        "inputs_map",
        "new_inputs",
    ]


def test_custom_model_parametrized_decorator():
    """Tests using custom_model as a decorator with parameters."""

    def cosine(x, amplitude=1):
        return [amplitude * np.cos(x)]

    @custom_model(fit_deriv=cosine)
    def sine(x, amplitude=1):
        return amplitude * np.sin(x)

    assert issubclass(sine, Model)
    s = sine(2)
    assert_allclose(s(np.pi / 2), 2)
    assert_allclose(s.fit_deriv(0, 2), 2)


def test_custom_model_n_outputs():
    """
    Test creating a custom_model which has more than one output, which
    requires special handling.
        Demonstrates issue #11791's ``n_outputs`` error has been solved
    """

    @custom_model
    def model(x, y, n_outputs=2):
        return x + 1, y + 1

    m = model()
    assert not isinstance(m.n_outputs, Parameter)
    assert isinstance(m.n_outputs, int)
    assert m.n_outputs == 2
    assert m.outputs == ("x0", "x1")
    assert (
        separability_matrix(m)
        == [
            [True, True],
            [True, True],
        ]
    ).all()

    @custom_model
    def model(x, y, z, n_outputs=3):
        return x + 1, y + 1, z + 1

    m = model()
    assert not isinstance(m.n_outputs, Parameter)
    assert isinstance(m.n_outputs, int)
    assert m.n_outputs == 3
    assert m.outputs == ("x0", "x1", "x2")
    assert (
        separability_matrix(m)
        == [
            [True, True, True],
            [True, True, True],
            [True, True, True],
        ]
    ).all()


def test_custom_model_settable_parameters():
    """
    Test creating a custom_model which specifically sets adjustable model
    parameters.
        Demonstrates part of issue #11791's notes about what passed parameters
        should/shouldn't be allowed. In this case, settable parameters
        should be allowed to have defaults set.
    """

    @custom_model
    def model(x, y, n_outputs=2, bounding_box=((1, 2), (3, 4))):
        return x + 1, y + 1

    m = model()
    assert m.n_outputs == 2
    assert m.bounding_box == ((1, 2), (3, 4))
    m.bounding_box = ((9, 10), (11, 12))
    assert m.bounding_box == ((9, 10), (11, 12))
    m = model(bounding_box=((5, 6), (7, 8)))
    assert m.n_outputs == 2
    assert m.bounding_box == ((5, 6), (7, 8))
    m.bounding_box = ((9, 10), (11, 12))
    assert m.bounding_box == ((9, 10), (11, 12))

    @custom_model
    def model(x, y, n_outputs=2, outputs=("z0", "z1")):
        return x + 1, y + 1

    m = model()
    assert m.n_outputs == 2
    assert m.outputs == ("z0", "z1")
    m.outputs = ("a0", "a1")
    assert m.outputs == ("a0", "a1")
    m = model(outputs=("w0", "w1"))
    assert m.n_outputs == 2
    assert m.outputs == ("w0", "w1")
    m.outputs = ("a0", "a1")
    assert m.outputs == ("a0", "a1")


def test_custom_model_regected_parameters():
    """
    Test creating a custom_model which attempts to override non-overridable
    parameters.
        Demonstrates part of issue #11791's notes about what passed parameters
        should/shouldn't be allowed. In this case, non-settable parameters
        should raise an error (unexpected behavior may occur).
    """

    with pytest.raises(
        ValueError, match=r"Parameter 'n_inputs' cannot be a model property: *"
    ):

        @custom_model
        def model1(x, y, n_outputs=2, n_inputs=3):
            return x + 1, y + 1

    with pytest.raises(
        ValueError, match=r"Parameter 'uses_quantity' cannot be a model property: *"
    ):

        @custom_model
        def model2(x, y, n_outputs=2, uses_quantity=True):
            return x + 1, y + 1


def test_custom_inverse():
    """Test setting a custom inverse on a model."""

    p = models.Polynomial1D(1, c0=-2, c1=3)
    # A trivial inverse for a trivial polynomial
    inv = models.Polynomial1D(1, c0=(2.0 / 3.0), c1=(1.0 / 3.0))

    MESSAGE = (
        r"No analytical or user-supplied inverse transform has been implemented for"
        r" this model"
    )
    with pytest.raises(NotImplementedError, match=MESSAGE):
        p.inverse

    p.inverse = inv

    x = np.arange(100)

    assert_allclose(x, p(p.inverse(x)))
    assert_allclose(x, p.inverse(p(x)))

    p.inverse = None

    with pytest.raises(NotImplementedError, match=MESSAGE):
        p.inverse


def test_custom_inverse_reset():
    """Test resetting a custom inverse to the model's default inverse."""

    class TestModel(Model):
        n_inputs = 0
        outputs = ("y",)

        @property
        def inverse(self):
            return models.Shift()

        @staticmethod
        def evaluate():
            return 0

    # The above test model has no meaning, nor does its inverse--this just
    # tests that setting an inverse and resetting to the default inverse works

    m = TestModel()
    assert isinstance(m.inverse, models.Shift)

    m.inverse = models.Scale()
    assert isinstance(m.inverse, models.Scale)

    del m.inverse
    assert isinstance(m.inverse, models.Shift)


def test_render_model_2d():
    imshape = (71, 141)
    image = np.zeros(imshape)
    coords = y, x = np.indices(imshape)

    model = models.Gaussian2D(x_stddev=6.1, y_stddev=3.9, theta=np.pi / 3)

    # test points for edges
    ye, xe = [0, 35, 70], [0, 70, 140]
    # test points for floating point positions
    yf, xf = [35.1, 35.5, 35.9], [70.1, 70.5, 70.9]

    test_pts = [(a, b) for a in xe for b in ye]
    test_pts += [(a, b) for a in xf for b in yf]

    for x0, y0 in test_pts:
        model.x_mean = x0
        model.y_mean = y0
        expected = model(x, y)
        for xy in [coords, None]:
            for im in [image.copy(), None]:
                if (im is None) & (xy is None):
                    # this case is tested in Fittable2DModelTester
                    continue
                actual = model.render(out=im, coords=xy)
                if im is None:
                    assert_allclose(actual, model.render(coords=xy))
                # assert images match
                assert_allclose(expected, actual, atol=3e-7)
                # assert model fully captured
                if (x0, y0) == (70, 35):
                    boxed = model.render()
                    flux = np.sum(expected)
                    assert ((flux - np.sum(boxed)) / flux) < 1e-7
    # test an error is raised when the bounding box is larger than the input array
    try:
        actual = model.render(out=np.zeros((1, 1)))
    except ValueError:
        pass


def test_render_model_1d():
    npix = 101
    image = np.zeros(npix)
    coords = np.arange(npix)

    model = models.Gaussian1D()

    # test points
    test_pts = [0, 49.1, 49.5, 49.9, 100]

    # test widths
    test_stdv = np.arange(5.5, 6.7, 0.2)

    for x0, stdv in [(p, s) for p in test_pts for s in test_stdv]:
        model.mean = x0
        model.stddev = stdv
        expected = model(coords)
        for x in [coords, None]:
            for im in [image.copy(), None]:
                if (im is None) & (x is None):
                    # this case is tested in Fittable1DModelTester
                    continue
                actual = model.render(out=im, coords=x)
                # assert images match
                assert_allclose(expected, actual, atol=3e-7)
                # assert model fully captured
                if (x0, stdv) == (49.5, 5.5):
                    boxed = model.render()
                    flux = np.sum(expected)
                    assert ((flux - np.sum(boxed)) / flux) < 1e-7


@pytest.mark.filterwarnings("ignore:invalid value encountered in less")
def test_render_model_3d():
    imshape = (17, 21, 27)
    image = np.zeros(imshape)
    coords = np.indices(imshape)

    def ellipsoid(x, y, z, x0=13.0, y0=10.0, z0=8.0, a=4.0, b=3.0, c=2.0, amp=1.0):
        rsq = ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 + ((z - z0) / c) ** 2
        val = (rsq < 1) * amp
        return val

    class Ellipsoid3D(custom_model(ellipsoid)):
        @property
        def bounding_box(self):
            return (
                (self.z0 - self.c, self.z0 + self.c),
                (self.y0 - self.b, self.y0 + self.b),
                (self.x0 - self.a, self.x0 + self.a),
            )

    model = Ellipsoid3D()

    # test points for edges
    ze, ye, xe = [0, 8, 16], [0, 10, 20], [0, 13, 26]
    # test points for floating point positions
    zf, yf, xf = [8.1, 8.5, 8.9], [10.1, 10.5, 10.9], [13.1, 13.5, 13.9]

    test_pts = [(x, y, z) for x in xe for y in ye for z in ze]
    test_pts += [(x, y, z) for x in xf for y in yf for z in zf]

    for x0, y0, z0 in test_pts:
        model.x0 = x0
        model.y0 = y0
        model.z0 = z0
        expected = model(*coords[::-1])
        for c in [coords, None]:
            for im in [image.copy(), None]:
                if (im is None) & (c is None):
                    continue
                actual = model.render(out=im, coords=c)
                boxed = model.render()
                # assert images match
                assert_allclose(expected, actual)
                # assert model fully captured
                if (z0, y0, x0) == (8, 10, 13):
                    boxed = model.render()
                    assert (np.sum(expected) - np.sum(boxed)) == 0


def test_render_model_out_dtype():
    """Test different out.dtype for model.render."""
    MESSAGE = (
        r"Cannot cast ufunc 'add' output from .* to .* with casting rule 'same_kind"
    )
    for model in [models.Gaussian2D(), models.Gaussian2D() + models.Planar2D()]:
        for dtype in [np.float64, np.float32, np.complex64]:
            im = np.zeros((40, 40), dtype=dtype)
            imout = model.render(out=im)
            assert imout is im
            assert imout.sum() != 0

        with pytest.raises(TypeError, match=MESSAGE):
            im = np.zeros((40, 40), dtype=np.int32)
            imout = model.render(out=im)


def test_custom_bounding_box_1d():
    """
    Tests that the bounding_box setter works.
    """
    # 1D models
    g1 = models.Gaussian1D()
    bb = g1.bounding_box
    expected = g1.render()

    # assign the same bounding_box, now through the bounding_box setter
    g1.bounding_box = bb
    assert_allclose(g1.render(), expected)

    # 2D models
    g2 = models.Gaussian2D()
    bb = g2.bounding_box
    expected = g2.render()

    # assign the same bounding_box, now through the bounding_box setter
    g2.bounding_box = bb
    assert_allclose(g2.render(), expected)


def test_n_submodels_in_single_models():
    assert models.Gaussian1D().n_submodels == 1
    assert models.Gaussian2D().n_submodels == 1


def test_compound_deepcopy():
    model = (models.Gaussian1D(10, 2, 3) | models.Shift(2)) & models.Rotation2D(21.3)
    new_model = model.deepcopy()
    assert id(model) != id(new_model)
    assert id(model._leaflist) != id(new_model._leaflist)
    assert id(model[0]) != id(new_model[0])
    assert id(model[1]) != id(new_model[1])
    assert id(model[2]) != id(new_model[2])


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_units_with_bounding_box():
    points = np.arange(10, 20)
    table = np.arange(10) * u.Angstrom
    t = models.Tabular1D(points, lookup_table=table)

    assert isinstance(t(10), u.Quantity)
    assert isinstance(t(10, with_bounding_box=True), u.Quantity)

    assert_quantity_allclose(t(10), t(10, with_bounding_box=True))


RENAMED_MODEL = models.Gaussian1D.rename("CustomGaussian")

MODEL_RENAME_CODE = """
from astropy.modeling.models import Gaussian1D
print(repr(Gaussian1D))
print(repr(Gaussian1D.rename('CustomGaussian')))
""".strip()

MODEL_RENAME_EXPECTED = b"""
<class 'astropy.modeling.functional_models.Gaussian1D'>
Name: Gaussian1D
N_inputs: 1
N_outputs: 1
Fittable parameters: ('amplitude', 'mean', 'stddev')
<class '__main__.CustomGaussian'>
Name: CustomGaussian (Gaussian1D)
N_inputs: 1
N_outputs: 1
Fittable parameters: ('amplitude', 'mean', 'stddev')
""".strip()


def test_rename_path(tmp_path):
    # Regression test for a bug that caused the path to the class to be
    # incorrect in a renamed model's __repr__.

    assert (
        repr(RENAMED_MODEL).splitlines()[0]
        == "<class 'astropy.modeling.tests.test_core.CustomGaussian'>"
    )

    # Make sure that when called from a user script, the class name includes
    # __main__.

    env = os.environ.copy()
    paths = [os.path.dirname(astropy.__path__[0])] + sys.path
    env["PYTHONPATH"] = os.pathsep.join(paths)

    script = tmp_path / "rename.py"
    with open(script, "w") as f:
        f.write(MODEL_RENAME_CODE)

    output = subprocess.check_output([sys.executable, script], env=env)
    assert output.splitlines() == MODEL_RENAME_EXPECTED.splitlines()


@pytest.mark.parametrize(
    "model_class",
    [models.Gaussian1D, models.Polynomial1D, models.Shift, models.Tabular1D],
)
def test_rename_1d(model_class):
    new_model = model_class.rename(name="Test1D")
    assert new_model.name == "Test1D"


@pytest.mark.parametrize(
    "model_class", [models.Gaussian2D, models.Polynomial2D, models.Tabular2D]
)
def test_rename_2d(model_class):
    new_model = model_class.rename(name="Test2D")
    assert new_model.name == "Test2D"


def test_fix_inputs_integer():
    """
    Tests that numpy integers can be passed as dictionary keys to fix_inputs
    Issue #11358
    """
    m = models.Identity(2)

    mf = models.fix_inputs(m, {1: 22})
    assert mf(1) == (1, 22)

    mf_int32 = models.fix_inputs(m, {np.int32(1): 33})
    assert mf_int32(1) == (1, 33)

    mf_int64 = models.fix_inputs(m, {np.int64(1): 44})
    assert mf_int64(1) == (1, 44)


def test_fix_inputs_empty_dict():
    """
    Tests that empty dictionary can be passed to fix_inputs
    Issue #11355
    """
    m = models.Identity(2)

    mf = models.fix_inputs(m, {})
    assert mf(1, 2) == (1, 2)


def test_rename_inputs_outputs():
    g2 = models.Gaussian2D(10, 2, 3, 1, 2)
    assert g2.inputs == ("x", "y")
    assert g2.outputs == ("z",)

    MESSAGE = r"Expected .* number of .*, got .*"
    with pytest.raises(ValueError, match=MESSAGE):
        g2.inputs = ("w",)

    with pytest.raises(ValueError, match=MESSAGE):
        g2.outputs = ("w", "e")


def test__prepare_output_single_model():
    model = models.Gaussian1D()

    # No broadcast
    assert (
        np.array([1, 2]) == model._prepare_output_single_model(np.array([1, 2]), None)
    ).all()

    # Broadcast to scalar
    assert model._prepare_output_single_model(np.array([1]), ()) == 1
    assert model._prepare_output_single_model(np.asanyarray(2), ()) == 2

    # Broadcast reshape
    output = np.array([[1, 2, 3], [4, 5, 6]])
    reshape = np.array([[1, 2], [3, 4], [5, 6]])
    assert (output == model._prepare_output_single_model(output, (2, 3))).all()
    assert (reshape == model._prepare_output_single_model(output, (3, 2))).all()

    # Broadcast reshape scalar
    assert model._prepare_output_single_model(np.array([1]), (1, 2)) == 1
    assert model._prepare_output_single_model(np.asanyarray(2), (3, 4)) == 2

    # Fail to broadcast
    assert (output == model._prepare_output_single_model(output, (1, 2))).all()
    assert (output == model._prepare_output_single_model(output, (3, 4))).all()


def test_prepare_outputs_mixed_broadcast():
    """
    Tests that _prepare_outputs_single_model does not fail when a smaller
    array is passed as first input, but output is broadcast to larger
    array.
    Issue #10170
    """

    model = models.Gaussian2D(1, 2, 3, 4, 5)

    output = model([1, 2], 3)
    assert output.shape == (2,)
    np.testing.assert_array_equal(output, [0.9692332344763441, 1.0])

    output = model(4, [5, 6])
    assert output.shape == (2,)
    np.testing.assert_allclose(output, [0.8146473164114145, 0.7371233743916278])


def test_prepare_outputs_complex_reshape():
    x = np.array(
        [
            [1, 2, 3, 4, 5],
            [6, 7, 8, 9, 10],
            [11, 12, 13, 14, 15],
        ]
    )
    y = np.array(
        [
            [16, 17, 18, 19, 20],
            [21, 22, 23, 24, 25],
            [26, 27, 28, 29, 30],
        ]
    )

    m = models.Identity(3) | models.Mapping((2, 1, 0))
    m.bounding_box = ((0, 100), (0, 200), (0, 50))
    mf = models.fix_inputs(m, {2: 22})
    t = mf | models.Mapping((2, 1), n_inputs=3)

    output = mf(1, 2)
    assert output == (22, 2, 1)

    output = t(1, 2)
    assert output == (1, 2)

    output = t(x, y)
    assert len(output) == 2
    np.testing.assert_array_equal(output[0], x)
    np.testing.assert_array_equal(output[1], y)

    m = models.Identity(3) | models.Mapping((0, 1, 2))
    m.bounding_box = ((0, 100), (0, 200), (0, 50))
    mf = models.fix_inputs(m, {2: 22})
    t = mf | models.Mapping((0, 1), n_inputs=3)

    output = mf(1, 2)
    assert output == (1, 2, 22)

    output = t(1, 2)
    assert output == (1, 2)

    output = t(x, y)
    assert len(output) == 2
    np.testing.assert_array_equal(output[0], x)
    np.testing.assert_array_equal(output[1], y)


def test_prepare_outputs_single_entry_vector():
    """
    jwst and gwcs both require that single entry vectors produce single
    entry output vectors, not scalars. This tests for that behavior.
    """

    model = models.Gaussian2D(1, 2, 3, 4, 5)

    output = model(np.array([1]), np.array([2]))
    assert output.shape == (1,)
    np.testing.assert_allclose(output, [0.9500411305585278])


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.filterwarnings("ignore: Using a non-tuple")
def test_prepare_outputs_sparse_grid():
    """
    Test to show that #11060 has been solved.
    """

    shape = (3, 3)
    data = np.arange(np.prod(shape)).reshape(shape) * u.m / u.s

    points_unit = u.pix
    points = [np.arange(size) * points_unit for size in shape]

    kwargs = {
        "bounds_error": False,
        "fill_value": np.nan,
        "method": "nearest",
    }

    transform = models.Tabular2D(points, data, **kwargs)
    truth = (
        np.array(
            [
                [0.0, 1.0, 2.0],
                [3.0, 4.0, 5.0],
                [6.0, 7.0, 8.0],
            ]
        )
        * u.m
        / u.s
    )

    points = np.meshgrid(np.arange(3), np.arange(3), indexing="ij", sparse=True)
    x = points[0] * u.pix
    y = points[1] * u.pix
    value = transform(x, y)
    assert (value == truth).all()

    points = (
        np.meshgrid(np.arange(3), np.arange(3), indexing="ij", sparse=False) * u.pix
    )
    value = transform(*points)
    assert (value == truth).all()


def test_coerce_units():
    model = models.Polynomial1D(1, c0=1, c1=2)

    MESSAGE = r"Can only apply 'add' function to dimensionless quantities when other .*"
    with pytest.raises(u.UnitsError, match=MESSAGE):
        model(u.Quantity(10, u.m))

    with_input_units = model.coerce_units({"x": u.m})
    result = with_input_units(u.Quantity(10, u.m))
    assert np.isclose(result, 21.0)

    with_input_units_tuple = model.coerce_units((u.m,))
    result = with_input_units_tuple(u.Quantity(10, u.m))
    assert np.isclose(result, 21.0)

    with_return_units = model.coerce_units(return_units={"y": u.s})
    result = with_return_units(10)
    assert np.isclose(result.value, 21.0)
    assert result.unit == u.s

    with_return_units_tuple = model.coerce_units(return_units=(u.s,))
    result = with_return_units_tuple(10)
    assert np.isclose(result.value, 21.0)
    assert result.unit == u.s

    with_both = model.coerce_units({"x": u.m}, {"y": u.s})

    result = with_both(u.Quantity(10, u.m))
    assert np.isclose(result.value, 21.0)
    assert result.unit == u.s

    with pytest.raises(
        ValueError, match=r"input_units keys.*do not match model inputs"
    ):
        model.coerce_units({"q": u.m})

    with pytest.raises(ValueError, match=r"input_units length does not match n_inputs"):
        model.coerce_units((u.m, u.s))

    model_with_existing_input_units = models.BlackBody()
    with pytest.raises(
        ValueError,
        match=r"Cannot specify input_units for model with existing input units",
    ):
        model_with_existing_input_units.coerce_units({"x": u.m})

    with pytest.raises(
        ValueError, match=r"return_units keys.*do not match model outputs"
    ):
        model.coerce_units(return_units={"q": u.m})

    with pytest.raises(
        ValueError, match=r"return_units length does not match n_outputs"
    ):
        model.coerce_units(return_units=(u.m, u.s))


def test_bounding_box_general_inverse():
    model = NonFittableModel(42.5)

    MESSAGE = r"No bounding box is defined for this model"
    with pytest.raises(NotImplementedError, match=MESSAGE):
        model.bounding_box
    model.bounding_box = ()
    assert model.bounding_box.bounding_box() == ()

    model.inverse = NonFittableModel(3.14)
    inverse_model = model.inverse
    with pytest.raises(NotImplementedError, match=MESSAGE):
        inverse_model.bounding_box


def test__add_special_operator():
    sop_name = "name"
    sop = "value"

    key = _add_special_operator(sop_name, "value")
    assert key[0] == sop_name
    assert key[1] == SPECIAL_OPERATORS._unique_id

    assert key in SPECIAL_OPERATORS
    assert SPECIAL_OPERATORS[key] == sop


def test_print_special_operator_CompoundModel(capsys):
    """
    Test that issue #11310 has been fixed
    """

    model = convolve_models(models.Sersic2D(), models.Gaussian2D())
    with astropy.conf.set_temp("max_width", 80):
        # fmt: off
        assert str(model) == (
            "Model: CompoundModel\n"
            "Inputs: ('x', 'y')\n"
            "Outputs: ('z',)\n"
            "Model set size: 1\n"
            "Expression: convolve_fft (([0]), ([1]))\n"
            "Components: \n"
            "    [0]: <Sersic2D(amplitude=1., r_eff=1., n=4., "
            "x_0=0., y_0=0., ellip=0., theta=0.)>\n"
            "\n"
            "    [1]: <Gaussian2D(amplitude=1., x_mean=0., y_mean=0., "
            "x_stddev=1., y_stddev=1., theta=0.)>\n"
            "Parameters:\n"
            "    amplitude_0 r_eff_0 n_0 x_0_0 y_0_0 ... y_mean_1 x_stddev_1 y_stddev_1 theta_1\n"
            "    ----------- ------- --- ----- ----- ... -------- ---------- ---------- -------\n"
            "            1.0     1.0 4.0   0.0   0.0 ...      0.0        1.0        1.0     0.0"
        )
        # fmt: on


def test__validate_input_shape():
    model = models.Gaussian1D()
    model._n_models = 2

    _input = np.array(
        [
            [1, 2, 3],
            [4, 5, 6],
        ]
    )

    # Successful validation
    assert model._validate_input_shape(_input, 0, model.inputs, 1, False) == (2, 3)

    # Fail number of axes
    MESSAGE = r"For model_set_axis=2, all inputs must be at least 3-dimensional"
    with pytest.raises(ValueError, match=MESSAGE):
        model._validate_input_shape(_input, 0, model.inputs, 2, True)

    # Fail number of models (has argname)
    MESSAGE = r"Input argument '.*' does not have the correct dimensions in .*"
    with pytest.raises(ValueError, match=MESSAGE):
        model._validate_input_shape(_input, 0, model.inputs, 1, True)

    # Fail number of models  (no argname)
    with pytest.raises(ValueError, match=MESSAGE):
        model._validate_input_shape(_input, 0, [], 1, True)


def test__validate_input_shapes():
    model = models.Gaussian1D()
    model._n_models = 2
    inputs = [mk.MagicMock() for _ in range(3)]
    argnames = mk.MagicMock()
    model_set_axis = mk.MagicMock()
    all_shapes = [mk.MagicMock() for _ in inputs]

    # Successful validation
    with mk.patch.object(
        Model, "_validate_input_shape", autospec=True, side_effect=all_shapes
    ) as mkValidate:
        with mk.patch.object(core, "check_broadcast", autospec=True) as mkCheck:
            assert mkCheck.return_value == model._validate_input_shapes(
                inputs, argnames, model_set_axis
            )
            assert mkCheck.call_args_list == [mk.call(*all_shapes)]
            assert mkValidate.call_args_list == [
                mk.call(model, _input, idx, argnames, model_set_axis, True)
                for idx, _input in enumerate(inputs)
            ]

    # Fail check_broadcast
    MESSAGE = r"All inputs must have identical shapes or must be scalars"
    with mk.patch.object(
        Model, "_validate_input_shape", autospec=True, side_effect=all_shapes
    ) as mkValidate:
        with mk.patch.object(
            core, "check_broadcast", autospec=True, return_value=None
        ) as mkCheck:
            with pytest.raises(ValueError, match=MESSAGE):
                model._validate_input_shapes(inputs, argnames, model_set_axis)
            assert mkCheck.call_args_list == [mk.call(*all_shapes)]
            assert mkValidate.call_args_list == [
                mk.call(model, _input, idx, argnames, model_set_axis, True)
                for idx, _input in enumerate(inputs)
            ]


def test__remove_axes_from_shape():
    model = models.Gaussian1D()

    # len(shape) == 0
    assert model._remove_axes_from_shape((), mk.MagicMock()) == ()

    # axis < 0
    assert model._remove_axes_from_shape((1, 2, 3), -1) == (1, 2)
    assert model._remove_axes_from_shape((1, 2, 3), -2) == (1, 3)
    assert model._remove_axes_from_shape((1, 2, 3), -3) == (2, 3)

    # axis >= len(shape)
    assert model._remove_axes_from_shape((1, 2, 3), 3) == ()
    assert model._remove_axes_from_shape((1, 2, 3), 4) == ()

    # 0 <= axis < len(shape)
    assert model._remove_axes_from_shape((1, 2, 3), 0) == (2, 3)
    assert model._remove_axes_from_shape((1, 2, 3), 1) == (3,)
    assert model._remove_axes_from_shape((1, 2, 3), 2) == ()


def test_get_bounding_box():
    model = models.Const2D(2)

    # No with_bbox
    assert model.get_bounding_box(False) is None

    # No bounding_box
    MESSAGE = r"No bounding box is defined for this model"
    with pytest.raises(NotImplementedError, match=MESSAGE):
        model.bounding_box
    assert model.get_bounding_box(True) is None

    # Normal bounding_box
    model.bounding_box = ((0, 1), (0, 1))
    assert not isinstance(model.bounding_box, CompoundBoundingBox)
    assert model.get_bounding_box(True) == ((0, 1), (0, 1))

    # CompoundBoundingBox with no removal
    bbox = CompoundBoundingBox.validate(
        model,
        {(1,): ((-1, 0), (-1, 0)), (2,): ((0, 1), (0, 1))},
        selector_args=[("y", False)],
    )
    model.bounding_box = bbox
    assert isinstance(model.bounding_box, CompoundBoundingBox)
    # Get using argument not with_bbox
    assert model.get_bounding_box(True) == bbox
    # Get using with_bbox not argument
    assert model.get_bounding_box((1,)) == ((-1, 0), (-1, 0))
    assert model.get_bounding_box((2,)) == ((0, 1), (0, 1))


def test_compound_bounding_box():
    model = models.Gaussian1D()
    truth = models.Gaussian1D()
    bbox1 = CompoundBoundingBox.validate(
        model, {(1,): (-1, 0), (2,): (0, 1)}, selector_args=[("x", False)]
    )
    bbox2 = CompoundBoundingBox.validate(
        model, {(-0.5,): (-1, 0), (0.5,): (0, 1)}, selector_args=[("x", False)]
    )

    # Using with_bounding_box to pass a selector
    model.bounding_box = bbox1
    assert model(-0.5) == truth(-0.5)
    assert model(-0.5, with_bounding_box=(1,)) == truth(-0.5)
    assert np.isnan(model(-0.5, with_bounding_box=(2,)))
    assert model(0.5) == truth(0.5)
    assert model(0.5, with_bounding_box=(2,)) == truth(0.5)
    assert np.isnan(model(0.5, with_bounding_box=(1,)))

    # Using argument value to pass bounding_box
    model.bounding_box = bbox2
    assert model(-0.5) == truth(-0.5)
    assert model(-0.5, with_bounding_box=True) == truth(-0.5)
    assert model(0.5) == truth(0.5)
    assert model(0.5, with_bounding_box=True) == truth(0.5)

    MESSAGE = r"No bounding box is defined for selector: .*"
    with pytest.raises(RuntimeError, match=MESSAGE):
        model(0, with_bounding_box=True)

    model1 = models.Gaussian1D()
    truth1 = models.Gaussian1D()
    model2 = models.Const1D(2)
    truth2 = models.Const1D(2)
    model = model1 + model2
    truth = truth1 + truth2
    assert isinstance(model, CompoundModel)

    model.bounding_box = bbox1
    assert model(-0.5) == truth(-0.5)
    assert model(-0.5, with_bounding_box=1) == truth(-0.5)
    assert np.isnan(model(-0.5, with_bounding_box=(2,)))
    assert model(0.5) == truth(0.5)
    assert model(0.5, with_bounding_box=2) == truth(0.5)
    assert np.isnan(model(0.5, with_bounding_box=(1,)))

    model.bounding_box = bbox2
    assert model(-0.5) == truth(-0.5)
    assert model(-0.5, with_bounding_box=True) == truth(-0.5)
    assert model(0.5) == truth(0.5)
    assert model(0.5, with_bounding_box=True) == truth(0.5)
    with pytest.raises(RuntimeError, match=MESSAGE):
        model(0, with_bounding_box=True)


def test_bind_bounding_box():
    model = models.Polynomial2D(3)
    bbox = ((-1, 1), (-2, 2))

    bind_bounding_box(model, bbox)
    assert model.get_bounding_box() is not None
    assert model.bounding_box == bbox
    assert model.bounding_box["x"] == (-2, 2)
    assert model.bounding_box["y"] == (-1, 1)

    bind_bounding_box(model, bbox, order="F")
    assert model.get_bounding_box() is not None
    assert model.bounding_box == bbox
    assert model.bounding_box["x"] == (-1, 1)
    assert model.bounding_box["y"] == (-2, 2)


def test_bind_compound_bounding_box_using_with_bounding_box_select():
    """
    This demonstrates how to bind multiple bounding_boxes which are
    selectable using the `with_bounding_box`, note there must be a
    fall-back to implicit.
    """
    model = models.Gaussian1D()
    truth = models.Gaussian1D()

    bbox = (0, 1)
    MESSAGE = r"'tuple' object has no attribute 'items"
    with pytest.raises(AttributeError, match=MESSAGE):
        bind_compound_bounding_box(model, bbox, "x")

    bbox = {0: (-1, 0), 1: (0, 1)}
    bind_compound_bounding_box(model, bbox, [("x", False)])

    # No bounding box
    assert model(-0.5) == truth(-0.5)
    assert model(0.5) == truth(0.5)
    assert model(0) == truth(0)
    assert model(1) == truth(1)

    # `with_bounding_box` selects as `-0.5` will not be a key
    assert model(-0.5, with_bounding_box=0) == truth(-0.5)
    assert np.isnan(model(-0.5, with_bounding_box=1))

    # `with_bounding_box` selects as `0.5` will not be a key
    assert model(0.5, with_bounding_box=1) == truth(0.5)
    assert np.isnan(model(0.5, with_bounding_box=(0,)))

    # Fall back onto implicit selector
    assert model(0, with_bounding_box=True) == truth(0)
    assert model(1, with_bounding_box=True) == truth(1)

    # Attempt to fall-back on implicit selector, but no bounding_box
    MESSAGE = r"No bounding box is defined for selector: .*"
    with pytest.raises(RuntimeError, match=MESSAGE):
        model(0.5, with_bounding_box=True)

    # Override implicit selector
    assert np.isnan(model(1, with_bounding_box=0))


def test_fix_inputs_compound_bounding_box():
    base_model = models.Gaussian2D(1, 2, 3, 4, 5)
    bbox = {2.5: (-1, 1), 3.14: (-7, 3)}

    model = fix_inputs(base_model, {"y": 2.5}, bounding_boxes=bbox)
    assert model.bounding_box == (-1, 1)
    model = fix_inputs(base_model, {"x": 2.5}, bounding_boxes=bbox)
    assert model.bounding_box == (-1, 1)

    model = fix_inputs(
        base_model, {"y": 2.5}, bounding_boxes=bbox, selector_args=(("y", True),)
    )
    assert model.bounding_box == (-1, 1)
    model = fix_inputs(
        base_model, {"x": 2.5}, bounding_boxes=bbox, selector_args=(("x", True),)
    )
    assert model.bounding_box == (-1, 1)
    model = fix_inputs(
        base_model, {"x": 2.5}, bounding_boxes=bbox, selector_args=((0, True),)
    )
    assert model.bounding_box == (-1, 1)

    base_model = models.Identity(4)
    bbox = {(2.5, 1.3): ((-1, 1), (-3, 3)), (2.5, 2.71): ((-3, 3), (-1, 1))}

    model = fix_inputs(base_model, {"x0": 2.5, "x1": 1.3}, bounding_boxes=bbox)
    assert model.bounding_box == ((-1, 1), (-3, 3))

    model = fix_inputs(
        base_model,
        {"x0": 2.5, "x1": 1.3},
        bounding_boxes=bbox,
        selector_args=(("x0", True), ("x1", True)),
    )
    assert model.bounding_box == ((-1, 1), (-3, 3))
    model = fix_inputs(
        base_model,
        {"x0": 2.5, "x1": 1.3},
        bounding_boxes=bbox,
        selector_args=((0, True), (1, True)),
    )
    assert model.bounding_box == ((-1, 1), (-3, 3))


def test_model_copy_with_bounding_box():
    model = models.Polynomial2D(2)
    bbox = ModelBoundingBox.validate(model, ((-0.5, 1047.5), (-0.5, 2047.5)), order="F")

    # No bbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert model_copy.get_bounding_box() is None
    assert model.get_bounding_box() is None

    # with bbox
    model.bounding_box = bbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert id(model_copy.bounding_box) != id(model.bounding_box)
    for index, interval in model.bounding_box.intervals.items():
        interval_copy = model_copy.bounding_box.intervals[index]
        assert interval == interval_copy
        assert id(interval) != interval_copy

    # add model to compound model
    model1 = model | models.Identity(1)
    model_copy = model1.copy()
    assert id(model_copy) != id(model1)
    assert model_copy.get_bounding_box() is None
    assert model1.get_bounding_box() is None


def test_compound_model_copy_with_bounding_box():
    model = models.Shift(1) & models.Shift(2) & models.Identity(1)
    model.inputs = ("x", "y", "slit_id")
    bbox = ModelBoundingBox.validate(
        model, ((-0.5, 1047.5), (-0.5, 2047.5), (-np.inf, np.inf)), order="F"
    )

    # No bbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert model_copy.get_bounding_box() is None
    assert model.get_bounding_box() is None

    # with bbox
    model.bounding_box = bbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert id(model_copy.bounding_box) != id(model.bounding_box)
    for index, interval in model.bounding_box.intervals.items():
        interval_copy = model_copy.bounding_box.intervals[index]
        assert interval == interval_copy
        assert id(interval) != interval_copy

    # add model to compound model
    model1 = model | models.Identity(3)
    model_copy = model1.copy()
    assert id(model_copy) != id(model1)
    assert model_copy.get_bounding_box() is None
    assert model1.get_bounding_box() is None


def test_model_copy_with_compound_bounding_box():
    model = models.Polynomial2D(2)
    bbox = {(0,): (-0.5, 1047.5), (1,): (-0.5, 3047.5)}
    cbbox = CompoundBoundingBox.validate(
        model, bbox, selector_args=[("x", True)], order="F"
    )

    # No cbbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert model_copy.get_bounding_box() is None
    assert model.get_bounding_box() is None

    # with cbbox
    model.bounding_box = cbbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert id(model_copy.bounding_box) != id(model.bounding_box)
    assert model_copy.bounding_box.selector_args == model.bounding_box.selector_args
    assert id(model_copy.bounding_box.selector_args) != id(
        model.bounding_box.selector_args
    )
    for selector, bbox in model.bounding_box.bounding_boxes.items():
        for index, interval in bbox.intervals.items():
            interval_copy = model_copy.bounding_box.bounding_boxes[selector].intervals[
                index
            ]
            assert interval == interval_copy
            assert id(interval) != interval_copy

    # add model to compound model
    model1 = model | models.Identity(1)
    model_copy = model1.copy()
    assert id(model_copy) != id(model1)
    assert model_copy.get_bounding_box() is None
    assert model1.get_bounding_box() is None


def test_compound_model_copy_with_compound_bounding_box():
    model = models.Shift(1) & models.Shift(2) & models.Identity(1)
    model.inputs = ("x", "y", "slit_id")
    bbox = {
        (0,): ((-0.5, 1047.5), (-0.5, 2047.5)),
        (1,): ((-0.5, 3047.5), (-0.5, 4047.5)),
    }
    cbbox = CompoundBoundingBox.validate(
        model, bbox, selector_args=[("slit_id", True)], order="F"
    )

    # No cbbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert model_copy.get_bounding_box() is None
    assert model.get_bounding_box() is None

    # with cbbox
    model.bounding_box = cbbox
    model_copy = model.copy()
    assert id(model_copy) != id(model)
    assert id(model_copy.bounding_box) != id(model.bounding_box)
    assert model_copy.bounding_box.selector_args == model.bounding_box.selector_args
    assert id(model_copy.bounding_box.selector_args) != id(
        model.bounding_box.selector_args
    )
    for selector, bbox in model.bounding_box.bounding_boxes.items():
        for index, interval in bbox.intervals.items():
            interval_copy = model_copy.bounding_box.bounding_boxes[selector].intervals[
                index
            ]
            assert interval == interval_copy
            assert id(interval) != interval_copy

    # add model to compound model
    model1 = model | models.Identity(3)
    model_copy = model1.copy()
    assert id(model_copy) != id(model1)
    assert model_copy.get_bounding_box() is None
    assert model1.get_bounding_box() is None


def test_compound_model_copy_user_attribute():
    """Regression test for issue #12370"""

    model = models.Gaussian2D(100, 25, 25, 5, 5) | models.Identity(1)
    model.xname = "x_mean"  # user-defined attribute
    assert hasattr(model, "xname")
    assert model.xname == "x_mean"

    model_copy = model.copy()
    model_copy.xname
    assert hasattr(model_copy, "xname")
    assert model_copy.xname == "x_mean"


def test_model_mixed_array_scalar_bounding_box():
    """Regression test for issue #12319"""

    model = models.Gaussian2D()
    bbox = ModelBoundingBox.validate(model, ((-1, 1), (-np.inf, np.inf)), order="F")
    model.bounding_box = bbox

    x = np.array([-0.5, 0.5])
    y = 0

    # Everything works when its all in the bounding box
    assert (model(x, y) == (model(x, y, with_bounding_box=True))).all()


def test_compound_model_mixed_array_scalar_bounding_box():
    """Regression test for issue #12319"""

    model = models.Shift(1) & models.Shift(2) & models.Identity(1)
    model.inputs = ("x", "y", "slit_id")
    bbox = ModelBoundingBox.validate(
        model, ((-0.5, 1047.5), (-0.5, 2047.5), (-np.inf, np.inf)), order="F"
    )
    model.bounding_box = bbox
    x = np.array([1000, 1001])
    y = np.array([2000, 2001])
    slit_id = 0

    # Everything works when its all in the bounding box
    value0 = model(x, y, slit_id)
    value1 = model(x, y, slit_id, with_bounding_box=True)
    assert_equal(value0, value1)


def test_model_with_bounding_box_true_and_single_output():
    """Regression test for issue #12373"""

    model = models.Mapping((1,))
    x = [1, 2]
    y = [3, 4]

    # Check baseline
    assert_equal(model(x, y), [3, 4])
    # Check with_bounding_box=True should be the same
    assert_equal(model(x, y, with_bounding_box=True), [3, 4])

    model.bounding_box = ((-np.inf, np.inf), (-np.inf, np.inf))
    # Check baseline
    assert_equal(model(x, y), [3, 4])
    # Check with_bounding_box=True should be the same
    assert_equal(model(x, y, with_bounding_box=True), [3, 4])


def test_compound_model_with_bounding_box_true_and_single_output():
    """Regression test for issue #12373"""

    model = models.Mapping((1,)) | models.Shift(1)
    x = [1, 2]
    y = [3, 4]

    # Check baseline
    assert_equal(model(x, y), [4, 5])
    # Check with_bounding_box=True should be the same
    assert_equal(model(x, y, with_bounding_box=True), [4, 5])

    model.bounding_box = ((-np.inf, np.inf), (-np.inf, np.inf))
    # Check baseline
    assert_equal(model(x, y), [4, 5])
    # Check with_bounding_box=True should be the same
    assert_equal(model(x, y, with_bounding_box=True), [4, 5])


def test_bounding_box_pass_with_ignored():
    """Test the possibility of setting ignored variables in bounding box"""

    model = models.Polynomial2D(2)
    bbox = ModelBoundingBox.validate(model, (-1, 1), ignored=["y"])
    model.bounding_box = bbox
    assert model.bounding_box.bounding_box() == (-1, 1)
    assert model.bounding_box == bbox

    model = models.Polynomial2D(2)
    bind_bounding_box(model, (-1, 1), ignored=["y"])
    assert model.bounding_box.bounding_box() == (-1, 1)
    assert model.bounding_box == bbox


def test_compound_bounding_box_pass_with_ignored():
    model = models.Shift(1) & models.Shift(2) & models.Identity(1)
    model.inputs = ("x", "y", "slit_id")
    bbox = {
        (0,): (-0.5, 1047.5),
        (1,): (-0.5, 2047.5),
    }
    cbbox = CompoundBoundingBox.validate(
        model, bbox, selector_args=[("slit_id", True)], ignored=["y"], order="F"
    )
    model.bounding_box = cbbox

    model = models.Shift(1) & models.Shift(2) & models.Identity(1)
    model.inputs = ("x", "y", "slit_id")
    bind_compound_bounding_box(
        model, bbox, selector_args=[("slit_id", True)], ignored=["y"], order="F"
    )
    assert model.bounding_box == cbbox


@pytest.mark.parametrize("int_type", [int, np.int32, np.int64, np.uint32, np.uint64])
def test_model_integer_indexing(int_type):
    """Regression for PR 12561; verify that compound model components
    can be accessed by integer index"""
    gauss = models.Gaussian2D()
    airy = models.AiryDisk2D()
    compound = gauss + airy

    assert compound[int_type(0)] == gauss
    assert compound[int_type(1)] == airy


def test_model_string_indexing():
    """Regression for PR 12561; verify that compound model components
    can be accessed by indexing with model name"""
    gauss = models.Gaussian2D()
    gauss.name = "Model1"
    airy = models.AiryDisk2D()
    airy.name = "Model2"
    compound = gauss + airy

    assert compound["Model1"] == gauss
    assert compound["Model2"] == airy
