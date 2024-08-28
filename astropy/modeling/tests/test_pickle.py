"""Tests that models are picklable."""

from pickle import dumps, loads

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.modeling import (
    functional_models,
    mappings,
    math_functions,
    physical_models,
    polynomial,
    powerlaws,
    projections,
    rotations,
    spline,
    tabular,
)
from astropy.modeling.math_functions import ArctanhUfunc
from astropy.utils.compat.optional_deps import HAS_SCIPY

MATH_FUNCTIONS = (func for func in math_functions.__all__ if func != "ArctanhUfunc")


PROJ_TO_REMOVE = (
    [
        "Projection",
        "Pix2SkyProjection",
        "Sky2PixProjection",
        "Zenithal",
        "Conic",
        "Cylindrical",
        "PseudoCylindrical",
        "PseudoConic",
        "QuadCube",
        "HEALPix",
        "AffineTransformation2D",
        "projcodes",
        "Pix2Sky_ZenithalPerspective",
    ]
    + [f"Pix2Sky_{code}" for code in projections.projcodes]
    + [f"Sky2Pix_{code}" for code in projections.projcodes]
)

PROJECTIONS = (func for func in projections.__all__ if func not in PROJ_TO_REMOVE)

OTHER_MODELS = [
    mappings.Mapping((1, 0)),
    mappings.Identity(2),
    ArctanhUfunc(),
    rotations.Rotation2D(23),
    tabular.Tabular1D(lookup_table=[1, 2, 3, 4]),
    tabular.Tabular2D(lookup_table=[[1, 2, 3, 4], [5, 6, 7, 8]]),
]

POLYNOMIALS_1D = ["Chebyshev1D", "Hermite1D", "Legendre1D", "Polynomial1D"]
POLYNOMIALS_2D = ["Chebyshev2D", "Hermite2D", "Legendre2D", "InverseSIP"]

ROTATIONS = [
    rotations.RotateCelestial2Native(12, 23, 34),
    rotations.RotateNative2Celestial(12, 23, 34),
    rotations.EulerAngleRotation(12, 23, 34, "xyz"),
    rotations.RotationSequence3D([12, 23, 34], axes_order="xyz"),
    rotations.SphericalRotationSequence([12, 23, 34], "xyz"),
    rotations.Rotation2D(12),
]


@pytest.fixture
def inputs():
    return 0.3, 0.4


@pytest.fixture
def inputs_math():
    return 1, -0.5


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("model", functional_models.__all__)
def test_pickle_functional(inputs, model):
    m = getattr(functional_models, model)()
    mp = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(inputs[0]), mp(inputs[0]))
    else:
        assert_allclose(m(*inputs), mp(*inputs))


@pytest.mark.parametrize("model", MATH_FUNCTIONS)
def test_pickle_math_functions(inputs_math, model):
    m = getattr(math_functions, model)()
    mp = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(inputs_math[0]), mp(inputs_math[0]))
    else:
        assert_allclose(m(*inputs_math), mp(*inputs_math))


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("m", OTHER_MODELS)
def test_pickle_other(inputs, m):
    mp = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(inputs[0]), mp(inputs[0]))
    else:
        assert_allclose(m(*inputs), mp(*inputs))


def test_pickle_units_mapping(inputs):
    m = mappings.UnitsMapping(((u.m, None),))
    mp = loads(dumps(m))
    assert_allclose(m(inputs[0] * u.km), mp(inputs[0] * u.km))


def test_pickle_affine_transformation_2D(inputs):
    m = projections.AffineTransformation2D(matrix=[[1, 1], [1, 1]], translation=[1, 1])
    m.matrix.fixed = True
    mp = loads(dumps(m))
    assert_allclose(m(*inputs), mp(*inputs))
    assert m.matrix.fixed is True


@pytest.mark.parametrize("model", physical_models.__all__)
def test_pickle_physical_models(inputs, model):
    m = getattr(physical_models, model)()
    m1 = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(inputs[0]), m1(inputs[0]))
    else:
        assert_allclose(m(*inputs), m1(*inputs))


@pytest.mark.parametrize("model", POLYNOMIALS_1D)
def test_pickle_1D_polynomials(inputs, model):
    m = getattr(polynomial, model)
    m = m(2)
    m1 = loads(dumps(m))
    assert_allclose(m(inputs[1]), m1(inputs[0]))


@pytest.mark.parametrize("model", POLYNOMIALS_2D)
def test_pickle_2D_polynomials(inputs, model):
    m = getattr(polynomial, model)
    m = m(2, 3)
    m1 = loads(dumps(m))
    assert_allclose(m(*inputs), m1(*inputs))


def test_pickle_polynomial_2D(inputs):
    # Polynomial2D is initialized with 1 degree but
    # requires 2 inputs
    m = polynomial.Polynomial2D
    m = m(2)
    m1 = loads(dumps(m))
    assert_allclose(m(*inputs), m1(*inputs))


def test_pickle_sip(inputs):
    m = polynomial.SIP
    m = m((21, 23), 2, 3)
    m1 = loads(dumps(m))
    assert_allclose(m(*inputs), m1(*inputs))


@pytest.mark.parametrize("model", powerlaws.__all__)
def test_pickle_powerlaws(inputs, model):
    m = getattr(powerlaws, model)()
    m1 = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(inputs[0]), m1(inputs[0]))
    else:
        assert_allclose(m(*inputs), m1(*inputs))


@pytest.mark.parametrize("model", PROJECTIONS)
def test_pickle_projections(inputs, model):
    m = getattr(projections, model)()
    m1 = loads(dumps(m))
    assert_allclose(m(*inputs), m1(*inputs))


@pytest.mark.parametrize("m", ROTATIONS)
def test_pickle_rotations(inputs, m):
    mp = loads(dumps(m))

    if m.n_inputs == 2:
        assert_allclose(m(*inputs), mp(*inputs))
    else:
        assert_allclose(m(inputs[0], *inputs), mp(inputs[0], *inputs))


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_pickle_spline(inputs):
    def func(x, noise):
        return np.exp(-(x**2)) + 0.1 * noise

    noise = np.random.randn(50)
    x = np.linspace(-3, 3, 50)
    y = func(x, noise)

    fitter = spline.SplineInterpolateFitter()
    spl = spline.Spline1D(degree=3)
    m = fitter(spl, x, y)

    mp = loads(dumps(m))
    assert_allclose(m(inputs[0]), mp(inputs[0]))
