"""Tests that models are picklable."""

from pickle import dumps, loads
from numpy.testing import assert_allclose
import pytest

from astropy import units as u
from astropy.utils.compat.optional_deps import HAS_SCIPY

from astropy.modeling import functional_models
from astropy.modeling import mappings
from astropy.modeling import math_functions
from astropy.modeling import physical_models
from astropy.modeling import polynomial
from astropy.modeling import powerlaws
from astropy.modeling import projections
from astropy.modeling import rotations
from astropy.modeling import spline
from astropy.modeling import tabular

from astropy.modeling.math_functions import ArctanhUfunc

math_functions_all = math_functions.__all__[:]
math_functions_all.remove("ArctanhUfunc")

polynomial_all = polynomial.__all__[:]
# remove base classes
polynomial_all.remove("PolynomialModel")
polynomial_all.remove("OrthoPolynomialBase")

projections_all = projections.__all__[:]
proj_to_remove = [
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
for p in proj_to_remove:
    projections_all.remove(p)

for code in projections.projcodes:
    projections_all.remove(f"Pix2Sky_{code}")
    projections_all.remove(f"Sky2Pix_{code}")


# can parametrize on x,y as well
x, y = 0.3, 0.4
x_math, y_math = 1, -0.5


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("model", functional_models.__all__[:])
def test_pickle_functional(model):
    m = getattr(functional_models, model)()
    m1 = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(x), m1(x))
    else:
        assert_allclose(m(x, y), m1(x, y))


@pytest.mark.parametrize("model", math_functions_all)
def test_pickle_math_functions(model):
    m = getattr(math_functions, model)()
    m1 = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(x_math), m1(x_math))
    else:
        assert_allclose(m(x_math, y_math), m1(x_math, y_math))


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_pickle_other():
    # Test models which don't fit in the other test functions.

    x = 0.5
    y = 1

    model = mappings.Mapping((1, 0))
    modelp = loads(dumps(model))
    assert_allclose(model(x, y), modelp(x, y))

    model = mappings.Identity(2)
    modelp = loads(dumps(model))
    assert_allclose(model(x, y), modelp(x, y))

    model = mappings.UnitsMapping(((u.m, None),))
    modelp = loads(dumps(model))
    assert_allclose(model(x * u.km), modelp(x * u.km))

    model = ArctanhUfunc()
    modelp = loads(dumps(model))
    assert_allclose(model(x), modelp(x))

    model = rotations.Rotation2D(23)
    modelp = loads(dumps(model))
    assert_allclose(model(x, y), modelp(x, y))

    model = tabular.Tabular1D(lookup_table=[1, 2, 3, 4])
    modelp = loads(dumps(model))
    assert_allclose(model(y), modelp(y))

    model = tabular.Tabular2D(lookup_table=[[1, 2, 3, 4], [5, 6, 7, 8]])
    modelp = loads(dumps(model))
    assert_allclose(model(x, y), modelp(x, y))

    model = spline.Spline1D()
    modelp = loads(dumps(model))
    assert_allclose(model(x), modelp(x))

    model = projections.AffineTransformation2D(
        matrix=[[1, 1], [1, 1]], translation=[1, 1]
    )
    model.matrix.fixed = True
    modelp = loads(dumps(model))
    assert_allclose(model(x, y), modelp(x, y))
    assert model.matrix.fixed is True


@pytest.mark.parametrize("model", physical_models.__all__)
def test_pickle_physical_models(model):
    m = getattr(physical_models, model)()
    m1 = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(x), m1(x))
    else:
        assert_allclose(m(x, y), m1(x, y))


def test_pickle_polynomial():
    # models initialized with 1 degree
    models2d = ["Chebyshev2D", "Hermite2D", "Legendre2D", "InverseSIP"]
    # models initialized with 2 degree
    models1d = ["Chebyshev1D", "Hermite1D", "Legendre1D", "Polynomial1D"]

    sip = ["SIP", "InverseSIP"]

    for model in models1d:
        m = getattr(polynomial, model)
        m = m(2)
        m1 = loads(dumps(m))
        assert_allclose(m(x), m1(x))

    for model in models2d:
        m = getattr(polynomial, model)
        m = m(2, 3)
        m1 = loads(dumps(m))
        assert_allclose(m(x, y), m1(x, y))

    # Polynomial2D is initialized with 1 degree but
    # requires 2 inputs
    m = polynomial.Polynomial2D
    m = m(2)
    m1 = loads(dumps(m))
    assert_allclose(m(x, y), m1(x, y))

    m = polynomial.SIP
    m = m((21, 23), 2, 3)
    m1 = loads(dumps(m))
    assert_allclose(m(x, y), m1(x, y))


@pytest.mark.parametrize("model", powerlaws.__all__)
def test_pickle_powerlaws(model):
    m = getattr(powerlaws, model)()
    m1 = loads(dumps(m))
    if m.n_inputs == 1:
        assert_allclose(m(x), m1(x))
    else:
        assert_allclose(m(x, y), m1(x, y))


@pytest.mark.parametrize("model", projections_all)
def test_pickle_projections(model):
    m = getattr(projections, model)()
    m1 = loads(dumps(m))
    assert_allclose(m(x, y), m1(x, y))


def test_pickle_rotations():
    for model in ["RotateCelestial2Native", "RotateNative2Celestial"]:
        m = getattr(rotations, model)(12, 23, 34)
        m1 = loads(dumps(m))
        assert_allclose(m(x, y), m1(x, y))

    m = rotations.EulerAngleRotation(12, 23, 34, "xyz")
    m1 = loads(dumps(m))
    assert_allclose(m(x, y), m1(x, y))

    m = rotations.RotationSequence3D([12, 23, 34], axes_order="xyz")
    m1 = loads(dumps(m))
    assert_allclose(m(x, y, y), m1(x, y, y))

    m = rotations.SphericalRotationSequence([12, 23, 34], "xyz")
    m1 = loads(dumps(m))
    assert_allclose(m(x, y), m1(x, y))

    m = rotations.Rotation2D(12)
    m1 = loads(dumps(m))
    assert_allclose(m(x, y), m1(x, y))
