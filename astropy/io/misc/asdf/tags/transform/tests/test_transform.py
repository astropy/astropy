# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

asdf = pytest.importorskip("asdf")

import warnings

import asdf
import numpy as np
from asdf import AsdfFile, util
from asdf.exceptions import AsdfDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        category=AsdfDeprecationWarning,
        message=r"asdf.tests.helpers is deprecated.*",
    )
    from asdf.tests.helpers import assert_roundtrip_tree, yaml_to_asdf

from packaging.version import Version

import astropy.units as u
from astropy.modeling import models as astmodels
from astropy.modeling.core import fix_inputs
from astropy.utils.compat.optional_deps import HAS_SCIPY


def custom_and_analytical_inverse():
    p1 = astmodels.Polynomial1D(1)
    p2 = astmodels.Polynomial1D(1)
    p3 = astmodels.Polynomial1D(1)
    p4 = astmodels.Polynomial1D(1)
    m1 = p1 & p2
    m2 = p3 & p4
    m1.inverse = m2
    return m1


def custom_inputs_outputs():
    m = astmodels.Gaussian2D()
    m.inputs = ("a", "b")
    m.outputs = ("c",)
    return m


test_models = [
    astmodels.Identity(2),
    astmodels.Polynomial1D(2, c0=1, c1=2, c2=3),
    astmodels.Polynomial2D(1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Shift(2.0),
    astmodels.Hermite1D(2, c0=2, c1=3, c2=0.5),
    astmodels.Legendre1D(2, c0=2, c1=3, c2=0.5),
    astmodels.Chebyshev1D(2, c0=2, c1=3, c2=0.5),
    astmodels.Chebyshev2D(1, 1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Legendre2D(1, 1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Hermite2D(1, 1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Scale(3.4),
    astmodels.RotateNative2Celestial(5.63, -72.5, 180),
    astmodels.Multiply(3),
    astmodels.Multiply(10 * u.m),
    astmodels.RotateCelestial2Native(5.63, -72.5, 180),
    astmodels.EulerAngleRotation(23, 14, 2.3, axes_order="xzx"),
    astmodels.Mapping((0, 1), n_inputs=3),
    astmodels.Shift(2.0 * u.deg),
    astmodels.Scale(3.4 * u.deg),
    astmodels.RotateNative2Celestial(5.63 * u.deg, -72.5 * u.deg, 180 * u.deg),
    astmodels.RotateCelestial2Native(5.63 * u.deg, -72.5 * u.deg, 180 * u.deg),
    astmodels.RotationSequence3D([1.2, 2.3, 3.4, 0.3], "xyzx"),
    astmodels.SphericalRotationSequence([1.2, 2.3, 3.4, 0.3], "xyzy"),
    astmodels.AiryDisk2D(amplitude=10.0, x_0=0.5, y_0=1.5),
    astmodels.Box1D(amplitude=10.0, x_0=0.5, width=5.0),
    astmodels.Box2D(amplitude=10.0, x_0=0.5, x_width=5.0, y_0=1.5, y_width=7.0),
    astmodels.Const1D(amplitude=5.0),
    astmodels.Const2D(amplitude=5.0),
    astmodels.Disk2D(amplitude=10.0, x_0=0.5, y_0=1.5, R_0=5.0),
    astmodels.Ellipse2D(amplitude=10.0, x_0=0.5, y_0=1.5, a=2.0, b=4.0, theta=0.1),
    astmodels.Exponential1D(amplitude=10.0, tau=3.5),
    astmodels.Gaussian1D(amplitude=10.0, mean=5.0, stddev=3.0),
    astmodels.Gaussian2D(
        amplitude=10.0, x_mean=5.0, y_mean=5.0, x_stddev=3.0, y_stddev=3.0
    ),
    astmodels.KingProjectedAnalytic1D(amplitude=10.0, r_core=5.0, r_tide=2.0),
    astmodels.Logarithmic1D(amplitude=10.0, tau=3.5),
    astmodels.Lorentz1D(amplitude=10.0, x_0=0.5, fwhm=2.5),
    astmodels.Moffat1D(amplitude=10.0, x_0=0.5, gamma=1.2, alpha=2.5),
    astmodels.Moffat2D(amplitude=10.0, x_0=0.5, y_0=1.5, gamma=1.2, alpha=2.5),
    astmodels.Planar2D(slope_x=0.5, slope_y=1.2, intercept=2.5),
    astmodels.RedshiftScaleFactor(z=2.5),
    astmodels.RickerWavelet1D(amplitude=10.0, x_0=0.5, sigma=1.2),
    astmodels.RickerWavelet2D(amplitude=10.0, x_0=0.5, y_0=1.5, sigma=1.2),
    astmodels.Ring2D(amplitude=10.0, x_0=0.5, y_0=1.5, r_in=5.0, width=10.0),
    astmodels.Sersic1D(amplitude=10.0, r_eff=1.0, n=4.0),
    astmodels.Sersic2D(
        amplitude=10.0, r_eff=1.0, n=4.0, x_0=0.5, y_0=1.5, ellip=0.0, theta=0.0
    ),
    astmodels.Sine1D(amplitude=10.0, frequency=0.5, phase=1.0),
    astmodels.Cosine1D(amplitude=10.0, frequency=0.5, phase=1.0),
    astmodels.Tangent1D(amplitude=10.0, frequency=0.5, phase=1.0),
    astmodels.ArcSine1D(amplitude=10.0, frequency=0.5, phase=1.0),
    astmodels.ArcCosine1D(amplitude=10.0, frequency=0.5, phase=1.0),
    astmodels.ArcTangent1D(amplitude=10.0, frequency=0.5, phase=1.0),
    astmodels.Trapezoid1D(amplitude=10.0, x_0=0.5, width=5.0, slope=1.0),
    astmodels.TrapezoidDisk2D(amplitude=10.0, x_0=0.5, y_0=1.5, R_0=5.0, slope=1.0),
    astmodels.Voigt1D(x_0=0.55, amplitude_L=10.0, fwhm_L=0.5, fwhm_G=0.9),
    astmodels.BlackBody(scale=10.0, temperature=6000.0 * u.K),
    astmodels.Drude1D(amplitude=10.0, x_0=0.5, fwhm=2.5),
    astmodels.Plummer1D(mass=10.0, r_plum=5.0),
    astmodels.BrokenPowerLaw1D(amplitude=10, x_break=0.5, alpha_1=2.0, alpha_2=3.5),
    astmodels.ExponentialCutoffPowerLaw1D(10, 0.5, 2.0, 7.0),
    astmodels.LogParabola1D(
        amplitude=10,
        x_0=0.5,
        alpha=2.0,
        beta=3.0,
    ),
    astmodels.PowerLaw1D(amplitude=10.0, x_0=0.5, alpha=2.0),
    astmodels.SmoothlyBrokenPowerLaw1D(
        amplitude=10.0, x_break=5.0, alpha_1=2.0, alpha_2=3.0, delta=0.5
    ),
    custom_and_analytical_inverse(),
    custom_inputs_outputs(),
]

if HAS_SCIPY:
    test_models.append(
        astmodels.Spline1D(
            np.array([-3.0, -3.0, -3.0, -3.0, -1.0, 0.0, 1.0, 3.0, 3.0, 3.0, 3.0]),
            np.array(
                [
                    0.10412331,
                    0.07013616,
                    -0.18799552,
                    1.35953147,
                    -0.15282581,
                    0.03923,
                    -0.04297299,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]
            ),
            3,
        )
    )

math_models = []

for kl in astmodels.math.__all__:
    klass = getattr(astmodels.math, kl)
    math_models.append(klass())

test_models.extend(math_models)

test_models_with_constraints = [
    astmodels.Legendre2D(
        x_degree=1,
        y_degree=1,
        c0_0=1,
        c0_1=2,
        c1_0=3,
        fixed={"c1_0": True, "c0_1": True},
        bounds={"c0_0": (-10, 10)},
    )
]
test_models.extend(test_models_with_constraints)


def test_transforms_compound(tmpdir):
    tree = {
        "compound": astmodels.Shift(1) & astmodels.Shift(2)
        | astmodels.Sky2Pix_TAN()
        | astmodels.Rotation2D()
        | astmodels.AffineTransformation2D([[2, 0], [0, 2]], [42, 32])
        + astmodels.Rotation2D(32)
    }

    assert_roundtrip_tree(tree, tmpdir)


def test_inverse_transforms(tmpdir):
    rotation = astmodels.Rotation2D(32)
    rotation.inverse = astmodels.Rotation2D(45)

    real_rotation = astmodels.Rotation2D(32)

    tree = {"rotation": rotation, "real_rotation": real_rotation}

    def check(ff):
        assert ff.tree["rotation"].inverse.angle == 45

    assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check)


@pytest.mark.parametrize("model", test_models)
def test_single_model(tmpdir, model):
    with warnings.catch_warnings():
        # Some schema files are missing from asdf<=2.6.0 which causes warnings
        if Version(asdf.__version__) <= Version("2.6.0"):
            warnings.filterwarnings("ignore", "Unable to locate schema file")
        tree = {"single_model": model}
        assert_roundtrip_tree(tree, tmpdir)


def test_name(tmpdir):
    def check(ff):
        assert ff.tree["rot"].name == "foo"

    tree = {"rot": astmodels.Rotation2D(23, name="foo")}
    assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check)


def test_zenithal_with_arguments(tmpdir):
    tree = {"azp": astmodels.Sky2Pix_AZP(0.5, 0.3)}

    assert_roundtrip_tree(tree, tmpdir)


def test_naming_of_compound_model(tmpdir):
    """Issue #87"""

    def asdf_check(ff):
        assert ff.tree["model"].name == "compound_model"

    offx = astmodels.Shift(1)
    scl = astmodels.Scale(2)
    model = (offx | scl).rename("compound_model")
    tree = {"model": model}
    assert_roundtrip_tree(tree, tmpdir, asdf_check_func=asdf_check)


@pytest.mark.slow
def test_generic_projections(tmpdir):
    from astropy.io.misc.asdf.tags.transform import projections

    for tag_name, (name, params, version) in projections._generic_projections.items():
        tree = {
            "forward": util.resolve_name(
                f"astropy.modeling.projections.Sky2Pix_{name}"
            )(),
            "backward": util.resolve_name(
                f"astropy.modeling.projections.Pix2Sky_{name}"
            )(),
        }
        with warnings.catch_warnings():
            # Some schema files are missing from asdf<=2.4.2 which causes warnings
            if Version(asdf.__version__) <= Version("2.5.1"):
                warnings.filterwarnings("ignore", "Unable to locate schema file")
            assert_roundtrip_tree(tree, tmpdir)


def test_tabular_model(tmpdir):
    points = np.arange(0, 5)
    values = [1.0, 10, 2, 45, -3]
    model = astmodels.Tabular1D(points=points, lookup_table=values)
    tree = {"model": model}
    assert_roundtrip_tree(tree, tmpdir)
    table = np.array([[3.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]])
    points = ([1, 2, 3], [1, 2, 3])
    model2 = astmodels.Tabular2D(
        points,
        lookup_table=table,
        bounds_error=False,
        fill_value=None,
        method="nearest",
    )
    tree = {"model": model2}
    assert_roundtrip_tree(tree, tmpdir)


def test_bounding_box(tmpdir):
    model = astmodels.Shift(1) & astmodels.Shift(2)
    model.bounding_box = ((1, 3), (2, 4))
    tree = {"model": model}
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.parametrize("standard_version", asdf.versioning.supported_versions)
def test_const1d(tmpdir, standard_version):
    assert_roundtrip_tree(
        {"model": astmodels.Const1D(amplitude=5.0)},
        tmpdir,
        init_options={"version": standard_version},
    )


@pytest.mark.parametrize("standard_version", asdf.versioning.supported_versions)
@pytest.mark.parametrize(
    "model",
    [
        astmodels.Polynomial1D(1, c0=5, c1=17),
        astmodels.Polynomial1D(1, c0=5, c1=17, domain=[-5, 4], window=[-2, 3]),
        astmodels.Polynomial2D(2, c0_0=3, c1_0=5, c0_1=7),
        astmodels.Polynomial2D(
            2,
            c0_0=3,
            c1_0=5,
            c0_1=7,
            x_domain=[-2, 2],
            y_domain=[-4, 4],
            x_window=[-6, 6],
            y_window=[-8, 8],
        ),
    ],
)
def test_polynomial(tmpdir, standard_version, model):
    assert_roundtrip_tree(
        {"model": model}, tmpdir, init_options={"version": standard_version}
    )


def test_domain_orthopoly(tmpdir):
    model1d = astmodels.Chebyshev1D(2, c0=2, c1=3, c2=0.5, domain=[-2, 2])
    model2d = astmodels.Chebyshev2D(
        1, 1, c0_0=1, c0_1=2, c1_0=3, x_domain=[-2, 2], y_domain=[-2, 2]
    )
    fa = AsdfFile()
    fa.tree["model1d"] = model1d
    fa.tree["model2d"] = model2d

    file_path = str(tmpdir.join("orthopoly_domain.asdf"))
    fa.write_to(file_path)
    with asdf.open(file_path) as f:
        assert f.tree["model1d"](1.8) == model1d(1.8)
        assert f.tree["model2d"](1.8, -1.5) == model2d(1.8, -1.5)


def test_window_orthopoly(tmpdir):
    model1d = astmodels.Chebyshev1D(
        2, c0=2, c1=3, c2=0.5, domain=[-2, 2], window=[-0.5, 0.5]
    )
    model2d = astmodels.Chebyshev2D(
        1,
        1,
        c0_0=1,
        c0_1=2,
        c1_0=3,
        x_domain=[-2, 2],
        y_domain=[-2, 2],
        x_window=[-0.5, 0.5],
        y_window=[-0.1, 0.5],
    )
    fa = AsdfFile()
    fa.tree["model1d"] = model1d
    fa.tree["model2d"] = model2d

    file_path = str(tmpdir.join("orthopoly_window.asdf"))
    fa.write_to(file_path)
    with asdf.open(file_path) as f:
        assert f.tree["model1d"](1.8) == model1d(1.8)
        assert f.tree["model2d"](1.8, -1.5) == model2d(1.8, -1.5)


def test_linear1d(tmpdir):
    model = astmodels.Linear1D()
    tree = {"model": model}
    assert_roundtrip_tree(tree, tmpdir)


def test_linear1d_quantity(tmpdir):
    model = astmodels.Linear1D(1 * u.nm, 1 * (u.nm / u.pixel))
    tree = {"model": model}
    assert_roundtrip_tree(tree, tmpdir)


def test_tabular_model_units(tmpdir):
    points = np.arange(0, 5) * u.pix
    values = [1.0, 10, 2, 45, -3] * u.nm
    model = astmodels.Tabular1D(points=points, lookup_table=values)
    tree = {"model": model}
    assert_roundtrip_tree(tree, tmpdir)
    table = np.array([[3.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]]) * u.nm
    points = ([1, 2, 3], [1, 2, 3]) * u.pix
    model2 = astmodels.Tabular2D(
        points,
        lookup_table=table,
        bounds_error=False,
        fill_value=None,
        method="nearest",
    )
    tree = {"model": model2}
    assert_roundtrip_tree(tree, tmpdir)


def test_fix_inputs(tmpdir):
    with warnings.catch_warnings():
        # Some schema files are missing from asdf<=2.4.2 which causes warnings
        if Version(asdf.__version__) <= Version("2.5.1"):
            warnings.filterwarnings("ignore", "Unable to locate schema file")

        model0 = astmodels.Pix2Sky_TAN()
        model0.input_units_equivalencies = {
            "x": u.dimensionless_angles(),
            "y": u.dimensionless_angles(),
        }
        model1 = astmodels.Rotation2D()
        model = model0 | model1

        tree = {
            "compound": fix_inputs(model, {"x": 45}),
            "compound1": fix_inputs(model, {0: 45}),
        }

        assert_roundtrip_tree(tree, tmpdir)


def test_fix_inputs_type(tmpdir):
    with pytest.raises(TypeError):
        tree = {"compound": fix_inputs(3, {"x": 45})}
        assert_roundtrip_tree(tree, tmpdir)

    with pytest.raises(AttributeError):
        tree = {"compound": astmodels.Pix2Sky_TAN() & {"x": 45}}
        assert_roundtrip_tree(tree, tmpdir)


comp_model = custom_and_analytical_inverse()


@pytest.mark.parametrize(
    "model",
    [
        astmodels.Shift(1) & astmodels.Shift(2) | comp_model,
        comp_model | astmodels.Shift(1) & astmodels.Shift(2),
        astmodels.Shift(1) & comp_model,
        comp_model & astmodels.Shift(1),
    ],
)
def test_custom_and_analytical(model, tmpdir):
    fa = AsdfFile()
    fa.tree["model"] = model
    file_path = str(tmpdir.join("custom_and_analytical_inverse.asdf"))
    fa.write_to(file_path)
    with asdf.open(file_path) as f:
        assert f.tree["model"].inverse is not None


def test_deserialize_compound_user_inverse(tmpdir):
    """
    Confirm that we are able to correctly reconstruct a
    compound model with a user inverse set on one of its
    component models.

    Due to code in TransformType that facilitates circular
    inverses, the user inverse of the component model is
    not available at the time that the CompoundModel is
    constructed.
    """

    yaml = """
model: !transform/concatenate-1.2.0
  forward:
  - !transform/shift-1.2.0
    inverse: !transform/shift-1.2.0 {offset: 5.0}
    offset: -10.0
  - !transform/shift-1.2.0 {offset: -20.0}
  """
    buff = yaml_to_asdf(yaml)
    with asdf.open(buff) as af:
        model = af["model"]
        assert model.has_inverse()
        assert model.inverse(-5, -20) == (0, 0)


# test some models and compound models with some input unit equivalencies
def models_with_input_eq():
    # 1D model
    m1 = astmodels.Shift(1 * u.kg)
    m1.input_units_equivalencies = {"x": u.mass_energy()}

    # 2D model
    m2 = astmodels.Const2D(10 * u.Hz)
    m2.input_units_equivalencies = {
        "x": u.dimensionless_angles(),
        "y": u.dimensionless_angles(),
    }

    # 2D model with only one input equivalencies
    m3 = astmodels.Const2D(10 * u.Hz)
    m3.input_units_equivalencies = {"x": u.dimensionless_angles()}

    # model using equivalency that has args using units
    m4 = astmodels.PowerLaw1D(amplitude=1 * u.m, x_0=10 * u.pix, alpha=7)
    m4.input_units_equivalencies = {
        "x": u.equivalencies.pixel_scale(0.5 * u.arcsec / u.pix)
    }

    return [m1, m2, m3, m4]


def compound_models_with_input_eq():
    m1 = astmodels.Gaussian1D(10 * u.K, 11 * u.arcsec, 12 * u.arcsec)
    m1.input_units_equivalencies = {"x": u.parallax()}
    m2 = astmodels.Gaussian1D(5 * u.s, 2 * u.K, 3 * u.K)
    m2.input_units_equivalencies = {"x": u.temperature()}

    return [m1 | m2, m1 & m2, m1 + m2]


test_models.extend(models_with_input_eq())
test_models.extend(compound_models_with_input_eq())
