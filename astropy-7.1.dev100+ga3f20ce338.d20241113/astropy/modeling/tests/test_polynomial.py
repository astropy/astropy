# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for polynomial models."""

# pylint: disable=invalid-name
import os
import unittest.mock as mk
import warnings
from itertools import product

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import conf, wcs
from astropy.io import fits
from astropy.modeling import fitting
from astropy.modeling.functional_models import Linear1D
from astropy.modeling.mappings import Identity
from astropy.modeling.polynomial import (
    SIP,
    Chebyshev1D,
    Chebyshev2D,
    Hermite1D,
    Hermite2D,
    Legendre1D,
    Legendre2D,
    OrthoPolynomialBase,
    Polynomial1D,
    Polynomial2D,
    PolynomialBase,
)
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyUserWarning

linear1d = {
    Chebyshev1D: {
        "args": (3,),
        "kwargs": {"domain": [1, 10]},
        "parameters": {"c0": 1.2, "c1": 2, "c2": 2.3, "c3": 0.2},
        "constraints": {"fixed": {"c0": True}},
    },
    Hermite1D: {
        "args": (3,),
        "kwargs": {"domain": [1, 10]},
        "parameters": {"c0": 1.2, "c1": 2, "c2": 2.3, "c3": 0.2},
        "constraints": {"fixed": {"c0": True}},
    },
    Legendre1D: {
        "args": (3,),
        "kwargs": {"domain": [1, 10]},
        "parameters": {"c0": 1.2, "c1": 2, "c2": 2.3, "c3": 0.2},
        "constraints": {"fixed": {"c0": True}},
    },
    Polynomial1D: {
        "args": (3,),
        "kwargs": {"domain": [1, 10]},
        "parameters": {"c0": 1.2, "c1": 2, "c2": 2.3, "c3": 0.2},
        "constraints": {"fixed": {"c0": True}},
    },
    Linear1D: {
        "args": (),
        "kwargs": {},
        "parameters": {"intercept": 1.2, "slope": 23.1},
        "constraints": {"fixed": {"intercept": True}},
    },
}


linear2d = {
    Chebyshev2D: {
        "args": (1, 1),
        "kwargs": {"x_domain": [0, 99], "y_domain": [0, 82]},
        "parameters": {"c0_0": 1.2, "c1_0": 2, "c0_1": 2.3, "c1_1": 0.2},
        "constraints": {"fixed": {"c0_0": True}},
    },
    Hermite2D: {
        "args": (1, 1),
        "kwargs": {"x_domain": [0, 99], "y_domain": [0, 82]},
        "parameters": {"c0_0": 1.2, "c1_0": 2, "c0_1": 2.3, "c1_1": 0.2},
        "constraints": {"fixed": {"c0_0": True}},
    },
    Legendre2D: {
        "args": (1, 1),
        "kwargs": {"x_domain": [0, 99], "y_domain": [0, 82]},
        "parameters": {"c0_0": 1.2, "c1_0": 2, "c0_1": 2.3, "c1_1": 0.2},
        "constraints": {"fixed": {"c0_0": True}},
    },
    Polynomial2D: {
        "args": (1,),
        "kwargs": {},
        "parameters": {"c0_0": 1.2, "c1_0": 2, "c0_1": 2.3},
        "constraints": {"fixed": {"c0_0": True}},
    },
}

fitters = [
    fitting.LevMarLSQFitter,
    fitting.TRFLSQFitter,
    fitting.LMLSQFitter,
    fitting.DogBoxLSQFitter,
]


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
class TestFitting:
    """Test linear fitter with polynomial models."""

    def setup_class(self):
        self.N = 100
        self.M = 100
        self.x1 = np.linspace(1, 10, 100)
        self.y2, self.x2 = np.mgrid[:100, :83]
        rsn = np.random.default_rng(0)
        self.n1 = rsn.standard_normal(self.x1.size) * 0.1
        self.n2 = rsn.standard_normal(self.x2.size)
        self.n2.shape = self.x2.shape
        self.linear_fitter = fitting.LinearLSQFitter()

    # TODO: Most of these test cases have some pretty repetitive setup that we
    # could probably factor out

    @pytest.mark.parametrize(
        ("model_class", "constraints"),
        list(product(sorted(linear1d, key=str), (False, True))),
    )
    def test_linear_fitter_1D(self, model_class, constraints):
        """Test fitting with LinearLSQFitter"""

        model_args = linear1d[model_class]
        kwargs = {}
        kwargs.update(model_args["kwargs"])
        kwargs.update(model_args["parameters"])

        if constraints:
            kwargs.update(model_args["constraints"])

        model = model_class(*model_args["args"], **kwargs)

        y1 = model(self.x1)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r"The fit may be poorly conditioned",
                category=AstropyUserWarning,
            )
            model_lin = self.linear_fitter(model, self.x1, y1 + self.n1)

        if constraints:
            # For the constraints tests we're not checking the overall fit,
            # just that the constraint was maintained
            fixed = model_args["constraints"].get("fixed", None)
            if fixed:
                for param in fixed:
                    expected = model_args["parameters"][param]
                    assert getattr(model_lin, param).value == expected
        else:
            assert_allclose(model_lin.parameters, model.parameters, atol=0.2)

    @pytest.mark.parametrize(
        ("model_class", "constraints"),
        list(product(sorted(linear1d, key=str), (False, True))),
    )
    @pytest.mark.parametrize("fitter", fitters)
    def test_non_linear_fitter_1D(self, model_class, constraints, fitter):
        """Test fitting with non-linear LevMarLSQFitter"""
        fitter = fitter()

        model_args = linear1d[model_class]
        kwargs = {}
        kwargs.update(model_args["kwargs"])
        kwargs.update(model_args["parameters"])

        if constraints:
            kwargs.update(model_args["constraints"])

        model = model_class(*model_args["args"], **kwargs)

        y1 = model(self.x1)
        with pytest.warns(AstropyUserWarning, match="Model is linear in parameters"):
            model_nlin = fitter(model, self.x1, y1 + self.n1)

        if constraints:
            fixed = model_args["constraints"].get("fixed", None)
            if fixed:
                for param in fixed:
                    expected = model_args["parameters"][param]
                    assert getattr(model_nlin, param).value == expected
        else:
            assert_allclose(model_nlin.parameters, model.parameters, atol=0.2)

    @pytest.mark.parametrize(
        ("model_class", "constraints"),
        list(product(sorted(linear2d, key=str), (False, True))),
    )
    def test_linear_fitter_2D(self, model_class, constraints):
        """Test fitting with LinearLSQFitter"""

        model_args = linear2d[model_class]
        kwargs = {}
        kwargs.update(model_args["kwargs"])
        kwargs.update(model_args["parameters"])

        if constraints:
            kwargs.update(model_args["constraints"])

        model = model_class(*model_args["args"], **kwargs)

        z = model(self.x2, self.y2)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r"The fit may be poorly conditioned",
                category=AstropyUserWarning,
            )
            model_lin = self.linear_fitter(model, self.x2, self.y2, z + self.n2)

        if constraints:
            fixed = model_args["constraints"].get("fixed", None)
            if fixed:
                for param in fixed:
                    expected = model_args["parameters"][param]
                    assert getattr(model_lin, param).value == expected
        else:
            assert_allclose(model_lin.parameters, model.parameters, atol=0.2)

    @pytest.mark.parametrize(
        ("model_class", "constraints"),
        list(product(sorted(linear2d, key=str), (False, True))),
    )
    @pytest.mark.parametrize("fitter", fitters)
    def test_non_linear_fitter_2D(self, model_class, constraints, fitter):
        """Test fitting with non-linear LevMarLSQFitter"""
        fitter = fitter()

        model_args = linear2d[model_class]
        kwargs = {}
        kwargs.update(model_args["kwargs"])
        kwargs.update(model_args["parameters"])

        if constraints:
            kwargs.update(model_args["constraints"])

        model = model_class(*model_args["args"], **kwargs)

        z = model(self.x2, self.y2)
        with pytest.warns(AstropyUserWarning, match="Model is linear in parameters"):
            model_nlin = fitter(model, self.x2, self.y2, z + self.n2)

        if constraints:
            fixed = model_args["constraints"].get("fixed", None)
            if fixed:
                for param in fixed:
                    expected = model_args["parameters"][param]
                    assert getattr(model_nlin, param).value == expected
        else:
            assert_allclose(model_nlin.parameters, model.parameters, atol=0.2)


@pytest.mark.parametrize("model_class", list(list(linear1d) + list(linear2d)))
def test_polynomial_init_with_constraints(model_class):
    """
    Test that polynomial models can be instantiated with constraints, but no
    parameters specified.

    Regression test for https://github.com/astropy/astropy/issues/3606
    """

    # Just determine which parameter to place a constraint on; it doesn't
    # matter which parameter it is to exhibit the problem so long as it's a
    # valid parameter for the model
    if "1D" in model_class.__name__:
        param = "c0"
    else:
        param = "c0_0"

    if issubclass(model_class, Linear1D):
        param = "intercept"

    if issubclass(model_class, OrthoPolynomialBase):
        degree = (2, 2)
    else:
        degree = (2,)

    m = model_class(*degree, fixed={param: True})

    assert m.fixed[param] is True
    assert getattr(m, param).fixed is True

    if issubclass(model_class, OrthoPolynomialBase):
        assert (
            repr(m)
            == f"<{model_class.__name__}(2, 2, c0_0=0., c1_0=0., c2_0=0., c0_1=0., "
            "c1_1=0., c2_1=0., c0_2=0., c1_2=0., c2_2=0.)>"
        )
        assert (
            str(m) == f"Model: {model_class.__name__}\n"
            "Inputs: ('x', 'y')\n"
            "Outputs: ('z',)\n"
            "Model set size: 1\n"
            "X_Degree: 2\n"
            "Y_Degree: 2\n"
            "Parameters:\n"
            "    c0_0 c1_0 c2_0 c0_1 c1_1 c2_1 c0_2 c1_2 c2_2\n"
            "    ---- ---- ---- ---- ---- ---- ---- ---- ----\n"
            "     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0"
        )
    else:
        if model_class.__name__ == "Polynomial2D":
            assert (
                repr(m) == "<Polynomial2D(2, c0_0=0., c1_0=0., c2_0=0., "
                "c0_1=0., c0_2=0., c1_1=0.)>"
            )
            assert (
                str(m) == "Model: Polynomial2D\n"
                "Inputs: ('x', 'y')\n"
                "Outputs: ('z',)\n"
                "Model set size: 1\n"
                "Degree: 2\n"
                "Parameters:\n"
                "    c0_0 c1_0 c2_0 c0_1 c0_2 c1_1\n"
                "    ---- ---- ---- ---- ---- ----\n"
                "     0.0  0.0  0.0  0.0  0.0  0.0"
            )
        elif model_class.__name__ == "Linear1D":
            assert repr(m) == "<Linear1D(slope=2., intercept=0.)>"
            assert (
                str(m) == "Model: Linear1D\n"
                "Inputs: ('x',)\n"
                "Outputs: ('y',)\n"
                "Model set size: 1\n"
                "Parameters:\n"
                "    slope intercept\n"
                "    ----- ---------\n"
                "      2.0       0.0"
            )
        else:
            assert repr(m) == f"<{model_class.__name__}(2, c0=0., c1=0., c2=0.)>"
            assert (
                str(m) == f"Model: {model_class.__name__}\n"
                "Inputs: ('x',)\n"
                "Outputs: ('y',)\n"
                "Model set size: 1\n"
                "Degree: 2\n"
                "Parameters:\n"
                "     c0  c1  c2\n"
                "    --- --- ---\n"
                "    0.0 0.0 0.0"
            )


def test_sip_hst():
    """Test SIP against astropy.wcs"""

    test_file = get_pkg_data_filename(os.path.join("data", "hst_sip.hdr"))
    hdr = fits.Header.fromtextfile(test_file)
    crpix1 = hdr["CRPIX1"]
    crpix2 = hdr["CRPIX2"]
    wobj = wcs.WCS(hdr)
    a_pars = dict(**hdr["A_*"])
    b_pars = dict(**hdr["B_*"])
    a_order = a_pars.pop("A_ORDER")
    b_order = b_pars.pop("B_ORDER")
    sip = SIP([crpix1, crpix2], a_order, b_order, a_pars, b_pars)
    coords = [1, 1]
    rel_coords = [1 - crpix1, 1 - crpix2]
    astwcs_result = wobj.sip_pix2foc([coords], 1)[0] - rel_coords
    assert_allclose(sip(1, 1), astwcs_result)

    # Test changing of inputs and calling it with keyword argumenrts.
    sip.inputs = ("r", "t")
    assert_allclose(sip(r=1, t=1), astwcs_result)
    assert_allclose(sip(1, t=1), astwcs_result)

    # Test representations
    assert (
        repr(sip) == "<SIP([<Shift(offset=-2048.)>, <Shift(offset=-1024.)>, "
        "<_SIP1D(4, 'A', A_2_0=0.00000855, A_3_0=-0., A_4_0=0., A_0_2=0.00000217, "
        "A_0_3=0., A_0_4=0., A_1_1=-0.0000052, A_1_2=-0., A_1_3=-0., "
        "A_2_1=-0., A_2_2=0., A_3_1=0.)>, "
        "<_SIP1D(4, 'B', B_2_0=-0.00000175, B_3_0=0., B_4_0=-0., B_0_2=-0.00000722, "
        "B_0_3=-0., B_0_4=-0., B_1_1=0.00000618, B_1_2=-0., B_1_3=0., "
        "B_2_1=-0., B_2_2=-0., B_3_1=-0.)>])>"
    )
    with conf.set_temp("max_width", 80):
        # fmt: off
        assert str(sip) == (
            "Model: SIP\n"
            "    Model: Shift\n"
            "    Inputs: ('x',)\n"
            "    Outputs: ('y',)\n"
            "    Model set size: 1\n"
            "    Parameters:\n"
            "         offset\n"
            "        -------\n"
            "        -2048.0\n"
            "\n"
            "    Model: Shift\n"
            "    Inputs: ('x',)\n"
            "    Outputs: ('y',)\n"
            "    Model set size: 1\n"
            "    Parameters:\n"
            "         offset\n"
            "        -------\n"
            "        -1024.0\n"
            "\n"
            "    Model: _SIP1D\n"
            "    Inputs: ('x', 'y')\n"
            "    Outputs: ('z',)\n"
            "    Model set size: 1\n"
            "    Order: 4\n"
            "    Coeff. Prefix: A\n"
            "    Parameters:\n"
            "                A_2_0                 A_3_0          ...         A_3_1        \n"
            "        --------------------- ---------------------- ... ---------------------\n"
            "        8.551277582556502e-06 -4.730444829222791e-10 ... 1.971022971660309e-15\n"
            "\n"
            "    Model: _SIP1D\n"
            "    Inputs: ('x', 'y')\n"
            "    Outputs: ('z',)\n"
            "    Model set size: 1\n"
            "    Order: 4\n"
            "    Coeff. Prefix: B\n"
            "    Parameters:\n"
            "                B_2_0                  B_3_0         ...         B_3_1         \n"
            "        ---------------------- --------------------- ... ----------------------\n"
            "        -1.746491877058669e-06 8.567635427816317e-11 ... -3.779506805487476e-15\n"
        )
        # fmt: on

    # Test get num of coeffs
    assert sip.sip1d_a.get_num_coeff(1) == 6
    # Test error
    MESSAGE = "Degree of polynomial must be 2< deg < 9"
    sip.sip1d_a.order = 1
    with pytest.raises(ValueError, match=MESSAGE):
        sip.sip1d_a.get_num_coeff(1)
    sip.sip1d_a.order = 10
    with pytest.raises(ValueError, match=MESSAGE):
        sip.sip1d_a.get_num_coeff(1)


def test_sip_irac():
    """Test forward and inverse SIP against astropy.wcs"""

    test_file = get_pkg_data_filename(os.path.join("data", "irac_sip.hdr"))
    hdr = fits.Header.fromtextfile(test_file)
    crpix1 = hdr["CRPIX1"]
    crpix2 = hdr["CRPIX2"]
    wobj = wcs.WCS(hdr)
    a_pars = dict(**hdr["A_*"])
    b_pars = dict(**hdr["B_*"])
    ap_pars = dict(**hdr["AP_*"])
    bp_pars = dict(**hdr["BP_*"])
    a_order = a_pars.pop("A_ORDER")
    b_order = b_pars.pop("B_ORDER")
    ap_order = ap_pars.pop("AP_ORDER")
    bp_order = bp_pars.pop("BP_ORDER")
    del a_pars["A_DMAX"]
    del b_pars["B_DMAX"]
    pix = [200, 200]
    rel_pix = [200 - crpix1, 200 - crpix2]
    sip = SIP(
        [crpix1, crpix2],
        a_order,
        b_order,
        a_pars,
        b_pars,
        ap_order=ap_order,
        ap_coeff=ap_pars,
        bp_order=bp_order,
        bp_coeff=bp_pars,
    )

    foc = wobj.sip_pix2foc([pix], 1)
    newpix = wobj.sip_foc2pix(foc, 1)[0]
    assert_allclose(sip(*pix), foc[0] - rel_pix)
    assert_allclose(sip.inverse(*foc[0]) + foc[0] - rel_pix, newpix - pix)

    # Test inverse representations
    assert (
        repr(sip.inverse)
        == "<InverseSIP([<Polynomial2D(2, c0_0=0., c1_0=0.0000114, c2_0=0.00002353, "
        "c0_1=-0.00000546, c0_2=-0.00000667, c1_1=-0.00001801)>, "
        "<Polynomial2D(2, c0_0=0., c1_0=-0.00001495, c2_0=0.00000122, c0_1=0.00001975, "
        "c0_2=-0.00002601, c1_1=0.00002944)>])>"
    )
    assert (
        str(sip.inverse) == "Model: InverseSIP\n"
        "    Model: Polynomial2D\n"
        "    Inputs: ('x', 'y')\n"
        "    Outputs: ('z',)\n"
        "    Model set size: 1\n"
        "    Degree: 2\n"
        "    Parameters:\n"
        "        c0_0   c1_0      c2_0      c0_1       c0_2       c1_1   \n"
        "        ---- -------- --------- ---------- ---------- ----------\n"
        "         0.0 1.14e-05 2.353e-05 -5.463e-06 -6.666e-06 -1.801e-05\n"
        "\n"
        "    Model: Polynomial2D\n"
        "    Inputs: ('x', 'y')\n"
        "    Outputs: ('z',)\n"
        "    Model set size: 1\n"
        "    Degree: 2\n"
        "    Parameters:\n"
        "        c0_0    c1_0       c2_0      c0_1      c0_2       c1_1  \n"
        "        ---- ---------- --------- --------- ---------- ---------\n"
        "         0.0 -1.495e-05 1.225e-06 1.975e-05 -2.601e-05 2.944e-05\n"
    )


def test_sip_no_coeff():
    sip = SIP([10, 12], 2, 2)
    assert_allclose(sip.sip1d_a.parameters, [0.0, 0.0, 0])
    assert_allclose(sip.sip1d_b.parameters, [0.0, 0.0, 0])
    MESSAGE = r"SIP inverse coefficients are not available"
    with pytest.raises(NotImplementedError, match=MESSAGE):
        sip.inverse

    # Test model set
    sip = SIP([10, 12], 2, 2, n_models=2)
    assert sip.sip1d_a.model_set_axis == 0
    assert sip.sip1d_b.model_set_axis == 0


@pytest.mark.parametrize(
    "cls",
    (Polynomial1D, Chebyshev1D, Legendre1D, Polynomial2D, Chebyshev2D, Legendre2D),
)
def test_zero_degree_polynomial(cls):
    """
    A few tests that degree=0 polynomials are correctly evaluated and
    fitted.

    Regression test for https://github.com/astropy/astropy/pull/3589
    """

    MESSAGE = "Degree of polynomial must be positive or null"

    if cls.n_inputs == 1:  # Test 1D polynomials
        p1 = cls(degree=0, c0=1)
        assert p1(0) == 1
        assert np.all(p1(np.zeros(5)) == np.ones(5))

        x = np.linspace(0, 1, 100)
        # Add a little noise along a straight line
        y = 1 + np.random.uniform(0, 0.1, len(x))

        p1_init = cls(degree=0)
        fitter = fitting.LinearLSQFitter()
        p1_fit = fitter(p1_init, x, y)

        # The fit won't be exact of course, but it should get close to within
        # 1%
        assert_allclose(p1_fit.c0, 1, atol=0.10)

        # Error from negative degree
        with pytest.raises(ValueError, match=MESSAGE):
            cls(degree=-1)
    elif cls.n_inputs == 2:  # Test 2D polynomials
        if issubclass(cls, OrthoPolynomialBase):
            p2 = cls(x_degree=0, y_degree=0, c0_0=1)

            # different shaped x and y inputs
            a = np.array([1, 2, 3])
            b = np.array([1, 2])
            with mk.patch.object(
                PolynomialBase,
                "prepare_inputs",
                autospec=True,
                return_value=((a, b), mk.MagicMock()),
            ):
                with pytest.raises(
                    ValueError, match=r"Expected input arrays to have the same shape"
                ):
                    p2.prepare_inputs(mk.MagicMock(), mk.MagicMock())

            # Error from negative degree
            with pytest.raises(ValueError, match=MESSAGE):
                cls(x_degree=-1, y_degree=0)
            with pytest.raises(ValueError, match=MESSAGE):
                cls(x_degree=0, y_degree=-1)
        else:
            p2 = cls(degree=0, c0_0=1)

            # Error from negative degree
            with pytest.raises(ValueError, match=MESSAGE):
                cls(degree=-1)

        assert p2(0, 0) == 1
        assert np.all(p2(np.zeros(5), np.zeros(5)) == np.ones(5))

        y, x = np.mgrid[0:1:100j, 0:1:100j]
        z = (1 + np.random.uniform(0, 0.1, x.size)).reshape(100, 100)

        if issubclass(cls, OrthoPolynomialBase):
            p2_init = cls(x_degree=0, y_degree=0)
        else:
            p2_init = cls(degree=0)
        fitter = fitting.LinearLSQFitter()
        p2_fit = fitter(p2_init, x, y, z)

        assert_allclose(p2_fit.c0_0, 1, atol=0.10)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.parametrize("fitter", fitters)
def test_2d_orthopolynomial_in_compound_model(fitter):
    """
    Ensure that OrthoPolynomialBase (ie. Chebyshev2D & Legendre2D) models get
    evaluated & fitted correctly when part of a compound model.

    Regression test for https://github.com/astropy/astropy/pull/6085.
    """
    fitter = fitter()

    y, x = np.mgrid[0:5, 0:5]
    z = x + y

    simple_model = Chebyshev2D(2, 2)
    with pytest.warns(AstropyUserWarning, match="Model is linear in parameters"):
        simple_fit = fitter(simple_model, x, y, z)

    compound_model = Identity(2) | Chebyshev2D(2, 2)
    compound_model.fittable = True
    compound_model.linear = True
    with pytest.warns(AstropyUserWarning, match="Model is linear in parameters"):
        compound_fit = fitter(compound_model, x, y, z)

    assert_allclose(simple_fit(x, y), compound_fit(x, y), atol=1e-11)


def test_Hermite1D_clenshaw():
    model = Hermite1D(degree=2)

    assert model.clenshaw(1, [3]) == 3
    assert model.clenshaw(1, [3, 4]) == 11
    assert model.clenshaw(1, [3, 4, 5]) == 21
    assert model.clenshaw(1, [3, 4, 5, 6]) == -3


def test__fcache():
    model = OrthoPolynomialBase(x_degree=2, y_degree=2)
    MESSAGE = r"Subclasses should implement this"
    with pytest.raises(NotImplementedError, match=MESSAGE):
        model._fcache(np.asanyarray(1), np.asanyarray(1))

    model = Hermite2D(x_degree=2, y_degree=2)
    assert model._fcache(np.asanyarray(1), np.asanyarray(1)) == {
        0: np.asanyarray(1),
        1: 2,
        3: np.asanyarray(1),
        4: 2,
        2: 2.0,
        5: -4.0,
    }

    model = Legendre2D(x_degree=2, y_degree=2)
    assert model._fcache(np.asanyarray(1), np.asanyarray(1)) == {
        0: np.asanyarray(1),
        1: np.asanyarray(1),
        2: 1.0,
        3: np.asanyarray(1),
        4: np.asanyarray(1),
        5: 1.0,
    }

    model = Chebyshev2D(x_degree=2, y_degree=2)
    assert model._fcache(np.asanyarray(1), np.asanyarray(1)) == {
        0: np.asanyarray(1),
        1: np.asanyarray(1),
        2: 1.0,
        3: np.asanyarray(1),
        4: np.asanyarray(1),
        5: 1.0,
    }


def test_fit_deriv_shape_error():
    model = Hermite2D(x_degree=2, y_degree=2)
    MESSAGE = r"x and y must have the same shape"
    with pytest.raises(ValueError, match=MESSAGE):
        model.fit_deriv(np.array([1, 2]), np.array([3, 4, 5]))

    model = Chebyshev2D(x_degree=2, y_degree=2)
    with pytest.raises(ValueError, match=MESSAGE):
        model.fit_deriv(np.array([1, 2]), np.array([3, 4, 5]))

    model = Legendre2D(x_degree=2, y_degree=2)
    with pytest.raises(ValueError, match=MESSAGE):
        model.fit_deriv(np.array([1, 2]), np.array([3, 4, 5]))

    model = Polynomial2D(degree=2)
    MESSAGE = r"Expected x and y to be of equal size"
    with pytest.raises(ValueError, match=MESSAGE):
        model.fit_deriv(np.array([1, 2]), np.array([3, 4, 5]))
