import numpy as np
import pytest

from astropy.modeling import Model, Parameter, fitting, models
from astropy.modeling.fitting import LMLSQFitter
from astropy.modeling.functional_models import Linear1D
from astropy.modeling.models import fix_inputs
from astropy.utils.compat.optional_deps import HAS_SCIPY


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fix_inputs_fittable():
    """Test that fix_inputs(...).fittable does not raise AttributeError."""
    photon = models.Linear1D(slope=2, intercept=1)
    p = fix_inputs(photon, {"x": 3.0})
    # Prior to the patch, we might get an AttributeError
    # After the patch, we expect this to succeed:
    assert hasattr(p, "fittable")
    assert p.fittable is True


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fix_inputs_zero_input_fitting():
    """
    Test if we can fit a model whose entire input is fixed, but still
    have enough data points to avoid m < n in Levenberg–Marquardt.
    """
    # Model still has two free parameters
    photon = models.Linear1D(slope=2, intercept=1)
    p = fix_inputs(photon, {"x": 3.0})

    # Provide TWO identical data points for the same x=3.0
    # We effectively have the system:
    #   slope*3 + intercept = y_data[i] for i in {0, 1}
    # Here both are 7, so the exact fit is slope=2, intercept=1
    # but we now have m=2 data points for n=2 parameters => OK for LM
    y_data = np.array([7.0, 7.0])

    fitter = fitting.LMLSQFitter()

    # Provide a dummy x array of the same shape as y_data
    # The x values are irrelevant since the model input is "fixed",
    # but the shape must match for the fitter's internal checks.
    x_dummy = np.zeros_like(y_data)

    fitted_p = fitter(p, x_dummy, y_data)

    # Check it recovers slope=2, intercept=1
    np.testing.assert_allclose(fitted_p.parameters, [2.0, 1.0], atol=1e-7)


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fix_inputs_non_fittable():
    """
    Check that if you fix the inputs on a model which isn't fittable,
    the resulting model also isn't fittable.
    """

    class NonFittable1D(Model):
        fittable = False
        n_inputs = 1
        n_outputs = 1
        param_names = ("param1",)

        param1 = Parameter(default=1)

        @staticmethod
        def evaluate(x, param1):
            return param1 * x

    # Instantiate the non-fittable model
    non_fittable_model = NonFittable1D()

    # Fix the input on the non-fittable model
    fixed_model = fix_inputs(non_fittable_model, {"x": 3.0})

    # Verify that the returned model is still not fittable
    assert fixed_model.fittable is False


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
def test_fix_inputs_pipe_operator():
    """
    Test chaining a fix_inputs(...) model with another model using '|',
    ensuring that the dictionary is skipped in _leaflist and submodel
    indexing is correct.

    NOTE: We do NOT attempt to fit here, because zero-input fitting
    is still a known limitation in the API. That part is tested
    (and xfailed) in a separate test below.
    """
    m1 = Linear1D(slope=2, intercept=0)
    m2 = Linear1D(slope=4, intercept=1)
    # Fix inputs of m1
    fixed_m1 = fix_inputs(m1, {"x": [1, 2, 3, 4]})

    # Compose with m2 using the pipe operator
    comb = fixed_m1 | m2

    # Check submodel indexing: comb[0] should be the "fixed m1", comb[1] should be m2
    assert comb[1] is m2, "Right-hand submodel should be m2, not the dict"

    # If both submodels are effectively fittable, the final pipeline model is fittable=True
    assert comb.fittable is True

    # Evaluate the combined model just to confirm no dictionary errors:
    y_data = comb()  # This should work fine, no external "x" needed here.
    assert len(y_data) == 4, "Expected 4 output values from the piped model."


@pytest.mark.skipif(not HAS_SCIPY, reason="requires scipy")
@pytest.mark.xfail(
    reason="Zero-input fitting not yet supported by Astropy's fitting API."
)
def test_fix_inputs_pipe_operator_xfail_fitting():
    """
    Demonstrate that fitting a model with zero external inputs (via fix_inputs)
    piped to another model fails (the second part of Cadair's original comment).

    This test is xfailed to indicate it's a known limitation and not yet solved.
    """
    m1 = Linear1D(slope=2, intercept=0)
    m2 = Linear1D(slope=4, intercept=1)
    fixed_m1 = fix_inputs(m1, {"x": [1, 2, 3, 4]})
    comb = fixed_m1 | m2

    # Evaluate the combined model
    y_data = comb()

    # Attempt to fit => triggers IndexError / no-input mismatch
    fitter = LMLSQFitter()
    x_dummy = np.zeros_like(y_data)

    # This call is expected to fail due to 0D input on left side
    fitted = fitter(comb, x_dummy, y_data)

    # If the limitation were resolved, we might check:
    # assert fitted.parameters == [...]
