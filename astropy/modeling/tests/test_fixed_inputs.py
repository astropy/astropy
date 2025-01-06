import numpy as np

from astropy.modeling import fitting, models
from astropy.modeling.models import fix_inputs


def test_fix_inputs_fittable():
    """Test that fix_inputs(...).fittable does not raise AttributeError."""
    photon = models.Linear1D(slope=2, intercept=1)
    p = fix_inputs(photon, {"x": 3.0})
    # Prior to the patch, we might get an AttributeError
    # After the patch, we expect this to succeed:
    assert hasattr(p, "fittable")
    assert p.fittable is True


def test_fix_inputs_zero_input_fitting():
    """
    Test if we can fit a model whose entire input is fixed.
    This is not fully supported by default; we pass a dummy `x` or
    see if the fitter can handle no x at all.
    """
    # Our model still has free parameters (slope, intercept),
    # but the x input is fixed to 3.
    photon = models.Linear1D(slope=2, intercept=1)
    p = fix_inputs(photon, {"x": 3.0})

    # Suppose we want to fit for intercept given y_data at that single point:
    y_data = np.array(
        [7.0]
    )  # we want slope*3 + intercept=7 => 2*3 + 1=7 => intercept=1

    fitter = fitting.LMLSQFitter()

    # Typically, `__call__` expects `fitter(model, x, y)`:
    # a minimal workaround is to pass a dummy array for x.
    x_dummy = np.zeros_like(y_data)  # shape matches y_data
    fitted_p = fitter(p, x_dummy, y_data)

    # Check that it recovers slope=2, intercept=1
    # param order for a Linear1D is slope, intercept
    np.testing.assert_allclose(fitted_p.parameters, [2.0, 1.0], rtol=1e-6)
