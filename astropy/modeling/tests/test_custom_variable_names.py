import pytest
from astropy.modeling import models

def test_custom_variable_names():
    g = models.Gaussian1D(amplitude=1.2, mean=0.5, stddev=0.3, name = "test_gaussian")

    model_str = str(g)

    expected_names = ["test_gaussian_amplitude", "test_gaussian_mean", "test_gaussian_stddev"]
    for name in expected_names:
        assert name in model_str, f"Custom variable name '{name}' not found in model string representation"
