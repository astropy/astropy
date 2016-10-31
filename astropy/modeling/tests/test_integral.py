# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from numpy.testing import assert_allclose

from .. import models
from ...tests.helper import pytest


def test_no_integral():
    """Test that undefined integral raises appropriate exception.

    .. note::

        If the model being tested here has integral in the future, change
        the model being tested to one without integral.

    """
    model = models.RedshiftScaleFactor()
    with pytest.raises(NotImplementedError):
        model.integral


def test_Box1D():
    """Test Box1D integration.

    This tests most of the integral functionality.
    Tests for other models can probably just test that the
    integration is working properly for the "most popular" use case.

    """
    model = models.Box1D(amplitude=3., x_0=0., width=10.)
    bbox = model.bounding_box
    intg_model = model.integral

    assert isinstance(intg_model, models.Ramp1D)

    # Integrate over the whole range.
    y = intg_model(bbox)
    intg_result = y[1] - y[0]
    assert_allclose(intg_result, 30)

    # Partial integration
    assert_allclose(intg_model(6) - intg_model(0), 15)
    assert_allclose(intg_model(20) - intg_model(10), 0)

    # User-defined integral
    model.integral = models.Const1D(100)
    assert model.has_user_integral
    assert model.integral(20) == 100

    model.integral = None
    with pytest.raises(NotImplementedError):
        model.integral

    del model.integral
    assert not model.has_user_integral
    assert model.integral(20) == 30

    # n_models > 1
    model = models.Box1D(amplitude=[3, 3], x_0=[0, 4], width=[10, 5])
    assert_allclose(model.integral([0, 4]), [15, 7.5])
    y = model.integral(model.bounding_box)
    assert_allclose(y[1] - y[0], [30, 15])
