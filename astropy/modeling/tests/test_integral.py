# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from numpy.testing import assert_allclose

from .. import models
from ...tests.helper import pytest


def test_no_primitive():
    """Test that undefined primitive raises appropriate exception.

    .. note::

        If the model being tested here has primitive in the future, change
        the model being tested to one without primitive.

    """
    model = models.RedshiftScaleFactor()
    with pytest.raises(NotImplementedError):
        model.primitive
    with pytest.raises(NotImplementedError):
        model.integrate(1, 100)


def test_Box1D():
    """Test Box1D integration.

    This tests most of the primitive and integrate functionality.
    Tests for other models can probably just test that the
    integration is working properly for the "most popular" use case.

    """
    model = models.Box1D(amplitude=3., x_0=0., width=10.)

    assert isinstance(model.primitive, models.Ramp1D)

    # Integrate over the whole range.
    assert_allclose(model.integrate(*model.bounding_box), 30)

    # Partial integration
    assert_allclose(model.integrate(0, 6), 15)
    assert_allclose(model.integrate(10, 20), 0)

    # User-defined primitive
    model.primitive = models.Const1D(100)
    assert model.has_user_primitive
    assert_allclose(model.primitive(20), 100)
    assert_allclose(model.integrate(0, 6), 0)

    # Disable primitive
    model.primitive = None
    with pytest.raises(NotImplementedError):
        model.primitive

    # Restore default primitive
    del model.primitive
    assert not model.has_user_primitive
    assert_allclose(model.primitive(20), 30)

    # n_models > 1
    model = models.Box1D(amplitude=[3, 3], x_0=[0, 4], width=[10, 5])
    assert_allclose(model.primitive([0, 4]), [15, 7.5])
    assert_allclose(model.integrate(*model.bounding_box), [30, 15])
