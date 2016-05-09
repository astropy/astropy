# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to evaluating models with quantity parameters
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from numpy.testing import assert_allclose


from ..core import Fittable1DModel
from ..parameters import Parameter
from ..models import Gaussian1D
from ... import units as u
from ...units import UnitsError
from ...tests.helper import pytest, assert_quantity_allclose


# We start off by taking some simple cases where the units are defined by
# whatever the model is initialized with, and we check that the model evaluation
# returns quantities.


def test_evaluate_with_quantities():
    """
    Test evaluation of a single model with Quantity parameters that do
    not explicitly require units.
    """

    # We create two models here - one with quantities, and one without. The one
    # without is used to create the reference values for comparison.

    g = Gaussian1D(1, 1, 0.1)
    gq = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    # We first check that calling the Gaussian with quantities returns the
    # expected result
    assert_quantity_allclose(gq(1 * u.m), g(1) * u.J)

    # Units have to be specified for the Gaussian with quantities - if not, an
    # error is raised
    with pytest.raises(UnitsError) as exc:
        gq(1)
    assert exc.value.args[0] == ("Units of input 'x', (dimensionless), could not be "
                                 "converted to required input units of m (length)")

    # However, zero is a special case
    assert_quantity_allclose(gq(0), g(0) * u.J)

    # We can also evaluate models with equivalent units
    assert_allclose(gq(0.0005 * u.km).value, g(0.5))

    # But not with incompatible units
    with pytest.raises(UnitsError) as exc:
        gq(3 * u.s)
    assert exc.value.args[0] == ("Units of input 'x', s (time), could not be "
                                 "converted to required input units of m (length)")

    # We also can't evaluate the model without quantities with a quantity
    with pytest.raises(UnitsError) as exc:
        g(3 * u.m)
    assert exc.value.args[0] == ("Units of input 'x', m (length), could not be "
                                 "converted to required dimensionless input")


def test_evaluate_with_quantities_and_equivalencies():
    """
    We now make sure that equivalencies are correctly taken into account
    """

    g = Gaussian1D(1 * u.Jy, 10 * u.nm, 2 * u.nm)

    # We aren't setting the equivalencies, so this won't work
    with pytest.raises(UnitsError) as exc:
        g(30 * u.PHz)
    assert exc.value.args[0] == ("Units of input 'x', PHz (frequency), could "
                                 "not be converted to required input units of "
                                 "nm (length)")

    # But it should now work if we pass equivalencies when evaluating
    assert_quantity_allclose(g(30 * u.PHz, equivalencies=u.spectral()),
                             g(9.993081933333332 * u.nm))
