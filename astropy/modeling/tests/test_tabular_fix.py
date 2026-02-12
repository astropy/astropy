import pytest
import numpy as np
from astropy import units as u
from astropy.modeling.tabular import Tabular1D, Tabular2D
from astropy.modeling.mappings import Identity

def test_tabular_unit_validation_fixes():
    # TEST 1: Unitless Tabular1D should strictly reject Quantities
    # (Previously this crashed inside scipy; now it should raise a clean UnitsError)
    points = np.array([0, 1, 2, 3])
    table = np.array([0, 10, 20, 30])
    t_unitless = Tabular1D(points, table)

    # Works with floats
    assert t_unitless(1.5) == 15.0
    
    # Should FAIL validation if passed a unit
    with pytest.raises(u.UnitsError):
        t_unitless(1.5 * u.m)

    # TEST 2: Tabular2D with mixed units (e.g. Distance vs Time)
    # (Previously this failed because it assumed both inputs were 'm')
    p_x = np.array([0, 10]) * u.m
    p_y = np.array([0, 5]) * u.s
    table_2d = np.array([[1, 2], [3, 4]]) * u.K
    
    t_mixed = Tabular2D((p_x, p_y), table_2d)

    # Evaluate at (5m, 2.5s) -> should be center of grid
    result = t_mixed(5 * u.m, 2.5 * u.s)
    
    assert result.unit == u.K
    # Previously, this line would crash claiming 's' is not compatible with 'm'
    assert np.isclose(result.value, 2.5) 

    # Verify input checks work
    with pytest.raises(u.UnitsError):
        t_mixed(5, 2.5) # Missing units

def test_tabular_return_units():
    # Verify return_units property is populated correctly
    t = Tabular1D([1, 2, 3], [1, 2, 3])
    assert t.return_units['y'] == u.dimensionless_unscaled
    
    t_q = Tabular1D([1, 2, 3]*u.m, [1, 2, 3]*u.s)
    assert t_q.return_units['y'] == u.s

def test_identity_units():
    # Identity models should simply pass through whatever units they get
    id_model = Identity(1)
    
    # Pass float -> return float
    assert id_model(10) == 10
    
    # Pass Quantity -> return Quantity (Should NOT crash)
    res = id_model(10 * u.m)
    assert res == 10 * u.m
    assert res.unit == u.m