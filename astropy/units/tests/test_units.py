"""
Regression tests for the units package
"""


from __future__ import absolute_import, division, print_function
from astropy import units as u
import pytest

def test_convert():
    assert u.hour.convert(u.s)(1) == 3600
    
def test_convert_fail():
    with pytest.raises(ValueError):
        u.cm.convert(u.s)
    with pytest.raises(ValueError):
        (u.cm/u.s).convert(u.meter)
            
def test_compound():
    assert (u.cm/u.s*u.hour).scale == 36
    assert (u.cm/u.s*u.hour).convert(u.meter)(1) == 36
    assert u.cm*u.cm == u.cm**2
    #assert u1.num_units == u2.num_units and u1.denom_units == u2.denom_units and u1.scale == u2.scale
    assert u.cm*u.cm*u.cm == u.cm**3
    #assert u1.num_units == u2.num_units and u1.denom_units == u2.denom_units and u1.scale == u2.scale
    assert (1 / (u.cm * u.cm)) == u.cm**-2
    #assert u1.num_units == u2.num_units and u1.denom_units == u2.denom_units and u1.scale == u2.scale
    
def test_str():
    assert str(u.cm) == "Units: centimeter"
    
def test_repr():
    assert repr(u.cm) == "Length(scale=1.000000e-02, name='centimeter')"