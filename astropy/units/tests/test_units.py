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
    
def test_str():
    assert str(u.cm) == "Units: centimeter"
    
def test_repr():
    assert repr(u.cm) == "Length(scale=1.000000e-02, name='centimeter')"