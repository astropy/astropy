"""
Regression tests for the units package
"""

from __future__ import absolute_import, division, print_function
from astropy import units as u
import pytest

def test_convert():
    assert u.hr.converter_to(u.s)(1) == 3600
    
def test_convert_fail():
    with pytest.raises(u.UnitsException):
        u.cm.convert_to(u.s, 1)
    with pytest.raises(u.UnitsException):
        (u.cm/u.s).convert_to(u.m,1)
            
def test_composite():
    assert (u.cm/u.s*u.hr).converter_to(u.m)(1) == 36
    assert u.cm*u.cm == u.cm**2

    assert u.cm*u.cm*u.cm == u.cm**3

    assert (1 / (u.cm * u.cm)) == u.cm**-2

    #assert u.hz.convert(1000 * u.hz)(1) == 0.001
def test_str():
    assert str(u.cm) == "cm"
    
def test_repr():
    assert repr(u.cm) == 'unit("cm")'
    
def test_spectral():
    assert abs(u.sp_A.convert_to(u.sp_Hz,1) - 2.9979245799999995e+18) < 1000
    assert abs(u.sp_A.converter_to(u.sp_Hz)(1) - 2.9979245799999995e+18) < 1000
    with pytest.raises(TypeError):
        u.sp_A * u.m
        
def test_spectraldensity():
    assert abs(u.sd_A.convert_to(u.Jy,1,u.sp_eV,2.2) - 1059416252057.8357) < 10
