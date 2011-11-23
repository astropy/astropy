# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ...tests.helper import pytest
from ..arrayastrotime import ArrayAstroTime
import datetime

def test_arrayastrotime_jd():
    from numpy import array,all,abs
    
    jdn = 2451545 #surround J2000
    jds1 = [.5,.2]
    jds2 = [.55,.9]
   
    aat1 = ArrayAstroTime(jdn,jds1)
    aat2 = ArrayAstroTime(jdn,jds2)
    
    d1 = aat1.jd - array([ 2451545.5,  2451545.2])
    d2 = aat2.mjd - array([ 51545.05,  51545.4 ])
    assert all(abs(d1)<1e-9)
    assert all(abs(d2)<1e-9)
    
def test_arrayastrotime_epoch():
    from numpy import array,all,abs
    
    jdn = 2451545 #surround J2000
    jds1 = [.5,.2]
    jds2 = [.55,.9]
   
    aat1 = ArrayAstroTime(jdn,jds1)
    aat2 = ArrayAstroTime(jdn,jds2)
    
    d1 = aat1.jepoch - array([ 2000.00136893,  2000.00054757])
    d2 = aat2.bepoch - array([ 2000.00278336,  2000.00374163])
    assert all(abs(d1)<1e-8)
    assert all(abs(d2)<1e-8)
    

def test_diffarrayastrotime():
    pytest.skip('diff of ArrayAstroTime broken - should be fixed')
    from numpy import array,all
    
    jdn = 2451545 #surround J2000
    jds1 = [.5,.2]
    jds2 = [.55,.9]
   
    aat1 = ArrayAstroTime(jdn,jds1)
    aat2 = ArrayAstroTime(jdn,jds2)
    
    dts = aat1 - aat2
    assert aat1 == dts + aat2