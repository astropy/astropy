# Licensed under a 3-clause BSD style license - see LICENSE.rst


from .. import Constant,si,cgs

def test_c():    
    #no one is going to be updating c to be *too* much off from 3e8...
    assert abs(si.c - 3e8)<1e7
    assert abs(cgs.c - 3e10)<1e9
    
    #make sure it has the necessary attributes and they're not blank
    assert si.c.error==0 #c is a *defined* quantity
    assert si.c.name
    assert si.c.origin
    assert si.c.units
    assert cgs.c.error==0 #c is a *defined* quantity
    assert cgs.c.name
    assert cgs.c.origin
    assert cgs.c.units
    
def test_h():    
    #check that the float-translation stuff works as expected
    assert si.h.real == si.h*1
    assert cgs.h.real == cgs.h-0
    
    #make sure it has the necessary attributes and they're not blank
    assert si.h.error
    assert si.h.name
    assert si.h.origin
    assert si.h.units
    assert cgs.h.error
    assert cgs.h.name
    assert cgs.h.origin
    assert cgs.h.units