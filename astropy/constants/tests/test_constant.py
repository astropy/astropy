# Licensed under a 3-clause BSD style license - see LICENSE.rst


from .. import Constant, si, cgs


def test_c():
    #c is an exactly defined constant, so it shouldn't be changing
    assert si.c == 2.99792458e8
    assert cgs.c == 2.99792458e10

    #check that the float-translation stuff works as expected
    assert si.c.real == si.c*1
    assert cgs.c.real == cgs.c-0

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
    #check that the value is fairly close to what it should be (not exactly
    #checking because this might get updated in the future)
    assert abs(si.h-6.626e-34) < 1e-38
    assert abs(cgs.h-6.626e-27) < 1e-31

    #make sure it has the necessary attributes and they're not blank
    assert si.h.error
    assert si.h.name
    assert si.h.origin
    assert si.h.units
    assert cgs.h.error
    assert cgs.h.name
    assert cgs.h.origin
    assert cgs.h.units

def test_info():
    # check the info functions can be successfully called and return
    # something.
    assert si.info()
    assert cgs.info()
