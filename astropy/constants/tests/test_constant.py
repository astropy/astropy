# Licensed under a 3-clause BSD style license - see LICENSE.rst


from .. import Constant, si, cgs


def test_c():
    #c is an exactly defined constant, so it shouldn't be changing
    assert si.c.value == 2.99792458e8
    assert cgs.c.value == 2.99792458e10

    #make sure it has the necessary attributes and they're not blank
    assert si.c.uncertainty == 0  # c is a *defined* quantity
    assert si.c.name
    assert si.c.reference
    assert si.c.unit
    assert cgs.c.uncertainty == 0  # c is a *defined* quantity
    assert cgs.c.name
    assert cgs.c.reference
    assert cgs.c.unit


def test_h():
    #check that the value is fairly close to what it should be (not exactly
    #checking because this might get updated in the future)
    assert abs(si.h.value - 6.626e-34) < 1e-38
    assert abs(cgs.h.value - 6.626e-27) < 1e-31

    #make sure it has the necessary attributes and they're not blank
    assert si.h.uncertainty
    assert si.h.name
    assert si.h.reference
    assert si.h.unit
    assert cgs.h.uncertainty
    assert cgs.h.name
    assert cgs.h.reference
    assert cgs.h.unit


def test_unit():

    from ... import units as u

    for key, val in si.__dict__.items() + cgs.__dict__.items():
        if isinstance(val, Constant):
            # Getting the unit forces the unit parser to run.  Confirm
            # that none of the constants defined in astropy have
            # invalid unit.
            assert not isinstance(val.unit, u.UnrecognizedUnit)
