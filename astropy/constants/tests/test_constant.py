# Licensed under a 3-clause BSD style license - see LICENSE.rst


from .. import Constant
from ...units import Quantity as Q
from ...tests.helper import pytest


def test_c():

    from .. import c

    # c is an exactly defined constant, so it shouldn't be changing
    assert c.value == 2.99792458e8  # default is S.I.
    assert c.si.value == 2.99792458e8
    assert c.cgs.value == 2.99792458e10

    # make sure it has the necessary attributes and they're not blank
    assert c.uncertainty == 0  # c is a *defined* quantity
    assert c.name
    assert c.reference
    assert c.unit


def test_h():

    from .. import h

    # check that the value is fairly close to what it should be (not exactly
    # checking because this might get updated in the future)
    assert abs(h.value - 6.626e-34) < 1e-38
    assert abs(h.si.value - 6.626e-34) < 1e-38
    assert abs(h.cgs.value - 6.626e-27) < 1e-31

    # make sure it has the necessary attributes and they're not blank
    assert h.uncertainty
    assert h.name
    assert h.reference
    assert h.unit


def test_e():
    """Tests for #572 demonstrating how EM constants should behave."""

    from .. import e

    # A test quantity
    E = Q(100, 'V/m')

    # Without specifying a system e should not combine with other quantities
    pytest.raises(TypeError, lambda: e * E)

    # e.cgs is too ambiguous and should not work at all
    pytest.raises(TypeError, lambda: e.cgs * E)

    assert isinstance(e.si, Q)
    assert isinstance(e.gauss, Q)
    assert isinstance(e.esu, Q)

    assert e.si * E == Q(100, 'eV/m')
    assert e.gauss * E == Q(100, 'Fr V/m')
    assert e.esu * E == Q(100, 'Fr V/m')


def test_unit():

    from ... import units as u

    from ... import constants as const

    for key, val in vars(const).iteritems():
        if isinstance(val, Constant):
            # Getting the unit forces the unit parser to run.  Confirm
            # that none of the constants defined in astropy have
            # invalid unit.
            assert not isinstance(val.unit, u.UnrecognizedUnit)
