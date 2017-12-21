# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy

import pytest

from .. import Constant
from ...units import Quantity as Q


def test_c():

    from ..codata2010 import c

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

    from ..codata2010 import h
    from .. import h as h_current

    # check that the value is the CODATA2010 value
    assert abs(h.value - 6.62606957e-34) < 1e-43
    assert abs(h.si.value - 6.62606957e-34) < 1e-43
    assert abs(h.cgs.value - 6.62606957e-27) < 1e-36

    # Check it is different than the current value
    assert abs(h.value - h_current.value) > 4e-42

    # make sure it has the necessary attributes and they're not blank
    assert h.uncertainty
    assert h.name
    assert h.reference
    assert h.unit


def test_e():

    from ..astropyconst13 import e

    # A test quantity
    E = Q(100.00000348276221, 'V/m')

    # e.cgs is too ambiguous and should not work at all
    with pytest.raises(TypeError):
        e.cgs * E

    assert isinstance(e.si, Q)
    assert isinstance(e.gauss, Q)
    assert isinstance(e.esu, Q)

    assert e.si * E == Q(100, 'eV/m')
    assert e.gauss * E == Q(e.gauss.value * E.value, 'Fr V/m')
    assert e.esu * E == Q(e.esu.value * E.value, 'Fr V/m')


def test_g0():
    """Tests for #1263 demonstrating how g0 constant should behave."""
    from ..astropyconst13 import g0

    # g0 is an exactly defined constant, so it shouldn't be changing
    assert g0.value == 9.80665  # default is S.I.
    assert g0.si.value == 9.80665
    assert g0.cgs.value == 9.80665e2

    # make sure it has the necessary attributes and they're not blank
    assert g0.uncertainty == 0  # g0 is a *defined* quantity
    assert g0.name
    assert g0.reference
    assert g0.unit

    # Check that its unit have the correct physical type
    assert g0.unit.physical_type == 'acceleration'


def test_b_wien():
    """b_wien should give the correct peak wavelength for
    given blackbody temperature. The Sun is used in this test.

    """
    from ..astropyconst13 import b_wien
    from ... import units as u
    t = 5778 * u.K
    w = (b_wien / t).to(u.nm)
    assert round(w.value) == 502


def test_unit():

    from ... import units as u

    from .. import astropyconst13 as const

    for key, val in vars(const).items():
        if isinstance(val, Constant):
            # Getting the unit forces the unit parser to run.  Confirm
            # that none of the constants defined in astropy have
            # invalid unit.
            assert not isinstance(val.unit, u.UnrecognizedUnit)


def test_copy():
    from ... import constants as const
    cc = copy.deepcopy(const.c)
    assert cc == const.c

    cc = copy.copy(const.c)
    assert cc == const.c


def test_view():
    """Check that Constant and Quantity views can be taken (#3537, #3538)."""
    from .. import c
    c2 = c.view(Constant)
    assert c2 == c
    assert c2.value == c.value
    # make sure it has the necessary attributes and they're not blank
    assert c2.uncertainty == 0  # c is a *defined* quantity
    assert c2.name == c.name
    assert c2.reference == c.reference
    assert c2.unit == c.unit

    q1 = c.view(Q)
    assert q1 == c
    assert q1.value == c.value
    assert type(q1) is Q
    assert not hasattr(q1, 'reference')

    q2 = Q(c)
    assert q2 == c
    assert q2.value == c.value
    assert type(q2) is Q
    assert not hasattr(q2, 'reference')

    c3 = Q(c, subok=True)
    assert c3 == c
    assert c3.value == c.value
    # make sure it has the necessary attributes and they're not blank
    assert c3.uncertainty == 0  # c is a *defined* quantity
    assert c3.name == c.name
    assert c3.reference == c.reference
    assert c3.unit == c.unit

    c4 = Q(c, subok=True, copy=False)
    assert c4 is c


def test_context_manager():
    from ... import constants as const

    with const.set_enabled_constants('astropyconst13'):
        assert const.h.value == 6.62606957e-34  # CODATA2010

    assert const.h.value == 6.626070040e-34  # CODATA2014

    with pytest.raises(ValueError):
        with const.set_enabled_constants('notreal'):
            const.h
