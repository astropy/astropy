# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy

import numpy as np
import pytest

import astropy.constants as const
import astropy.units as u
from astropy import astronomical_constants, physical_constants
from astropy.constants import Constant, EMConstant


def test_c():
    from astropy.constants.codata2010 import c

    # c is an exactly defined constant, so it shouldn't be changing
    assert c.value == 2.99792458e8  # default is S.I.
    assert c.si.value == 2.99792458e8
    assert c.cgs.value == 2.99792458e10

    # make sure it has the necessary attributes and they're not blank
    assert c.uncertainty == 0  # c is a *defined* quantity
    assert c.name
    assert c.reference
    assert c.unit


# Check that if a constant is defined (uncertainty 0) in some version,
# it is also defined in the current, most up-to-date one.
@pytest.mark.parametrize(
    "version",
    sorted(
        set(physical_constants._versions.values())
        | set(astronomical_constants._versions.values())
    ),
)
def test_defined_constants_do_not_change(version):
    """Check that any defined constants do not change their value

    (the exception is mu0, which went from derived to measured in CODATA2018).
    """
    version_mod = getattr(const, version)
    defined_constants = [
        (name, constant)
        for (name, constant) in version_mod.__dict__.items()
        if isinstance(constant, Constant)
        and not isinstance(constant, EMConstant)
        and constant.uncertainty == 0
        and name != "mu0"
    ]
    for name, constant in defined_constants:
        current_constant = getattr(const, name)
        assert current_constant.uncertainty == 0
        if version == "iau2012" and name in ("pc", "kpc"):
            pytest.xfail(
                reason="different conversion from error-less au to pc and kpc."
            )
        assert current_constant.value == constant.value


def test_h():
    from astropy.constants import h as h_current
    from astropy.constants.codata2010 import h

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
    from astropy.constants.astropyconst13 import e as e_13

    # A test quantity
    E = u.Quantity(100.00000348276221, "V/m")

    # e.cgs is too ambiguous and should not work at all
    with pytest.raises(TypeError):
        e_13.cgs * E

    assert isinstance(e_13.si, u.Quantity)
    assert isinstance(e_13.gauss, u.Quantity)
    assert isinstance(e_13.esu, u.Quantity)

    assert e_13.gauss * E == u.Quantity(e_13.gauss.value * E.value, "Fr V/m")
    assert e_13.esu * E == u.Quantity(e_13.esu.value * E.value, "Fr V/m")


def test_g0():
    """Tests for #1263 demonstrating how g0 constant should behave."""
    from astropy.constants.astropyconst13 import g0

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
    assert g0.unit.physical_type == "acceleration"


def test_b_wien():
    """b_wien should give the correct peak wavelength for
    given blackbody temperature. The Sun is used in this test.

    """
    from astropy import units as u
    from astropy.constants.astropyconst13 import b_wien

    t = 5778 * u.K
    w = (b_wien / t).to(u.nm)
    assert round(w.value) == 502


def test_pc():
    """Parsec is defined to use small-angle limit per IAU 2015 Resolution B 2.
    iau2012 version still uses tan(parallax).
    """
    from astropy import units as u
    from astropy.constants import iau2012

    plx = np.radians(1 / 3600)
    assert np.allclose(
        u.pc.to("m") / iau2012.pc.si.value, np.tan(plx) / plx, rtol=1.0e-14, atol=0
    )


def test_masses():
    """Ensure mass values are set up correctly.
    https://github.com/astropy/astropy/issues/8920
    """
    from astropy.constants import astropyconst13, astropyconst20, astropyconst40

    ref_text = "Allen's Astrophysical Quantities 4th Ed."
    assert (
        astropyconst13.M_sun.reference == ref_text
        and astropyconst13.M_jup.reference == ref_text
        and astropyconst13.M_earth.reference == ref_text
    )

    ref_text = "IAU 2015 Resolution B 3 + CODATA 2014"
    assert (
        astropyconst20.M_sun.reference == ref_text
        and astropyconst20.M_jup.reference == ref_text
        and astropyconst20.M_earth.reference == ref_text
    )

    ref_text = "IAU 2015 Resolution B 3 + CODATA 2018"
    assert (
        astropyconst40.M_sun.reference == ref_text
        and astropyconst40.M_jup.reference == ref_text
        and astropyconst40.M_earth.reference == ref_text
    )


def test_unit():
    from astropy import units as u
    from astropy.constants import astropyconst13 as const

    for val in vars(const).values():
        if isinstance(val, Constant):
            # Getting the unit forces the unit parser to run.  Confirm
            # that none of the constants defined in astropy have
            # invalid unit.
            assert not isinstance(val.unit, u.UnrecognizedUnit)


def test_copy():
    from astropy import constants as const

    cc = copy.deepcopy(const.c)
    assert cc == const.c

    cc = copy.copy(const.c)
    assert cc == const.c


def test_view():
    """Check that Constant and Quantity views can be taken (#3537, #3538)."""
    from astropy.constants import c

    c2 = c.view(Constant)
    assert c2 == c
    assert c2.value == c.value
    # make sure it has the necessary attributes and they're not blank
    assert c2.uncertainty == 0  # c is a *defined* quantity
    assert c2.name == c.name
    assert c2.reference == c.reference
    assert c2.unit == c.unit

    q1 = c.view(u.Quantity)
    assert q1 == c
    assert q1.value == c.value
    assert type(q1) is u.Quantity
    assert not hasattr(q1, "reference")

    q2 = u.Quantity(c)
    assert q2 == c
    assert q2.value == c.value
    assert type(q2) is u.Quantity
    assert not hasattr(q2, "reference")

    c3 = u.Quantity(c, subok=True)
    assert c3 == c
    assert c3.value == c.value
    # make sure it has the necessary attributes and they're not blank
    assert c3.uncertainty == 0  # c is a *defined* quantity
    assert c3.name == c.name
    assert c3.reference == c.reference
    assert c3.unit == c.unit

    c4 = u.Quantity(c, subok=True, copy=False)
    assert c4 is c
