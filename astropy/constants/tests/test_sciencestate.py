# Licensed under a 3-clause BSD style license - see LICENSE.rst
import subprocess
import sys
from textwrap import dedent

import pytest

import astropy.constants as const
from astropy import astronomical_constants, physical_constants


def test_version_match():
    pversion = physical_constants.get()
    refpversion = const.h.__class__.__name__.lower()
    assert pversion == refpversion
    aversion = astronomical_constants.get()
    refaversion = const.M_sun.__class__.__name__.lower()
    assert aversion == refaversion


def test_previously_imported():
    with pytest.raises(RuntimeError):
        physical_constants.set("codata2018")

    with pytest.raises(RuntimeError):
        astronomical_constants.set("iau2015")


@pytest.mark.parametrize("version", sorted(set(physical_constants._versions.values())))
def test_physical_constants_versions(version):
    """Spot check that setting the different physical constants actually works.

    Here, checks the Rydberg constant as one that seems to change in every
    version (at least up to codata2022), and uses the electron charge to test
    that units derived from constants change as well.
    """
    cmd = dedent(f"""
    import astropy
    astropy.physical_constants.set({version!r})
    import astropy.constants as const
    from astropy.constants import {version} as {version}, codata2018 as codata2018
    assert const.Ryd.value == {version}.Ryd.value
    assert const.Ryd.uncertainty == {version}.Ryd.uncertainty
    if {version} != codata2018:
        assert {version}.Ryd.uncertainty != codata2018.Ryd.uncertainty
    import astropy.units as u
    assert u.eV.to(u.J) == {version}.e.value
    """)
    cp = subprocess.check_call([sys.executable, "-c", cmd])
    assert cp == 0
