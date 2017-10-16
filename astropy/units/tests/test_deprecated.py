# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


"""Regression tests for deprecated units or those that are "soft" deprecated
because they are required for VOUnit support but are not in common use."""

import pytest

from .. import deprecated, required_by_vounit
from ... import units as u


def test_emu():
    with pytest.raises(AttributeError):
        u.emu

    assert u.Bi.to(deprecated.emu, 1) == 1

    with deprecated.enable():
        assert u.Bi.compose()[0] == deprecated.emu

    assert u.Bi.compose()[0] == u.Bi

    # test that the earth/jupiter mass/rad are also in the deprecated bunch
    for body in ('earth', 'jupiter'):
        for phystype in ('Mass', 'Rad'):
            # only test a couple prefixes to same time
            for prefix in ('n', 'y'):
                namewoprefix = body + phystype
                unitname = prefix + namewoprefix

                with pytest.raises(AttributeError):
                    getattr(u, unitname)

                assert (getattr(deprecated, unitname).represents.bases[0] ==
                        getattr(u, namewoprefix))


def test_required_by_vounit():
    # The tests below could be replicated with all the various prefixes, but it
    # seems unnecessary because they all come as a set.  So we only use nano for
    # the purposes of this test.

    with pytest.raises(AttributeError):
        # nano-solar mass/rad/lum shouldn't be in the base unit namespace
        u.nsolMass
        u.nsolRad
        u.nsolLum

    # but they should be enabled by default via required_by_vounit, to allow
    # the Unit constructor to accept them
    assert u.Unit('nsolMass') == required_by_vounit.nsolMass
    assert u.Unit('nsolRad') == required_by_vounit.nsolRad
    assert u.Unit('nsolLum') == required_by_vounit.nsolLum

    # but because they are prefixes, they shouldn't be in find_equivalent_units
    assert required_by_vounit.nsolMass not in u.solMass.find_equivalent_units()
    assert required_by_vounit.nsolRad not in u.solRad.find_equivalent_units()
    assert required_by_vounit.nsolLum not in u.solLum.find_equivalent_units()
