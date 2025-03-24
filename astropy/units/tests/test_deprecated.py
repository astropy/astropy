# Licensed under a 3-clause BSD style license - see LICENSE.rst


import pytest

from astropy import units as u
from astropy.units import deprecated


def test_emu():
    with pytest.raises(AttributeError):
        u.emu

    assert u.Bi.to(deprecated.emu, 1) == 1

    with deprecated.enable():
        assert u.Bi.compose()[0] == deprecated.emu

    assert u.Bi.compose()[0] == u.Bi

    # test that the earth/jupiter mass/rad are also in the deprecated bunch
    for body in ("earth", "jupiter"):
        for phystype in ("Mass", "Rad"):
            # only test a couple prefixes to same time
            for prefix in ("n", "y"):
                namewoprefix = body + phystype
                unitname = prefix + namewoprefix

                with pytest.raises(AttributeError):
                    getattr(u, unitname)

                assert getattr(deprecated, unitname).represents.bases[0] == getattr(
                    u, namewoprefix
                )
