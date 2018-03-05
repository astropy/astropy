# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

asdf = pytest.importorskip('asdf', minversion='2.0.0')
from asdf.tests.helpers import assert_roundtrip_tree

from astropy import units
from astropy.coordinates import ICRS, Longitude, Latitude, Angle

from ....extension import AstropyExtension


def test_icrs_basic(tmpdir):
    wrap_angle = Angle(1.5, unit=units.rad)
    ra = Longitude(25, unit=units.deg, wrap_angle=wrap_angle)
    dec = Latitude(45, unit=units.deg)

    tree = {'coord': ICRS(ra=ra, dec=dec)}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())


@pytest.mark.xfail(
    reason="Compound ICRS coordinates have not been implemented for ASDF yet")
def test_icrs_compound(tmpdir):

    icrs = ICRS(ra=[0, 1, 2]*units.deg, dec=[3, 4, 5]*units.deg)

    tree = {'coord': icrs}

    assert_roundtrip_tree(tree, tmpdir, extensions=AstropyExtension())
