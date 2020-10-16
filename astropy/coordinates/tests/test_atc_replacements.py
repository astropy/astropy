# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test replacements for ERFA functions atciqz and aticq."""

from itertools import product

import pytest

from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time
from astropy import _erfa as erfa
from .utils import randomly_sample_sphere
from astropy.coordinates.builtin_frames.utils import get_jd12, atciqz, aticq

times = [Time("2014-06-25T00:00"), Time(["2014-06-25T00:00", "2014-09-24"])]
ra, dec, _ = randomly_sample_sphere(2)
positions = ((ra[0], dec[0]), (ra, dec))
spacetimes = product(times, positions)


@pytest.mark.parametrize('st', spacetimes)
def test_atciqz_aticq(st):
    """Check replacements against erfa versions for consistency."""
    t, pos = st
    jd1, jd2 = get_jd12(t, 'tdb')
    astrom, _ = erfa.apci13(jd1, jd2)

    ra, dec = pos
    ra = ra.value
    dec = dec.value
    assert_allclose(erfa.atciqz(ra, dec, astrom), atciqz(ra, dec, astrom))
    assert_allclose(erfa.aticq(ra, dec, astrom), aticq(ra, dec, astrom))
