# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Property-based tests for unit decomposition.
"""

from hypothesis import given, strategies as st

import astropy.units as u

UNITS = sorted(u.get_current_unit_registry().all_units, key=repr)


@given(st.lists(st.sampled_from(UNITS), min_size=2))
def test_intermediate_decomposition_is_no_op(lst):
    # This test checks that when multiplying togther any sequence of units,
    # it makes no difference whether we decompose between every step or only
    # at the end.  This is true for rational numbers, but not for floats!
    decomp = compound = u.dimensionless_unscaled
    for unit in lst:
        decomp = (decomp * unit).decompose()
        compound *= unit
        assert decomp == compound.decompose(), (decomp, compound)
