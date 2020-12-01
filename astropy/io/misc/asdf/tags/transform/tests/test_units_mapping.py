# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

asdf = pytest.importorskip("asdf")

from asdf.tests.helpers import assert_roundtrip_tree

from astropy.modeling.models import UnitsMapping
from astropy import units as u


def assert_model_roundtrip(model, tmpdir):
    tree = {"model": model}
    assert_roundtrip_tree(tree, tmpdir, tree_match_func=assert_models_equal)


def assert_models_equal(a, b):
    assert a.name == b.name
    assert a.inputs == b.inputs
    assert a.input_units == b.input_units
    assert a.outputs == b.outputs
    assert a.mapping == b.mapping
    assert a.input_units_allow_dimensionless == b.input_units_allow_dimensionless

    for i in a.inputs:
        if a.input_units_equivalencies is None:
            a_equiv = None
        else:
            a_equiv = a.input_units_equivalencies.get(i)

        if b.input_units_equivalencies is None:
            b_equiv = None
        else:
            b_equiv = b.input_units_equivalencies.get(i, None)

        assert a_equiv == b_equiv


def test_basic(tmpdir):
    m = UnitsMapping(((u.m, u.dimensionless_unscaled),))
    assert_model_roundtrip(m, tmpdir)


def test_remove_units(tmpdir):
    m = UnitsMapping(((u.m, None),))
    assert_model_roundtrip(m, tmpdir)


def test_accept_any_units(tmpdir):
    m = UnitsMapping(((None, u.m),))
    assert_model_roundtrip(m, tmpdir)


def test_with_equivalencies(tmpdir):
    m = UnitsMapping(((u.m, u.dimensionless_unscaled),), input_units_equivalencies={"x": u.equivalencies.spectral()})
    assert_model_roundtrip(m, tmpdir)


def test_with_allow_dimensionless(tmpdir):
    m = UnitsMapping(((u.m, u.dimensionless_unscaled), (u.s, u.Hz)), input_units_allow_dimensionless=True)
    assert_model_roundtrip(m, tmpdir)

    m = UnitsMapping(((u.m, u.dimensionless_unscaled), (u.s, u.Hz)), input_units_allow_dimensionless={"x0": True, "x1": False})
    assert_model_roundtrip(m, tmpdir)
