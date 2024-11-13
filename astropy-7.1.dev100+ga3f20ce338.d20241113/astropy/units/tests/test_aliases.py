# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test setting and adding unit aliases."""

import pytest

import astropy.units as u

trials = [
    ({"Angstroms": u.AA}, "Angstroms", u.AA),
    ({"counts": u.count}, "counts/s", u.count / u.s),
    (
        {"ergs": u.erg, "Angstroms": u.AA},
        "ergs/(s cm**2 Angstroms)",
        u.erg / (u.s * u.cm**2 * u.AA),
    ),
]


class TestAliases:
    def teardown_method(self):
        u.set_enabled_aliases({})

    def teardown_class(self):
        assert u.get_current_unit_registry().aliases == {}

    @pytest.mark.parametrize("format_", [None, "fits", "ogip", "vounit", "cds"])
    @pytest.mark.parametrize("aliases,bad,unit", trials)
    def test_set_enabled_aliases_context_manager(self, aliases, bad, unit, format_):
        if format_ == "cds":
            bad = bad.replace(" ", ".").replace("**", "")

        with u.set_enabled_aliases(aliases):
            assert u.get_current_unit_registry().aliases == aliases
            assert u.Unit(bad) == unit

        assert u.get_current_unit_registry().aliases == {}
        with pytest.raises(ValueError):
            u.Unit(bad)

    @pytest.mark.parametrize("aliases,bad,unit", trials)
    def test_add_enabled_aliases_context_manager(self, aliases, bad, unit):
        with u.add_enabled_aliases(aliases):
            assert u.get_current_unit_registry().aliases == aliases
            assert u.Unit(bad) == unit

        assert u.get_current_unit_registry().aliases == {}
        with pytest.raises(ValueError):
            u.Unit(bad)

    def test_set_enabled_aliases(self):
        for i, (aliases, bad, unit) in enumerate(trials):
            u.set_enabled_aliases(aliases)

            assert u.get_current_unit_registry().aliases == aliases

            assert u.Unit(bad) == unit

            for _, bad2, unit2 in trials:
                if bad2 == bad or bad2 in aliases:
                    assert u.Unit(bad2) == unit2
                else:
                    with pytest.raises(ValueError):
                        u.Unit(bad2)

    def test_add_enabled_aliases(self):
        expected_aliases = {}
        for i, (aliases, bad, unit) in enumerate(trials):
            u.add_enabled_aliases(aliases)

            expected_aliases.update(aliases)
            assert u.get_current_unit_registry().aliases == expected_aliases

            assert u.Unit(bad) == unit

            for j, (_, bad2, unit2) in enumerate(trials):
                if j <= i:
                    assert u.Unit(bad2) == unit2
                else:
                    with pytest.raises(ValueError):
                        u.Unit(bad2)

    def test_cannot_alias_existing_unit(self):
        with pytest.raises(ValueError, match="already means"):
            u.set_enabled_aliases({"pct": u.Unit(1e-12 * u.count)})

    def test_cannot_alias_existing_alias_to_another_unit(self):
        u.set_enabled_aliases({"counts": u.count})
        with pytest.raises(ValueError, match="already is an alias"):
            u.add_enabled_aliases({"counts": u.adu})
