# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test setting and adding unit aliases."""
import pytest

import astropy.units as u


trials = [
    ({"counts": u.count}, "10**(-6) counts/s", "10**(-6) count/s"),
    ({"Angstroms": u.AA}, "Angstroms", "Angstrom"),
    ({"ergs": u.erg, "Angstroms": u.AA},
     "10^-17 ergs/(s.cm^2.Angstroms)", "10^-17 erg/(s.cm^2.Angstrom)")]


class TestAliases:
    def teardown_method(self):
        u.set_enabled_aliases({})

    def test_set_enabled_aliases(self):
        for i, (aliases, bad, good) in enumerate(trials):
            u.set_enabled_aliases(aliases)
            assert u.Unit(bad) == u.Unit(good)

            for _, bad2, good2 in trials:
                if bad2 == bad or bad2 in aliases:
                    assert u.Unit(bad2) == u.Unit(good2)
                else:
                    with pytest.raises(ValueError):
                        u.Unit(bad2)

    def test_add_enabled_aliases(self):
        for i, (aliases, bad, good) in enumerate(trials):
            u.add_enabled_aliases(aliases)
            assert u.Unit(bad) == u.Unit(good)

            for j, (_, bad2, good2) in enumerate(trials):
                if j <= i:
                    assert u.Unit(bad2) == u.Unit(good2)
                else:
                    with pytest.raises(ValueError):
                        u.Unit(bad2)

    @pytest.mark.parametrize('aliases,bad,good', trials)
    def test_set_enabled_aliases_context_manager(self, aliases, bad, good):
        with u.set_enabled_aliases(aliases):
            assert u.Unit(bad) == u.Unit(good)

        with pytest.raises(ValueError):
            u.Unit(bad)

    @pytest.mark.parametrize('aliases,bad,good', trials)
    def test_add_enabled_aliases_context_manager(self, aliases, bad, good):
        with u.add_enabled_aliases(aliases):
            assert u.Unit(bad) == u.Unit(good)

        with pytest.raises(ValueError):
            u.Unit(bad)

    def test_cannot_alias_existing_unit(self):
        with pytest.raises(ValueError, match='already means'):
            u.set_enabled_aliases({'pct': u.Unit(1e-12*u.count)})

    def test_cannot_alias_existing_alias_to_another_unit(self):
        u.set_enabled_aliases({'counts': u.count})
        with pytest.raises(ValueError, match='already means'):
            u.add_enabled_aliases({'counts': u.adu})
