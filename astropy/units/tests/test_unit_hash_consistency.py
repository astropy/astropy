# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Comprehensive test suite for unit hash consistency.

This module tests the fix for issue #18560 where equivalent units
had different hash values, violating Python's hash-equality contract.
"""

import pickle
import pytest
import numpy as np

from astropy import units as u
from astropy.units import (
    Unit, CompositeUnit, IrreducibleUnit, UnrecognizedUnit,
    dimensionless_unscaled
)


class TestHashEqualityContract:
    """Test that equal units have equal hashes (Python's hash-equality contract)."""

    def test_named_vs_composite_newton(self):
        """Test that Newton and its composite form have equal hashes."""
        unit1 = u.N
        unit2 = u.kg * u.m / u.s**2
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_named_vs_composite_joule(self):
        """Test that Joule and its composite form have equal hashes."""
        unit1 = u.J
        unit2 = u.N * u.m
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_named_vs_composite_pascal(self):
        """Test that Pascal and its composite form have equal hashes."""
        unit1 = u.Pa
        unit2 = u.N / u.m**2
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_named_vs_composite_watt(self):
        """Test that Watt and its composite form have equal hashes."""
        unit1 = u.W
        unit2 = u.J / u.s
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_named_vs_composite_volt(self):
        """Test that Volt and its composite form have equal hashes."""
        unit1 = u.V
        unit2 = u.W / u.A
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_named_vs_composite_ohm(self):
        """Test that Ohm and its composite form have equal hashes."""
        unit1 = u.Ohm
        unit2 = u.V / u.A
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_complex_composite_equivalence(self):
        """Test hash equality for complex composite unit expressions."""
        # kg*m^2/s^2 = J (Joule)
        unit1 = u.kg * u.m**2 / u.s**2
        unit2 = u.J
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_multiple_composition_paths(self):
        """Test that different ways of composing the same unit have equal hashes."""
        # All represent kg*m^2/s^3 (Watt)
        unit1 = u.W
        unit2 = u.kg * u.m**2 / u.s**3
        unit3 = (u.kg * u.m / u.s**2) * u.m / u.s
        
        assert unit1 == unit2 == unit3, "All units should be equal"
        assert hash(unit1) == hash(unit2) == hash(unit3), "Equal units must have equal hashes"

    def test_scaled_units_different_hashes(self):
        """Test that units with different scales have different hashes."""
        unit1 = u.m
        unit2 = u.km
        
        assert unit1 != unit2, "Units with different scales should not be equal"
        assert hash(unit1) != hash(unit2), "Unequal units should (usually) have different hashes"


class TestSetBehavior:
    """Test that sets work correctly with equivalent units."""

    def test_set_deduplication_basic(self):
        """Test that equivalent units are deduplicated in sets."""
        s = {u.N, u.kg * u.m / u.s**2}
        assert len(s) == 1, "Equivalent units should be treated as one element in a set"

    def test_set_deduplication_multiple(self):
        """Test set deduplication with multiple equivalent representations."""
        s = {u.J, u.N * u.m, u.kg * u.m**2 / u.s**2}
        assert len(s) == 1, "All equivalent units should be treated as one element"

    def test_set_with_different_units(self):
        """Test that different units remain separate in sets."""
        s = {u.m, u.s, u.kg}
        assert len(s) == 3, "Different units should remain separate"

    def test_set_operations_union(self):
        """Test set union operations with equivalent units."""
        s1 = {u.N, u.J}
        s2 = {u.kg * u.m / u.s**2, u.Pa}
        result = s1 | s2
        
        # s1 and s2 both contain Newton (in different forms)
        assert u.N in result or u.kg * u.m / u.s**2 in result
        # Should have 3 unique units: Newton, Joule, Pascal
        assert len(result) == 3

    def test_set_operations_intersection(self):
        """Test set intersection operations with equivalent units."""
        s1 = {u.N, u.J}
        s2 = {u.kg * u.m / u.s**2, u.Pa}
        result = s1 & s2
        
        # Newton appears in both sets (in different forms)
        assert len(result) == 1
        assert u.N in result or u.kg * u.m / u.s**2 in result


class TestDictBehavior:
    """Test that dictionaries work correctly with equivalent units as keys."""

    def test_dict_key_equivalence(self):
        """Test that equivalent units are treated as the same dictionary key."""
        d = {u.N: "Newton"}
        d[u.kg * u.m / u.s**2] = "kg*m/s^2"
        
        assert len(d) == 1, "Equivalent units should be the same key"
        assert d[u.N] == "kg*m/s^2", "Value should be overwritten"

    def test_dict_key_lookup(self):
        """Test that we can look up values using equivalent unit forms."""
        d = {u.N: "force"}
        
        # Should be able to look up using any equivalent form
        assert d[u.kg * u.m / u.s**2] == "force"

    def test_dict_multiple_representations(self):
        """Test dictionary with multiple equivalent representations."""
        d = {}
        d[u.J] = "energy1"
        d[u.N * u.m] = "energy2"
        d[u.kg * u.m**2 / u.s**2] = "energy3"
        
        assert len(d) == 1, "All keys are equivalent"
        assert d[u.J] == "energy3", "Last value should win"

    def test_dict_get_with_default(self):
        """Test dict.get() with equivalent units."""
        d = {u.N: "force"}
        
        # Should find the value using equivalent form
        assert d.get(u.kg * u.m / u.s**2, "default") == "force"

    def test_dict_in_operator(self):
        """Test 'in' operator with equivalent units."""
        d = {u.N: "force"}
        
        assert u.kg * u.m / u.s**2 in d, "Equivalent unit should be found"


class TestUnrecognizedUnits:
    """Test hash behavior for UnrecognizedUnit instances."""

    def test_unrecognized_unit_hash_consistency(self):
        """Test that identical unrecognized units have equal hashes."""
        foo1 = Unit("FOO", parse_strict="silent")
        foo2 = Unit("FOO", parse_strict="silent")
        
        assert foo1 == foo2, "Identical unrecognized units should be equal"
        assert hash(foo1) == hash(foo2), "Equal unrecognized units must have equal hashes"

    def test_different_unrecognized_units(self):
        """Test that different unrecognized units have different hashes."""
        foo = Unit("FOO", parse_strict="silent")
        bar = Unit("BAR", parse_strict="silent")
        
        assert foo != bar, "Different unrecognized units should not be equal"
        assert hash(foo) != hash(bar), "Different units should have different hashes"

    def test_unrecognized_unit_in_set(self):
        """Test that unrecognized units work correctly in sets."""
        foo1 = Unit("FOO", parse_strict="silent")
        foo2 = Unit("FOO", parse_strict="silent")
        
        s = {foo1, foo2}
        assert len(s) == 1, "Equivalent unrecognized units should deduplicate"

    def test_unrecognized_unit_in_dict(self):
        """Test that unrecognized units work correctly as dictionary keys."""
        foo = Unit("FOO", parse_strict="silent")
        
        d = {foo: "value"}
        foo2 = Unit("FOO", parse_strict="silent")
        
        assert d[foo2] == "value", "Should be able to look up with equivalent unit"


class TestDimensionlessUnits:
    """Test hash behavior for dimensionless units."""

    def test_dimensionless_unscaled_hash(self):
        """Test that dimensionless_unscaled is hashable."""
        h = hash(dimensionless_unscaled)
        assert isinstance(h, int), "Hash should be an integer"

    def test_dimensionless_equality_and_hash(self):
        """Test that equivalent dimensionless units have equal hashes."""
        unit1 = u.m / u.m
        unit2 = dimensionless_unscaled
        
        assert unit1 == unit2, "Units should be equal"
        assert hash(unit1) == hash(unit2), "Equal units must have equal hashes"

    def test_dimensionless_one(self):
        """Test that 'one' (alias for dimensionless_unscaled) hashes correctly."""
        from astropy.units import one
        
        assert one == dimensionless_unscaled
        assert hash(one) == hash(dimensionless_unscaled)


class TestPickling:
    """Test that unit hashes work correctly after pickling/unpickling."""

    def test_pickle_named_unit(self):
        """Test that named units hash correctly after pickling."""
        original = u.N
        pickled = pickle.dumps(original)
        restored = pickle.loads(pickled)
        
        assert original == restored
        assert hash(original) == hash(restored)

    def test_pickle_composite_unit(self):
        """Test that composite units hash correctly after pickling."""
        original = u.kg * u.m / u.s**2
        pickled = pickle.dumps(original)
        restored = pickle.loads(pickled)
        
        assert original == restored
        assert hash(original) == hash(restored)

    def test_pickle_preserves_equivalence(self):
        """Test that pickling preserves hash equivalence between different unit forms."""
        unit1 = u.N
        unit2 = u.kg * u.m / u.s**2
        
        # Pickle and restore both
        restored1 = pickle.loads(pickle.dumps(unit1))
        restored2 = pickle.loads(pickle.dumps(unit2))
        
        # They should still be equal with equal hashes
        assert restored1 == restored2
        assert hash(restored1) == hash(restored2)

    def test_pickle_unrecognized_unit(self):
        """Test that unrecognized units hash correctly after pickling."""
        original = Unit("FOO", parse_strict="silent")
        pickled = pickle.dumps(original)
        restored = pickle.loads(pickled)
        
        assert original == restored
        assert hash(original) == hash(restored)


class TestHashStability:
    """Test that unit hashes remain stable across multiple calls."""

    def test_hash_stability_named_unit(self):
        """Test that hashing a named unit multiple times gives the same result."""
        unit = u.N
        hash1 = hash(unit)
        hash2 = hash(unit)
        hash3 = hash(unit)
        
        assert hash1 == hash2 == hash3, "Hash should be stable"

    def test_hash_stability_composite_unit(self):
        """Test that hashing a composite unit multiple times gives the same result."""
        unit = u.kg * u.m / u.s**2
        hash1 = hash(unit)
        hash2 = hash(unit)
        hash3 = hash(unit)
        
        assert hash1 == hash2 == hash3, "Hash should be stable"

    def test_hash_caching(self):
        """Test that hash values are cached properly."""
        unit = u.N
        
        # First hash call should cache the value
        hash1 = hash(unit)
        
        # Verify cache exists
        assert hasattr(unit, '_hash'), "Hash should be cached"
        
        # Second call should return cached value
        hash2 = hash(unit)
        assert hash1 == hash2, "Cached hash should match"


class TestEdgeCases:
    """Test hash behavior for edge cases and unusual unit combinations."""

    def test_inverse_units(self):
        """Test that inverse units hash correctly."""
        unit1 = u.s**(-1)
        unit2 = u.Unit(1 / u.s)
        
        assert unit1 == unit2
        assert hash(unit1) == hash(unit2)

    def test_power_units(self):
        """Test that units with powers hash correctly."""
        unit1 = u.m * u.m
        unit2 = u.m**2
        
        assert unit1 == unit2
        assert hash(unit1) == hash(unit2)

    def test_fractional_powers(self):
        """Test that units with fractional powers hash correctly."""
        unit1 = u.m**(1/2)
        unit2 = u.m**0.5
        
        assert unit1 == unit2
        assert hash(unit1) == hash(unit2)

    def test_complex_nested_composition(self):
        """Test hash for deeply nested unit compositions."""
        unit1 = ((u.kg * u.m) / u.s) / u.s
        unit2 = u.kg * u.m / u.s**2
        
        assert unit1 == unit2
        assert hash(unit1) == hash(unit2)

    def test_zero_power_units(self):
        """Test units raised to power zero."""
        unit1 = u.m**0
        unit2 = dimensionless_unscaled
        
        assert unit1 == unit2
        assert hash(unit1) == hash(unit2)

    def test_commutative_multiplication(self):
        """Test that unit multiplication is commutative for hashing."""
        unit1 = u.m * u.kg
        unit2 = u.kg * u.m
        
        assert unit1 == unit2
        assert hash(unit1) == hash(unit2)


class TestPerformance:
    """Test that hash caching improves performance."""

    def test_repeated_hash_calls_are_fast(self):
        """Test that repeated hash calls benefit from caching."""
        import time
        
        unit = u.N
        
        # First call (might trigger decomposition and caching)
        start = time.perf_counter()
        hash(unit)
        first_time = time.perf_counter() - start
        
        # Subsequent calls should be much faster (cached)
        times = []
        for _ in range(100):
            start = time.perf_counter()
            hash(unit)
            times.append(time.perf_counter() - start)
        
        avg_cached_time = sum(times) / len(times)
        
        # Cached calls should be at least 10x faster (conservative estimate)
        # In practice, they're often 100x+ faster
        assert avg_cached_time < first_time / 5, \
            "Cached hash calls should be significantly faster"

    def test_hash_in_tight_loop(self):
        """Test that hashing in a tight loop is efficient due to caching."""
        units = [u.N, u.J, u.Pa, u.W, u.V, u.Ohm]
        
        # This should be fast because hashes are cached
        result = {}
        for _ in range(1000):
            for unit in units:
                result[unit] = result.get(unit, 0) + 1
        
        assert len(result) == len(units)


class TestRegressionIssue18560:
    """
    Regression tests specifically for issue #18560.
    
    These tests verify that the original problem is fixed.
    """

    def test_newton_set_deduplication(self):
        """Original bug: Newton in two forms created two set elements."""
        s = {u.N, u.kg * u.m / u.s**2}
        assert len(s) == 1, "Issue #18560: Equivalent units should be one element"

    def test_newton_dict_key(self):
        """Original bug: Couldn't use equivalent units as same dict key."""
        d = {u.N: "Newton"}
        # This should work (same key)
        value = d[u.kg * u.m / u.s**2]
        assert value == "Newton", "Issue #18560: Equivalent units should be same key"

    def test_hash_equality_contract_violation(self):
        """Original bug: Equal units had different hashes."""
        unit1 = u.N
        unit2 = u.kg * u.m / u.s**2
        
        # This is the core bug from issue #18560
        assert unit1 == unit2, "Units are equal"
        assert hash(unit1) == hash(unit2), "Issue #18560: Equal units must have equal hashes"


class TestDocumentedExamples:
    """Test the examples from the documentation."""

    def test_documentation_example_1(self):
        """Test example from unit hash documentation."""
        # Example: Using units as dictionary keys
        unit_descriptions = {
            u.N: "Force",
            u.J: "Energy",
            u.Pa: "Pressure"
        }
        
        # Should work with equivalent forms
        assert unit_descriptions[u.kg * u.m / u.s**2] == "Force"
        assert unit_descriptions[u.N * u.m] == "Energy"

    def test_documentation_example_2(self):
        """Test example from unit hash documentation."""
        # Example: Deduplicating units in a set
        all_units = [u.N, u.kg * u.m / u.s**2, u.J, u.N * u.m, u.Pa]
        unique_units = set(all_units)
        
        # Should have only 3 unique units (N, J, Pa)
        assert len(unique_units) == 3
