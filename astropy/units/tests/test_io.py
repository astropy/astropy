# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import ast
import json
import re

import numpy as np
import pytest

import astropy.units as u
from astropy.io.misc.json import JSONExtendedEncoder, JSONExtendedDecoder
from astropy.io.misc.json.tests.test_core import JSONExtendedTestBase


class TestJSONExtendedUnits(JSONExtendedTestBase):
    """Tests for serializing builtins with extended JSON encoders and decoders."""

    def test_unit_simple(self):
        """Test round-tripping `astropy.units.Unit`."""
        obj = u.Unit("km")

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match=re.escape("Object of type PrefixUnit is not JSON serializable")):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "astropy.units.core.PrefixUnit"
        assert d["value"] == "km"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, u.UnitBase)
        assert out == obj

    @pytest.mark.skip("TODO!")
    def test_unit_structured(self):
        """Test round-tripping structured `astropy.units.Unit`."""
        obj = u.Unit(("km", "eV"))

    def test_quantity_simple(self):
        """Test round-tripping `astropy.units.Quantity`."""
        obj = u.Quantity([3, 4], dtype=float, unit=u.km)

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match=re.escape("Object of type Quantity is not JSON serializable")):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "astropy.units.quantity.Quantity"
        assert d["value"] == [3.0, 4.0]
        assert d["unit"] == "km"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, u.Quantity)
        assert np.array_equal(out, obj)

    @pytest.mark.skip("TODO!")
    def test_quantity_structured(self):
        """Test round-tripping structured `astropy.units.Quantity`."""
        dt = np.unit([("f1", np.int64), ("f2", np.float16)])
        obj = u.Quantity((1, 3.0), dtype=dt, unit=u.Unit("(u.km, u.eV)"))
