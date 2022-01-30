# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import ast
import json
import re

import numpy as np
import pytest

import astropy.units as u
import astropy.coordinates as coord
from astropy.io.misc.json import JSONExtendedEncoder, JSONExtendedDecoder
from astropy.io.misc.json.tests.test_core import JSONExtendedTestBase


class TestJSONExtendedUnits(JSONExtendedTestBase):
    """Tests for serializing builtins with extended JSON encoders and decoders."""

    def test_longitude(self):
        """Test round-tripping `astropy.coordinates.Longitude`."""
        obj = coord.Longitude([3, 4], dtype=float, unit=u.deg, wrap_angle=180*u.deg)

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match=re.escape("Object of type Longitude is not JSON serializable")):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "astropy.coordinates.angles.Longitude"
        assert d["value"] == [3.0, 4.0]
        assert d["unit"] == "deg"
        assert d["wrap_angle"]["value"] == 180.0
        assert d["wrap_angle"]["unit"] == "deg"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, coord.Longitude)
        assert np.array_equal(out, obj)

    def test_latitude(self):
        """Test round-tripping `astropy.coordinates.Latitude`."""
        obj = coord.Latitude([3, 4], dtype=float, unit=u.deg)

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match=re.escape("Object of type Latitude is not JSON serializable")):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "astropy.coordinates.angles.Latitude"
        assert d["value"] == [3.0, 4.0]
        assert d["unit"] == "deg"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, coord.Latitude)
        assert np.array_equal(out, obj)
