# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import ast
import json
import re

import numpy as np
import pytest

from astropy.io.misc.json import JSONExtendedEncoder, JSONExtendedDecoder

from .test_core import JSONExtendedTestBase


class TestJSONExtendedNumPy(JSONExtendedTestBase):
    """Tests for serializing builtins with extended JSON encoders and decoders."""

    def test_number(self):
        """Test round-tripping `numpy.number`."""
        obj = np.int64(10)

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match="Object of type int64 is not JSON serializable"):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "numpy.int64"
        assert d["value"] == "10"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, np.int64)
        assert out == obj

    def test_dtype_simple(self):
        """Test round-tripping `numpy.dtype`."""
        obj = np.dtype("int64")

        # Raises errors without extended encoder
        # TODO! "Object of type dtype[int64] is not JSON serializable" when py3.9+
        with pytest.raises(TypeError, match=re.escape("Object of type")):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "numpy.dtype"
        assert d["value"] == "int64"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, np.dtype)
        assert out == obj

    @pytest.mark.skip("TODO!")
    def test_dtype_structured(self):
        """Test round-tripping structured `numpy.dtype`."""
        obj = np.dtype([("f1", np.int64), ("f2", np.float16)])

    def test_ndarray_simple(self):
        """Test round-tripping `numpy.ndarray`."""
        obj = np.array([3, 4], dtype=float)

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match=re.escape("Object of type ndarray is not JSON serializable")):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "numpy.ndarray"
        assert d["value"] == [3.0, 4.0]

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, np.ndarray)
        assert np.array_equal(out, obj)

    @pytest.mark.skip("TODO!")
    def test_ndarray_structured(self):
        """Test round-tripping structured `numpy.ndarray`."""
        dt = np.dtype([("f1", np.int64), ("f2", np.float16)])
        obj = np.array((1, 3.0), dtype=dt)
