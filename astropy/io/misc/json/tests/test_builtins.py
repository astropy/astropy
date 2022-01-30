# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import ast
import json

import pytest

from astropy.io.misc.json import JSONExtendedEncoder, JSONExtendedDecoder

from .test_core import JSONExtendedTestBase


class TestJSONExtendedBuiltins(JSONExtendedTestBase):
    """Tests for serializing builtins with extended JSON encoders and decoders."""

    def test_bytes(self):
        """Test round-tripping `bytes`."""
        obj = b"1234"

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match="Object of type bytes is not JSON serializable"):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "builtins.bytes"
        assert d["value"] == "1234"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, bytes)
        assert out == obj

    def test_complex(self):
        """Test round-tripping `complex`."""
        obj = 1 + 2j

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match="Object of type complex is not JSON serializable"):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "builtins.complex"
        assert d["value"] == [1, 2]

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, complex)
        assert out == obj

    def test_set(self):
        """Test round-tripping `set`."""
        obj = {1, 2, 3}

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match="Object of type set is not JSON serializable"):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "builtins.set"
        assert d["value"] == [1, 2, 3]

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, set)
        assert out == obj

    def test_NotImplemented(self):
        """Test round-tripping `NotImplemented`."""
        obj = NotImplemented

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match="Object of type NotImplementedType is not JSON serializable"):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "builtins.NotImplemented"
        assert d["value"] == "NotImplemented"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, type(NotImplemented))
        assert out == obj

    def test_Ellipsis(self):
        """Test round-tripping `Ellipsis`."""
        obj = Ellipsis

        # Raises errors without extended encoder
        with pytest.raises(TypeError, match="Object of type ellipsis is not JSON serializable"):
            json.dumps(obj)

        # Works with the extended encoder
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)
        d = ast.literal_eval(serialized)
        assert d["__class__"] == "builtins.Ellipsis"
        assert d["value"] == "Ellipsis"

        # Comes back partially processed without extended decoder
        out = json.loads(serialized)
        assert isinstance(out, dict)

        # Roundtrips
        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, type(Ellipsis))
        assert out == obj
