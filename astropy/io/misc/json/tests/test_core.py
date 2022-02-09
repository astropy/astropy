# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import abc
import json
from collections.abc import Mapping

import numpy as np
import pytest

from astropy.io.misc.json.core import JSONExtendedEncoder, JSONExtendedDecoder


def _array_close(obj1, obj2):
    try:
        eq = np.allclose(obj1, obj2)
    except (TypeError, DeprecationWarning, np.VisibleDeprecationWarning):
        try:
            eq = (obj1 == obj2)
        except ValueError:  # Some element-wise failures. Maybe mappings?
            eq = False
    return eq


def _recursive_eq(obj1, obj2):

    if not isinstance(obj1, Mapping):
        return _array_close(obj1, obj2)

    elif not isinstance(obj2, Mapping):
        return False

    elif set(obj1.keys()) != set(obj2.keys()):
        return False

    for k, v in obj1.items():
        eq = _array_close(obj1, obj2)
        if not eq and isinstance(v, Mapping) and isinstance(obj2[k], Mapping):
            eq = _recursive_eq(v, obj2[k])
    if not np.all(eq):
        return False

    return True


class JSONExtendedTestBase(metaclass=abc.ABCMeta):
    """Base for testing JSON extended encoders and decoders"""

    @abc.abstractmethod
    def setup_class(self):
        """Set up test.

        Should define ``_type``, ``_obj``, ``_serialized_value``
        """
        # self._type
        # self._obj
        # self._serialized_value

    @pytest.fixture(scope="class")
    def obj_type(self):
        """Pytest fixture returning the type of the object to be tested."""
        return self._type

    @pytest.fixture(scope="class")
    def obj(self):
        """Pytest fixture returning the object to be tested."""
        return self._obj

    # ===============================================================

    def test_errors_without_extended_encoder(self, obj, obj_type):
        """Test objects of the tested type cannot normally be JSON serialized."""
        with pytest.raises(TypeError, match="is not JSON serializable"):
            json.dumps(obj)

    def test_works_with_extended_encoder(self, obj, obj_type):
        # Serialize
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        assert isinstance(serialized, str)

        # Partially unserialize (without extended decoder)
        out = json.loads(serialized)
        scls = getattr(self, "_serialized_class", f"{obj_type.__module__}.{obj_type.__qualname__}")
        assert out["!"] == scls
        assert _recursive_eq(out["value"], self._serialized_value)

    def test_roundtrips_with_extended_decoder(self, obj, obj_type):

        serialized = json.dumps(obj, cls=JSONExtendedEncoder)

        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, obj_type)
        assert np.all(out == obj)
