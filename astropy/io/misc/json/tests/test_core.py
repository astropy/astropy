# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import abc
import json

import numpy as np
import pytest

from astropy.io.misc.json.core import JSONExtendedEncoder, JSONExtendedDecoder


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
        assert out["value"] == self._serialized_value

    def test_roundtrips_with_extended_decoder(self, obj, obj_type):

        serialized = json.dumps(obj, cls=JSONExtendedEncoder)

        out = json.loads(serialized, cls=JSONExtendedDecoder)
        assert isinstance(out, obj_type)
        assert np.all(out == obj)
