# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import json

import numpy as np
import pytest

from .test_core import JSONExtendedTestBase
from astropy.io.misc.json.core import JSONExtendedEncoder


class TestJSONExtended_DType(JSONExtendedTestBase):
    def setup_class(self):
        self._type = np.dtype
        self._obj = np.dtype("int64")
        self._serialized_value = "int64"


class TestJSONExtended_DType_Metadata(TestJSONExtended_DType):
    def setup_class(self):
        self._type = np.dtype
        self._obj = np.dtype("int64", metadata=dict(a=1))
        self._serialized_value = "int64"

    def test_metadata(self, obj, obj_type):
        # Serialize
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        # Partially unserialize (without extended decoder)
        out = json.loads(serialized)
        out["metadata"] == dict(a=1)


class TestJSONExtended_DType_Shaped(TestJSONExtended_DType):
    def setup_class(self):
        self._type = np.dtype
        self._obj = np.dtype("10float64", align=True)
        self._serialized_value = ["float64", [10]]


class TestJSONExtended_StructuredDType(TestJSONExtended_DType):
    def setup_class(self):
        self._type = np.dtype
        self._obj = np.dtype([("f1", float, (1, 2)), ("f2", [("n1", float), ("n2", int)])])
        self._serialized_value = {
            "f1": [{"!": "numpy.dtype", "value": ["float64", [1, 2]]}, 0],
            "f2": [
                {
                    "!": "numpy.dtype",
                    "value": {
                        "n1": [{"!": "numpy.dtype", "value": "float64"}, 0],
                        "n2": [{"!": "numpy.dtype", "value": "int64"}, 8],
                    },
                    "align": False,
                },
                16,
            ],
        }

    def test_align(self, obj, obj_type):
        # Serialize
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        # Partially unserialize (without extended decoder)
        out = json.loads(serialized)
        out["align"] == False


# -------------------------------------------------------------------
# Scalars


class TestJSONExtended_Bool(JSONExtendedTestBase):
    def setup_class(self):
        self._type = np.bool_
        self._obj = np.bool_(True)
        self._serialized_value = True


class TestJSONExtended_Number(JSONExtendedTestBase):
    def setup_class(self):
        self._type = np.int64
        self._obj = np.int64(10)
        self._serialized_value = "10"


@pytest.mark.skip("TODO")
class TestJSONExtended_Void_Bytes(JSONExtendedTestBase):
    def setup_class(self):
        self._type = np.void
        self._obj = np.void(b"abcd")
        self._serialized_value = "abcd"


class TestJSONExtended_Void_Structure(JSONExtendedTestBase):
    def setup_class(self):
        self._type = np.void
        self._obj = val = np.array(
            (0, 0.60), dtype=np.dtype([("nu1", float), ("nu2", np.float32)])
        )[()]
        self._serialized_value = [0.0, 0.6]

    def test_dtype(self, obj):
        # Serialize
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        # Partially unserialize (without extended decoder)
        out = json.loads(serialized)
        out["dtype"] == {
            "!": "numpy.dtype",
            "value": {
                "nu1": [{"!": "numpy.dtype", "value": "float64"}, 0],
                "nu2": [{"!": "numpy.dtype", "value": "float32"}, 8],
            },
            "align": False,
        }


# -------------------------------------------------------------------
# Arrays


class TestJSONExtended_NDArray(JSONExtendedTestBase):
    def setup_class(self):
        self._type = np.ndarray
        self._obj = np.array([3, 4], dtype=float)
        self._serialized_value = [3.0, 4.0]


class TestJSONExtended_StructuredNDArray(JSONExtendedTestBase):
    def setup_class(self):
        self._type = np.ndarray

        self._obj = np.array((0, 0.6), dtype=np.dtype([("nu1", float), ("nu2", np.float32)]))
        self._serialized_value = {
            "nu1": {"!": "numpy.ndarray", "value": 0.0, "dtype": "float64"},
            "nu2": {"!": "numpy.ndarray", "value": 0.6000000238418579, "dtype": "float32"},
        }
