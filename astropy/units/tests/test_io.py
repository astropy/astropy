# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import json

import numpy as np
import pytest

import astropy.units as u
from astropy.io.misc.json import JSONExtendedEncoder
from astropy.io.misc.json.tests.test_core import JSONExtendedTestBase


class TestJSONExtended_Unit(JSONExtendedTestBase):
    def setup_class(self):
        self._type = u.PrefixUnit
        self._obj = u.km
        self._serialized_value = "km"


class TestJSONExtended_StructuredUnit(JSONExtendedTestBase):
    def setup_class(self):
        self._type = u.StructuredUnit
        self._obj = u.StructuredUnit("(km, km, (eV^2, eV))")
        self._serialized_value = "(km, km, (eV2, eV))"


class TestJSONExtended_CompositeUnit(JSONExtendedTestBase):
    def setup_class(self):
        self._type = u.CompositeUnit
        self._obj = u.km * u.eV ** 2
        self._serialized_value = "eV2 km"


# -------------------------------------------------------------------


class TestJSONExtended_Quantity(JSONExtendedTestBase):
    def setup_class(self):
        self._type = u.Quantity
        self._obj = u.Quantity([3, 4], dtype=float, unit=u.km)
        self._serialized_value = {"!": "numpy.ndarray", "value": ["3.0", "4.0"], "dtype": "float64"}

    def test_unit(self, obj):
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        out = json.loads(serialized)["unit"]
        out = out if not isinstance(out, dict) else out["value"]
        assert out == str(self._obj.unit)


class TestJSONExtended_StructuredVoidQuantity(TestJSONExtended_Quantity):
    def setup_class(self):
        self._type = u.Quantity
        dt = np.dtype([("nu1", float), ("nu2", np.float32)])
        unit = u.Unit("(eV, eV)")
        obj = u.Quantity((0, 0.6), dtype=dt, unit=unit)
        self._obj = obj
        self._serialized_value = {
            "!": "numpy.void",
            "value": ["0.0", "0.6"],
            "dtype": {"value": {"nu1": [{"!": "numpy.dtype", "value": "float64"}, 0],
                                "nu2": [{"!": "numpy.dtype", "value": "float32"}, 8]},
                      "align": False}}


class TestJSONExtended_StructuredArrayQuantity(TestJSONExtended_Quantity):
    def setup_class(self):
        self._type = u.Quantity
        dt = np.dtype([("nu1", float), ("nu2", np.float32)])
        unit = u.Unit("(eV, eV)")
        obj = u.Quantity([(0, 0.6)], dtype=dt, unit=unit)
        self._obj = obj
        self._serialized_value = {
            "!": "numpy.ndarray",
            "value": {"nu1": {"!": "numpy.ndarray", "value": ["0.0"], "dtype": "float64"},
                      "nu2": {"!": "numpy.ndarray", "value": ["0.6"], "dtype": "float32"}
            },
            "dtype": {
                "value": {"nu1": [{"!": "numpy.dtype", "value": "float64"}, 0],
                          "nu2": [{"!": "numpy.dtype", "value": "float32"}, 8]},
                "align": False}
        }
