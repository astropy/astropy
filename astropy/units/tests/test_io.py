# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import numpy as np
import pytest

import astropy.units as u
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
        self._obj = u.km * u.eV**2
        self._serialized_value = 'eV2 km'

# -------------------------------------------------------------------


class TestJSONExtended_Quantity(JSONExtendedTestBase):

    def setup_class(self):
        self._type = u.Quantity
        self._obj = u.Quantity([3, 4], dtype=float, unit=u.km)
        self._serialized_value = {'!': 'numpy.ndarray', 'value': [3.0, 4.0], 'dtype': 'float64'}


@pytest.mark.skip("TODO!")
class TestJSONExtended_StructuredQuantity(JSONExtendedTestBase):

    def setup_class(self):
        self._type = u.Quantity
        dt = np.unit([("f1", np.int64), ("f2", np.float16)])
        obj = u.Quantity((1, 3.0), dtype=dt, unit=u.Unit("(u.km, u.eV)"))
        self._obj = obj
        # self._serialized_value = "int64"
