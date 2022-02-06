# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

import astropy.units as u
import astropy.coordinates as coord
from astropy.io.misc.json.tests.test_core import JSONExtendedTestBase


class TestJSONExtended_Longitude(JSONExtendedTestBase):

    def setup_class(self):
        self._type = coord.Longitude
        self._obj = coord.Longitude([3, 4], dtype=float, unit=u.deg, wrap_angle=180*u.deg)
        self._serialized_value = {'!': 'numpy.ndarray', 'value': [3.0, 4.0], 'dtype': 'float64'}

    # def test_
    #     assert d["unit"] == "deg"
    #     assert d["wrap_angle"]["value"] == 180.0
    #     assert d["wrap_angle"]["unit"] == "deg"


class TestJSONExtended_Latitude(JSONExtendedTestBase):

    def setup_class(self):
        self._type = coord.Latitude
        self._obj = coord.Latitude([3, 4], dtype=float, unit=u.deg)
        self._serialized_value = {'!': 'numpy.ndarray', 'value': [3.0, 4.0], 'dtype': 'float64'}
