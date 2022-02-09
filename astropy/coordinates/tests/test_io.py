# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import json

import pytest

import astropy.units as u
import astropy.coordinates as coord
from astropy.io.misc.json import JSONExtendedEncoder
from astropy.units.tests.test_io import TestJSONExtended_Quantity


class TestJSONExtended_Longitude(TestJSONExtended_Quantity):

    def setup_class(self):
        self._type = coord.Longitude
        self._obj = coord.Longitude([3, 4], dtype=float, unit=u.deg, wrap_angle=180*u.deg)
        self._serialized_value = {'!': 'numpy.ndarray', 'value': ["3.0", "4.0"], 'dtype': 'float64'}

    def test_wrap_angle(self, obj):
        serialized = json.dumps(obj, cls=JSONExtendedEncoder)
        out = json.loads(serialized)

        assert out["wrap_angle"]["value"] == 180.0
        assert out["wrap_angle"]["unit"] == "deg"


class TestJSONExtended_Latitude(TestJSONExtended_Quantity):

    def setup_class(self):
        self._type = coord.Latitude
        self._obj = coord.Latitude([3, 4], dtype=float, unit=u.deg)
        self._serialized_value = {'!': 'numpy.ndarray', 'value': ["3.0", "4.0"], 'dtype': 'float64'}
