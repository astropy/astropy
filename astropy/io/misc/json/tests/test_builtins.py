# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import ast
import json

from .test_core import JSONExtendedTestBase


class TestJSONExtended_Bytes(JSONExtendedTestBase):

    def setup_class(self):
        self._type = bytes
        self._obj = b"1234"
        self._serialized_value = "1234"


class TestJSONExtended_Complex(JSONExtendedTestBase):

    def setup_class(self):
        self._type = complex
        self._obj = 1 + 2j
        self._serialized_value = [1, 2]


class TestJSONExtended_Set(JSONExtendedTestBase):

    def setup_class(self):
        self._type = set
        self._obj = {1, 2, 3}
        self._serialized_value = [1, 2, 3]


class TestJSONExtended_NotImplemented(JSONExtendedTestBase):

    def setup_class(self):
        self._type = type(NotImplemented)
        self._obj = NotImplemented
        self._serialized_class = "builtins.NotImplemented"
        self._serialized_value = "NotImplemented"


class TestJSONExtended_Ellipsis(JSONExtendedTestBase):

    def setup_class(self):
        self._type = type(Ellipsis)
        self._obj = Ellipsis
        self._serialized_class = "builtins.Ellipsis"
        self._serialized_value = "Ellipsis"
