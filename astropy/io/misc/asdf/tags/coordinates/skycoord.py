# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import numpy as np

from asdf.yamlutil import custom_tree_to_tagged_tree, tagged_tree_to_custom_tree

from astropy.coordinates import SkyCoord
from astropy.table.tests.test_operations import skycoord_equal

from ...types import AstropyType


class SkyCoordType(AstropyType):
    name = 'coordinates/skycoord'
    types = [SkyCoord]
    version = "1.0.0"

    @classmethod
    def to_tree(cls, obj, ctx):
        return custom_tree_to_tagged_tree(obj.info._represent_as_dict(), ctx)

    @classmethod
    def from_tree(cls, tree, ctx):
        return SkyCoord.info._construct_from_dict(tree)

    @classmethod
    def assert_equal(cls, old, new):
        assert skycoord_equal(old, new)
