# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import numpy as np

from asdf.yamlutil import custom_tree_to_tagged_tree, tagged_tree_to_custom_tree

from astropy.coordinates import SkyCoord

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
        old_dict = old.info._represent_as_dict()
        new_dict = new.info._represent_as_dict()

        assert old_dict.keys() == new_dict.keys()

        for k in old_dict.keys():
            if isinstance(old_dict[k], np.ndarray):
                assert (old_dict[k] == new_dict[k]).all()
            else:
                assert old_dict[k] == new_dict[k]
