# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


from asdf.yamlutil import custom_tree_to_tagged_tree

from astropy.time import TimeDelta
from astropy.time.tests.test_delta import (allclose_jd, allclose_jd2,
                                           allclose_sec)

from ...types import AstropyType


__all__ = ['TimeDeltaType']


class TimeDeltaType(AstropyType):
    name = 'time/timedelta'
    types = [TimeDelta]
    version = '1.0.0'


    @classmethod
    def to_tree(cls, obj, ctx):
        return custom_tree_to_tagged_tree(obj.info._represent_as_dict(), ctx)


    @classmethod
    def from_tree(cls, node, ctx):
        return TimeDelta.info._construct_from_dict(node)


    @classmethod
    def assert_equal(cls, old, new):
        assert allclose_jd(old.jd, new.jd)
        assert allclose_jd2(old.jd2, new.jd2)
        assert allclose_sec(old.sec, new.sec)
