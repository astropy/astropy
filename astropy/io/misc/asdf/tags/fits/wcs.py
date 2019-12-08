# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from astropy.wcs import WCS
from astropy.io.misc.asdf.types import AstropyType


class WCSType(AstropyType):
    name = "fits/wcs"
    types = ["astropy.wcs.WCS"]
    requires = ["astropy"]

    @classmethod
    def from_tree(cls, data, ctx):
        return WCS(header=data)

    @classmethod
    def to_tree(cls, wcsobj, ctx):
        return dict(wcsobj.to_header())

    @classmethod
    def assert_equal(cls, old, new):
        for wcsa, wcsb in zip(old, new):
            assert dict(wcsa.to_header()) == dict(wcsb.to_header())

