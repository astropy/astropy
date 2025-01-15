# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np

from astropy.io.misc.asdf.types import AstropyType
from astropy.time import TimeDelta

__all__ = ["TimeDeltaType"]

allclose_jd = functools.partial(np.allclose, rtol=2.0**-52, atol=0)
allclose_jd2 = functools.partial(
    np.allclose, rtol=2.0**-52, atol=2.0**-52
)  # 20 ps atol
allclose_sec = functools.partial(
    np.allclose, rtol=2.0**-52, atol=2.0**-52 * 24 * 3600
)  # 20 ps atol


class TimeDeltaType(AstropyType):
    name = "time/timedelta"
    types = [TimeDelta]
    version = "1.0.0"

    @classmethod
    def to_tree(cls, obj, ctx):
        return obj.info._represent_as_dict()

    @classmethod
    def from_tree(cls, node, ctx):
        return TimeDelta.info._construct_from_dict(node)

    @classmethod
    def assert_equal(cls, old, new):
        assert allclose_jd(old.jd, new.jd)
        assert allclose_jd2(old.jd2, new.jd2)
        assert allclose_sec(old.sec, new.sec)
