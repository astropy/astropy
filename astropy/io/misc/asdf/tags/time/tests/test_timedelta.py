# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree

from astropy.time import Time, TimeDelta


@pytest.mark.parametrize('fmt', Time.FORMATS.keys())
def test_timedelta(fmt, tmpdir):

    t1 = Time(Time.now(), format=fmt)
    t2 = Time(Time.now(), format=fmt)

    td = t2 - t1
    tree = dict(timedelta=td)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.parametrize('scale', list(TimeDelta.SCALES) + [None])
def test_timedetal_scales(scale, tmpdir):

    tree = dict(timedelta=TimeDelta(0.125, scale=scale))
    assert_roundtrip_tree(tree, tmpdir)
