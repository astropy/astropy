# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

asdf = pytest.importorskip('asdf')
from asdf.tests.helpers import assert_roundtrip_tree

from astropy import units as u
from astropy.time import Time, TimeDelta


@pytest.mark.parametrize('fmt', TimeDelta.FORMATS.keys())
def test_timedelta(fmt, tmpdir):

    t1 = Time(Time.now())
    t2 = Time(Time.now())

    td = TimeDelta(t2 - t1, format=fmt)
    tree = dict(timedelta=td)
    assert_roundtrip_tree(tree, tmpdir)


@pytest.mark.parametrize('scale', list(TimeDelta.SCALES) + [None])
def test_timedelta_scales(scale, tmpdir):

    tree = dict(timedelta=TimeDelta(0.125, scale=scale, format="jd"))
    assert_roundtrip_tree(tree, tmpdir)


def test_timedelta_vector(tmpdir):

    tree = dict(timedelta=TimeDelta([1, 2] * u.day))
    assert_roundtrip_tree(tree, tmpdir)
