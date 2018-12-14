# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest

from asdf.tests.helpers import assert_roundtrip_tree

from astropy.time import Time


@pytest.mark.parametrize('fmt', Time.FORMATS.keys())
def test_timedelta(fmt, tmpdir):

    t1 = Time(Time.now(), format=fmt)
    t2 = Time(Time.now(), format=fmt)

    td = t2 - t1
    tree = dict(timedelta=td)
    assert_roundtrip_tree(tree, tmpdir)
