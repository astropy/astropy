# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import pytest

asdf = pytest.importorskip('asdf')

import numpy as np

import astropy.units as u

from astropy.coordinates import Angle
import astropy.coordinates.representation as r

from asdf.tests.helpers import assert_roundtrip_tree


@pytest.fixture(params=filter(lambda x: "Base" not in x, r.__all__))
def representation(request):
    rep = getattr(r, request.param)

    angle_unit = u.deg
    other_unit = u.km

    kwargs = {}
    arr_len = 16  # arbitrary number
    rng = np.random.default_rng(seed=42)
    for aname, atype in rep.attr_classes.items():
        if issubclass(atype, Angle):
            unit = angle_unit
        else:
            unit = other_unit

        kwargs[aname] = rng.random(size=arr_len) * unit

    return rep(**kwargs)


def test_representations(tmpdir, representation):
    tree = {'representation': representation}
    assert_roundtrip_tree(tree, tmpdir)
