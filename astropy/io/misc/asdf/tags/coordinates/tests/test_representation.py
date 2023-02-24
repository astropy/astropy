# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

import pytest

asdf = pytest.importorskip("asdf")

from asdf.exceptions import AsdfDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        category=AsdfDeprecationWarning,
        message=r"asdf.tests.helpers is deprecated.*",
    )
    from asdf.tests.helpers import assert_roundtrip_tree

from numpy.random import randint, random

import astropy.coordinates.representation as r
import astropy.units as u
from astropy.coordinates import Angle


@pytest.fixture(params=filter(lambda x: "Base" not in x, r.__all__))
def representation(request):
    rep = getattr(r, request.param)

    angle_unit = u.deg
    other_unit = u.km

    kwargs = {}
    arr_len = randint(1, 100)
    for aname, atype in rep.attr_classes.items():
        if issubclass(atype, Angle):
            value = ([random()] * arr_len) * angle_unit
        else:
            value = ([random()] * arr_len) * other_unit

        kwargs[aname] = value

    return rep(**kwargs)


def test_representations(tmpdir, representation):
    tree = {"representation": representation}
    assert_roundtrip_tree(tree, tmpdir)
