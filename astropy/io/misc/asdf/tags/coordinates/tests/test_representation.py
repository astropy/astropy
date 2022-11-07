# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

asdf = pytest.importorskip("asdf")

from asdf.tests.helpers import assert_roundtrip_tree  # noqa: E402
from numpy.random import randint, random  # noqa: E402

import astropy.coordinates.representation as r  # noqa: E402
import astropy.units as u  # noqa: E402
from astropy.coordinates import Angle  # noqa: E402


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
