# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

import pytest

asdf = pytest.importorskip("asdf")

from asdf.tests import helpers  # noqa: E402

from astropy import units as u  # noqa: E402

# TODO: Implement defunit


def test_unit():
    yaml = """
unit: !unit/unit-1.0.0 "2.1798721  10-18kg m2 s-2"
    """

    buff = helpers.yaml_to_asdf(yaml)
    with asdf.open(buff) as ff:
        assert ff.tree["unit"].is_equivalent(u.Ry)

        buff2 = io.BytesIO()
        ff.write_to(buff2)

    buff2.seek(0)
    with asdf.open(buff2) as ff:
        assert ff.tree["unit"].is_equivalent(u.Ry)
