# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
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
    from asdf.tests.helpers import yaml_to_asdf

from astropy import units as u

# TODO: Implement defunit


def test_unit():
    yaml = """
unit: !unit/unit-1.0.0 "2.1798721  10-18kg m2 s-2"
    """

    buff = yaml_to_asdf(yaml)
    with asdf.open(buff) as ff:
        assert ff.tree["unit"].is_equivalent(u.Ry)

        buff2 = io.BytesIO()
        ff.write_to(buff2)

    buff2.seek(0)
    with asdf.open(buff2) as ff:
        assert ff.tree["unit"].is_equivalent(u.Ry)
