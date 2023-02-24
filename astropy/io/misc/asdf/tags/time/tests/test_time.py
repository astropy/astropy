# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import pytest

asdf = pytest.importorskip("asdf")

import datetime

import asdf.schema as asdf_schema
import numpy as np
from asdf import AsdfFile, tagged, yamlutil
from asdf.exceptions import AsdfDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        category=AsdfDeprecationWarning,
        message=r"asdf.tests.helpers is deprecated.*",
    )
    from asdf.tests.helpers import assert_roundtrip_tree

from astropy import time


def _flatten_combiners(schema):
    newschema = dict()

    def add_entry(path, schema, combiner):
        # TODO: Simplify?
        cursor = newschema
        for i in range(len(path)):
            part = path[i]
            if isinstance(part, int):
                cursor = cursor.setdefault("items", [])
                while len(cursor) <= part:
                    cursor.append({})
                cursor = cursor[part]
            elif part == "items":
                cursor = cursor.setdefault("items", dict())
            else:
                cursor = cursor.setdefault("properties", dict())
                if i < len(path) - 1 and isinstance(path[i + 1], int):
                    cursor = cursor.setdefault(part, [])
                else:
                    cursor = cursor.setdefault(part, dict())

        cursor.update(schema)


def test_time(tmpdir):
    time_array = time.Time(np.arange(100), format="unix")

    tree = {"large_time_array": time_array}

    assert_roundtrip_tree(tree, tmpdir)


def test_time_with_location(tmpdir):
    # See https://github.com/spacetelescope/asdf/issues/341
    from astropy import units as u
    from astropy.coordinates.earth import EarthLocation

    location = EarthLocation(x=[1, 2] * u.m, y=[3, 4] * u.m, z=[5, 6] * u.m)

    t = time.Time([1, 2], location=location, format="cxcsec")

    tree = {"time": t}

    assert_roundtrip_tree(tree, tmpdir)


def test_time_with_location_1_0_0(tmpdir):
    from astropy import units as u
    from astropy.coordinates.earth import EarthLocation

    location = EarthLocation(x=6378100 * u.m, y=0 * u.m, z=0 * u.m)

    t = time.Time("J2000.000", location=location, format="jyear_str")

    tree = {"time": t}

    # The version refers to ASDF Standard 1.0.0, which includes time-1.0.0
    assert_roundtrip_tree(tree, tmpdir, init_options={"version": "1.0.0"})


def test_isot(tmpdir):
    isot = time.Time("2000-01-01T00:00:00.000")

    tree = {"time": isot}

    assert_roundtrip_tree(tree, tmpdir)

    ff = asdf.AsdfFile(tree)
    tree = yamlutil.custom_tree_to_tagged_tree(ff.tree, ff)
    if isinstance(tree["time"], str):
        assert str(tree["time"]) == isot.value
    elif isinstance(tree["time"], dict):
        assert str(tree["time"]["value"]) == isot.value
        assert str(tree["time"]["base_format"]) == "isot"
    else:
        assert False


def test_isot_array(tmpdir):
    tree = {"time": time.Time(["2001-01-02T12:34:56", "2001-02-03T00:01:02"])}

    assert_roundtrip_tree(tree, tmpdir)


def test_time_tag():
    schema = asdf_schema.load_schema(
        "http://stsci.edu/schemas/asdf/time/time-1.1.0", resolve_references=True
    )
    schema = _flatten_combiners(schema)

    date = time.Time(datetime.datetime.now())
    tree = {"date": date}
    asdf = AsdfFile(tree=tree)
    instance = yamlutil.custom_tree_to_tagged_tree(tree["date"], asdf)

    asdf_schema.validate(instance, schema=schema)

    tag = "tag:stsci.edu:asdf/time/time-1.1.0"
    date = tagged.tag_object(tag, date)
    tree = {"date": date}
    asdf = AsdfFile(tree=tree)
    instance = yamlutil.custom_tree_to_tagged_tree(tree["date"], asdf)

    asdf_schema.validate(instance, schema=schema)
