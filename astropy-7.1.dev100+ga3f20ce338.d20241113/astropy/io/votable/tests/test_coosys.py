# Licensed under a 3-clause BSD style license - see LICENSE.rst
import io

import pytest

from astropy.io.votable.exceptions import E16, W27, W57
from astropy.io.votable.table import parse
from astropy.utils.data import get_pkg_data_filename

COOSYS_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v{}" version="{}">
  <RESOURCE>
    {}
  </RESOURCE>
</VOTABLE>
"""
COOSYS_REFPOSITION = (
    '<COOSYS ID="my_coosys" epoch="J2015.5" system="ICRS" refposition="BARYCENTER"/>'
)
COOSYS_NO_REFPOSITION = '<COOSYS ID="my_coosys" epoch="J2015.5" system="ICRS"/>'

COOSYS_SYSTEM_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3" version="{}">
  <RESOURCE>
    <COOSYS ID="my_coosys" epoch="J2015.5" system="{}"/>
  </RESOURCE>
</VOTABLE>
"""


def test_accept_refposition():
    for ns_version, version, coosys, expected_warnings in [
        ("1.1", "1.1", COOSYS_NO_REFPOSITION, None),
        ("1.2", "1.2", COOSYS_NO_REFPOSITION, W27),
        ("1.3", "1.3", COOSYS_NO_REFPOSITION, None),
        ("1.3", "1.4", COOSYS_NO_REFPOSITION, None),
        ("1.3", "1.5", COOSYS_NO_REFPOSITION, None),
        ("1.1", "1.1", COOSYS_REFPOSITION, W57),
        ("1.2", "1.2", COOSYS_REFPOSITION, (W27, W57)),
        ("1.3", "1.3", COOSYS_REFPOSITION, W57),
        ("1.3", "1.4", COOSYS_REFPOSITION, W57),
        ("1.3", "1.5", COOSYS_REFPOSITION, None),
    ]:
        xml_string = COOSYS_TEMPLATE.format(ns_version, version, coosys)
        input = io.BytesIO(bytes(xml_string, "utf-8"))
        if expected_warnings is not None:
            with pytest.warns(expected_warnings) as record:
                vot = parse(input, verify="warn")
            # If expected_warnings is a tuple, pytest.warns checked for at least one of the warnings,
            # so make sure we got the right number of them.
            if isinstance(expected_warnings, tuple):
                assert len(record) == len(expected_warnings)
        else:
            vot = parse(input, verify="warn")


def test_coosys_system():
    """
    Check that prior to v1.5, an error is issues for a system not in the hard-coded list.
    With 1.5 the vocabulary is external and we're not checking that yet.
    """
    for version, system, expected_warning in [
        ("1.4", "ICRS", None),
        ("1.4", "supergalactic", None),
        ("1.4", "InvalidSystem", E16),
        ("1.5", "ICRS", None),
        ("1.5", "supergalactic", None),
        ("1.5", "InvalidSystem", None),
    ]:
        xml_string = COOSYS_SYSTEM_TEMPLATE.format(version, system)
        input = io.BytesIO(bytes(xml_string, "utf-8"))
        if expected_warning is not None:
            with pytest.warns(expected_warning):
                vot = parse(input, verify="warn")
        else:
            vot = parse(input, verify="warn")


def _coosys_tests(votable):
    assert len(list(votable.iter_coosys())) == 2

    coosys = votable.get_coosys_by_id("coosys1")
    assert coosys.system == "ICRS"
    assert coosys.equinox == "J2000"
    assert coosys.epoch == "J2015.5"
    assert coosys.refposition == "BARYCENTER"

    coosys = votable.get_coosys_by_id("coosys2")
    assert coosys.system == "FK4"
    assert coosys.equinox == "B1950"
    assert coosys.epoch == "B1950.0"
    assert coosys.refposition == "UNKNOWN"


def test_coosys():
    votable = parse(get_pkg_data_filename("data/coosys.xml"))
    _coosys_tests(votable)


def test_coosys_roundtrip():
    orig_votable = parse(get_pkg_data_filename("data/coosys.xml"))
    bio = io.BytesIO()
    orig_votable.to_xml(bio)
    bio.seek(0)
    votable = parse(bio)
    _coosys_tests(votable)
