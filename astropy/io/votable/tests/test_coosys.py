# Licensed under a 3-clause BSD style license - see LICENSE.rst
import io
from contextlib import nullcontext

import pytest

from astropy.io.votable import tree
from astropy.io.votable.exceptions import W07, W08, W21, W27, W41, W57
from astropy.io.votable.table import parse
from astropy.io.votable.tree import Resource, VOTableFile
from astropy.tests.helper import PYTEST_LT_8_0
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
NO_COOSYS = ""


def test_accept_refposition():
    for ns_version, version, coosys, err in [
        ("1.1", "1.1", COOSYS_NO_REFPOSITION, None),
        ("1.2", "1.2", COOSYS_NO_REFPOSITION, (W27)),
        ("1.3", "1.3", COOSYS_NO_REFPOSITION, None),
        ("1.3", "1.4", COOSYS_NO_REFPOSITION, None),
        ("1.3", "1.5", COOSYS_NO_REFPOSITION, None),
        ("1.1", "1.1", COOSYS_REFPOSITION, (W57)),
        ("1.2", "1.2", COOSYS_REFPOSITION, (W27, W57)),
        ("1.3", "1.3", COOSYS_REFPOSITION, (W57)),
        ("1.3", "1.4", COOSYS_REFPOSITION, (W57)),
        ("1.3", "1.5", COOSYS_REFPOSITION, None),
    ]:
        xml_string = COOSYS_TEMPLATE.format(ns_version, version, coosys)
        input = io.BytesIO(bytes(xml_string, "utf-8"))
        if err is not None:
            with pytest.warns() as record:
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


if __name__ == "__main__":
    test_accept_refposition()


def test_check_astroyear_fail():
    config = {"verify": "exception"}
    field = tree.Field(None, name="astroyear", arraysize="1")
    with pytest.raises(W07):
        tree.check_astroyear("X2100", field, config)


def test_string_fail():
    config = {"verify": "exception"}
    with pytest.raises(W08):
        tree.check_string(42, "foo", config)


def test_make_Fields():
    votable = VOTableFile()
    # ...with one resource...
    resource = Resource()
    votable.resources.append(resource)

    # ... with one table
    table = tree.TableElement(votable)
    resource.tables.append(table)

    table.fields.extend(
        [tree.Field(votable, name="Test", datatype="float", unit="mag")]
    )


def test_unit_format():
    data = parse(get_pkg_data_filename("data/irsa-nph-error.xml"))
    assert data._config["version"] == "1.0"
    assert tree._get_default_unit_format(data._config) == "cds"
    data = parse(get_pkg_data_filename("data/names.xml"))
    assert data._config["version"] == "1.1"
    assert tree._get_default_unit_format(data._config) == "cds"
    data = parse(get_pkg_data_filename("data/gemini.xml"))
    assert data._config["version"] == "1.2"
    assert tree._get_default_unit_format(data._config) == "cds"
    data = parse(get_pkg_data_filename("data/binary2_masked_strings.xml"))
    assert data._config["version"] == "1.3"
    assert tree._get_default_unit_format(data._config) == "cds"
    data = parse(get_pkg_data_filename("data/timesys.xml"))
    assert data._config["version"] == "1.4"
    assert tree._get_default_unit_format(data._config) == "vounit"


def test_namespace_warning():
    """
    A version 1.4 VOTable must use the same namespace as 1.3.
    (see https://www.ivoa.net/documents/VOTable/20191021/REC-VOTable-1.4-20191021.html#ToC16).
    """
    bad_namespace = b"""<?xml version="1.0" encoding="utf-8"?>
        <VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.4"
                               xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
          <RESOURCE/>
        </VOTABLE>
    """
    with pytest.warns(W41):
        parse(io.BytesIO(bad_namespace), verify="exception")

    good_namespace_14 = b"""<?xml version="1.0" encoding="utf-8"?>
        <VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"
                               xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
          <RESOURCE/>
        </VOTABLE>
    """
    parse(io.BytesIO(good_namespace_14), verify="exception")

    good_namespace_13 = b"""<?xml version="1.0" encoding="utf-8"?>
        <VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"
                               xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
          <RESOURCE/>
        </VOTABLE>
    """
    parse(io.BytesIO(good_namespace_13), verify="exception")


def test_version():
    """
    VOTableFile.__init__ allows versions of '1.1' through '1.5'.
    VOTableFile.__init__ does not allow version of '1.0' anymore and now raises a ValueError as it does to other versions not supported.
    """
    # Exercise the checks in __init__
    for version in ("1.1", "1.2", "1.3", "1.4", "1.5"):
        VOTableFile(version=version)
    for version in ("0.9", "1.0", "1.6", "2.0"):
        with pytest.raises(
            ValueError, match=r"should be in \('1.1', '1.2', '1.3', '1.4', '1.5'\)."
        ):
            VOTableFile(version=version)

    # Exercise the checks in the setter
    vot = VOTableFile()
    for version in ("1.1", "1.2", "1.3", "1.4"):
        vot.version = version
    for version in ("1.0", "1.6", "2.0"):
        with pytest.raises(
            ValueError,
            match=r"supports VOTable versions '1.1', '1.2', '1.3', '1.4', '1.5'$",
        ):
            vot.version = version

    # Exercise the checks in the parser.
    begin = b'<?xml version="1.0" encoding="utf-8"?><VOTABLE version="'
    middle = b'" xmlns="http://www.ivoa.net/xml/VOTable/v'
    end = (
        b'" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"><RESOURCE/></VOTABLE>'
    )

    # Valid versions
    for bversion in (b"1.1", b"1.2", b"1.3"):
        parse(
            io.BytesIO(begin + bversion + middle + bversion + end), verify="exception"
        )
    parse(io.BytesIO(begin + b"1.4" + middle + b"1.3" + end), verify="exception")

    if PYTEST_LT_8_0:
        ctx = nullcontext()
    else:
        ctx = pytest.warns(W41)

    # Invalid versions
    for bversion in (b"1.0", b"2.0"):
        with pytest.warns(W21), ctx:
            parse(
                io.BytesIO(begin + bversion + middle + bversion + end),
                verify="exception",
            )
