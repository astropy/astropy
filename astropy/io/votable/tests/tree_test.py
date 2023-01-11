# Licensed under a 3-clause BSD style license - see LICENSE.rst
import io

import pytest

from astropy.io.votable import tree
from astropy.io.votable.exceptions import W07, W08, W21, W41
from astropy.io.votable.table import parse
from astropy.io.votable.tree import Resource, VOTableFile
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyDeprecationWarning


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
    votable = tree.VOTableFile()
    # ...with one resource...
    resource = tree.Resource()
    votable.resources.append(resource)

    # ... with one table
    table = tree.Table(votable)
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
    VOTableFile.__init__ allows versions of '1.0', '1.1', '1.2', '1.3' and '1.4'.
    The '1.0' is curious since other checks in parse() and the version setter do not allow '1.0'.
    This test confirms that behavior for now.  A future change may remove the '1.0'.
    """
    # Exercise the checks in __init__
    with pytest.warns(AstropyDeprecationWarning):
        VOTableFile(version="1.0")
    for version in ("1.1", "1.2", "1.3", "1.4"):
        VOTableFile(version=version)
    for version in ("0.9", "2.0"):
        with pytest.raises(
            ValueError, match=r"should be in \('1.0', '1.1', '1.2', '1.3', '1.4'\)."
        ):
            VOTableFile(version=version)

    # Exercise the checks in the setter
    vot = VOTableFile()
    for version in ("1.1", "1.2", "1.3", "1.4"):
        vot.version = version
    for version in ("1.0", "2.0"):
        with pytest.raises(
            ValueError, match=r"supports VOTable versions '1.1', '1.2', '1.3', '1.4'$"
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

    # Invalid versions
    for bversion in (b"1.0", b"2.0"):
        with pytest.warns(W21):
            parse(
                io.BytesIO(begin + bversion + middle + bversion + end),
                verify="exception",
            )


def votable_xml_string(version):
    votable_file = VOTableFile(version=version)
    votable_file.resources.append(Resource())

    xml_bytes = io.BytesIO()
    votable_file.to_xml(xml_bytes)
    xml_bytes.seek(0)
    bstring = xml_bytes.read()
    s = bstring.decode("utf-8")
    return s


def test_votable_tag():
    xml = votable_xml_string("1.1")
    assert 'xmlns="http://www.ivoa.net/xml/VOTable/v1.1"' in xml
    assert 'xsi:noNamespaceSchemaLocation="http://www.ivoa.net/xml/VOTable/v1.1"' in xml

    xml = votable_xml_string("1.2")
    assert 'xmlns="http://www.ivoa.net/xml/VOTable/v1.2"' in xml
    assert 'xsi:noNamespaceSchemaLocation="http://www.ivoa.net/xml/VOTable/v1.2"' in xml

    xml = votable_xml_string("1.3")
    assert 'xmlns="http://www.ivoa.net/xml/VOTable/v1.3"' in xml
    assert 'xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 '
    assert 'http://www.ivoa.net/xml/VOTable/VOTable-1.3.xsd"' in xml

    xml = votable_xml_string("1.4")
    assert 'xmlns="http://www.ivoa.net/xml/VOTable/v1.3"' in xml
    assert 'xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 '
    assert 'http://www.ivoa.net/xml/VOTable/VOTable-1.4.xsd"' in xml
