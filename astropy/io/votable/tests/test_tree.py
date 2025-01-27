# Licensed under a 3-clause BSD style license - see LICENSE.rst
import filecmp
import io

import numpy as np
import pytest

from astropy.io.votable import tree
from astropy.io.votable.exceptions import E26, W07, W08, W21, W41
from astropy.io.votable.table import parse
from astropy.io.votable.tree import MivotBlock, Resource, VOTableFile
from astropy.utils.data import get_pkg_data_filename


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


def test_votable_values_empty_min_max():
    """Regression test for https://github.com/astropy/astropy/issues/16825"""
    with_empty_minmax = b"""<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3" version="1.4">
        <RESOURCE type="results">
          <TABLE name="main">
            <PARAM name="break" datatype="int" value=""/>
          <FIELD ID="hd" datatype="int" name="hd" ucd="meta.id;meta.main">
            <DESCRIPTION>HD number for this object</DESCRIPTION>
            <VALUES null="-2147483648">
              <MIN value=""/>
              <MAX value=""/>
            </VALUES>
          </FIELD>
          <DATA>
            <BINARY>
              <STREAM encoding="base64">AAMNIg==</STREAM>
            </BINARY>
          </DATA>
        </TABLE>
      </RESOURCE>
    </VOTABLE>
    """
    parse(io.BytesIO(with_empty_minmax), verify="exception")


def _assert_minmax_match(val, expected_val):
    """val can be number or numpy array; expected_val can be number or list"""
    if isinstance(val, np.ndarray):
        assert np.allclose(val, np.array(expected_val))
    else:
        assert val == expected_val


@pytest.mark.parametrize(
    ("testvals"),
    [
        {
            "dt": "float",
            "min": "1.2",
            "max": "3.4",
            "expected_min": 1.2,
            "expected_max": 3.4,
        },
        {
            "dt": "float",
            "min": "1.2 3.4",
            "max": "3.4 5.6",
            "expected_min": [1.2, 3.4],
            "expected_max": [3.4, 5.6],
        },
        {
            "dt": "double",
            "min": "1.2",
            "max": "3.4",
            "expected_min": 1.2,
            "expected_max": 3.4,
        },
        {
            "dt": "double",
            "min": "1.2 3.4",
            "max": "3.4 5.6",
            "expected_min": [1.2, 3.4],
            "expected_max": [3.4, 5.6],
        },
        {
            "dt": "unsignedByte",
            "min": "1",
            "max": "3",
            "expected_min": 1,
            "expected_max": 3,
        },
        {
            "dt": "unsignedByte",
            "min": "1 3",
            "max": "3 5",
            "expected_min": [1, 3],
            "expected_max": [3, 5],
        },
        {"dt": "short", "min": "1", "max": "3", "expected_min": 1, "expected_max": 3},
        {
            "dt": "short",
            "min": "1 3",
            "max": "3 5",
            "expected_min": [1, 3],
            "expected_max": [3, 5],
        },
        {"dt": "int", "min": "1", "max": "3", "expected_min": 1, "expected_max": 3},
        {
            "dt": "int",
            "min": "1 3",
            "max": "3 5",
            "expected_min": [1, 3],
            "expected_max": [3, 5],
        },
        {"dt": "long", "min": "1", "max": "3", "expected_min": 1, "expected_max": 3},
        {
            "dt": "long",
            "min": "1 3",
            "max": "3 5",
            "expected_min": [1, 3],
            "expected_max": [3, 5],
        },
    ],
)
def test_min_max_with_arrays(testvals):
    """When a FIELD/PARAM is an array type, ensure we accept (without warnings/exceptions)
    MIN and MAX values both as scalars (preferred except for some xtypes) or arrays."""
    minmax_template = """<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3" version="1.5">
        <RESOURCE type="results">
            <TABLE ID="result" name="result">
                <{param_or_field} arraysize="2" datatype="{dt}" name="param_or_field" xtype="interval">
                    <VALUES><MIN value="{min}"/><MAX value="{max}"/></VALUES>
                </{param_or_field}>
            </TABLE>
        </RESOURCE>
    </VOTABLE>
    """

    for param_or_field in ["PARAM", "FIELD"]:
        testvals["param_or_field"] = param_or_field
        xml_string = minmax_template.format(**testvals)
        xml_bytes = xml_string.encode("utf-8")

        # Parse without exceptions.
        vot = parse(io.BytesIO(xml_bytes), verify="exception")

        # Assert that the scalar or array values were parsed and stored as expected.
        param_or_field = vot.get_field_by_id_or_name("param_or_field")
        _assert_minmax_match(param_or_field.values.min, testvals["expected_min"])
        _assert_minmax_match(param_or_field.values.max, testvals["expected_max"])


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
    for version in ("1.1", "1.2", "1.3", "1.4", "1.5"):
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
    parse(io.BytesIO(begin + b"1.5" + middle + b"1.3" + end), verify="exception")

    # Invalid versions
    for bversion in (b"1.0", b"2.0"):
        with pytest.warns(W21), pytest.warns(W41):
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
    assert 'xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 ' in xml
    assert 'http://www.ivoa.net/xml/VOTable/VOTable-1.3.xsd"' in xml

    xml = votable_xml_string("1.4")
    assert 'xmlns="http://www.ivoa.net/xml/VOTable/v1.3"' in xml
    assert 'xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 ' in xml
    assert 'http://www.ivoa.net/xml/VOTable/VOTable-1.4.xsd"' in xml

    xml = votable_xml_string("1.5")
    assert 'xmlns="http://www.ivoa.net/xml/VOTable/v1.3"' in xml
    assert 'xsi:schemaLocation="http://www.ivoa.net/xml/VOTable/v1.3 ' in xml
    assert 'http://www.ivoa.net/xml/VOTable/VOTable-1.5.xsd"' in xml


def _squash_xml(data):
    """
    Utility squashing XML fragment to easier their comparison
    This function is only used in the test module. It was more convenient for comparing the xml.
    """
    return data.replace(" ", "").replace("\n", "").replace('"', "").replace("'", "")


def test_mivot_constructor():
    """
    Construct a MIVOT block with wrong tag to test the expected exception
    """
    with pytest.raises(ValueError, match="not well-formed"):
        MivotBlock(
            """
            <VODML xmlns="http://www.ivoa.net/xml/mivot" >
               <REPORT status="OK">Unit test mivot block1</REPORT>
               <WRONG TAG>
               </GLOBALS>
            </VODML>
            """
        )


def test_mivot_readout():
    """
    Test the MIVOT block extraction from a file against a reference block stored in data
    """
    votable = parse(get_pkg_data_filename("data/mivot_annotated_table.xml"))

    ref_data = ""
    for resource in votable.resources:
        with open(
            get_pkg_data_filename("data/mivot_block_custom_datatype.xml")
        ) as reference:
            ref_data = reference.read()
            assert _squash_xml(ref_data) == _squash_xml(resource.mivot_block.content)
        assert len(resource.tables) == 1


def test_mivot_write():
    """
    Build a VOTable, put a MIVOT block in the first resource, checks it can be retrieved
    as well as the following table
    """

    mivot_block = MivotBlock(
        """
    <VODML xmlns="http://www.ivoa.net/xml/mivot" >
       <REPORT status="OK">
       Unit test mivot block1
       </REPORT>
       <GLOBALS>
       </GLOBALS>
    </VODML>
    """
    )
    vtf = VOTableFile()
    mivot_resource = Resource()
    mivot_resource.type = "meta"
    mivot_resource.mivot_block = mivot_block
    # pack the meta resource in a top level resource
    r1 = Resource()
    r1.type = "results"
    r1.resources.append(mivot_resource)
    vtf.resources.append(r1)
    # Push the VOTable in an IOSTream (emulates a disk saving)
    buff = io.BytesIO()
    vtf.to_xml(buff)

    # Read the IOStream (emulates a disk readout)
    buff.seek(0)
    vtf2 = parse(buff)
    assert len(vtf2.resources) == 1
    for resource in vtf2.resources:
        assert _squash_xml(mivot_block.content) == _squash_xml(
            resource.mivot_block.content
        )
        assert len(resource.tables) == 0


def test_mivot_write_after_table():
    """
    Build a VOTable, put a MIVOT block and a table in the first resource, checks it can be retrieved
    as well as the following table
    """

    mivot_block = MivotBlock(
        """
    <VODML xmlns="http://www.ivoa.net/xml/mivot" >
      <REPORT status="OK">Unit test mivot block1</REPORT>
      <GLOBALS>
      </GLOBALS>
    </VODML>
    """
    )
    vtf = VOTableFile()
    mivot_resource = Resource()
    mivot_resource.type = "meta"
    mivot_resource.mivot_block = mivot_block
    # pack the meta resource in a top level resource
    r1 = Resource()
    r1.type = "results"
    i1 = tree.Info(name="test_name", value="test_value")
    r1.infos.append(i1)
    r1.resources.append(mivot_resource)
    t1 = tree.TableElement(vtf)
    t1.name = "t1"
    r1.tables.append(t1)
    vtf.resources.append(r1)
    # Push the VOTable in an IOSTream (emulates a disk saving)
    buff = io.BytesIO()
    vtf.to_xml(buff)

    # Read the IOStream (emulates a disk readout)
    buff.seek(0)
    vtf2 = parse(buff)
    assert len(vtf2.resources) == 1
    for resource in vtf2.resources:
        assert _squash_xml(mivot_block.content) == _squash_xml(
            resource.mivot_block.content
        )
        assert len(resource.tables) == 1


def test_write_no_mivot():
    """
    Build a VOTable, put an empty MIVOT block in the first resource, checks it can be retrieved
    as well as the following table
    """

    vtf = VOTableFile()
    mivot_resource = Resource()
    mivot_resource.type = "meta"
    # pack the meta resource in a top level resource
    r1 = Resource()
    r1.type = "results"
    r1.resources.append(mivot_resource)
    t1 = tree.TableElement(vtf)
    t1.name = "t1"
    r1.tables.append(t1)
    vtf.resources.append(r1)
    # Push the VOTable in an IOSTream (emulates a disk saving)
    buff = io.BytesIO()
    vtf.to_xml(buff)

    # Read the IOStream (emulates a disk readout)
    buff.seek(0)
    vtf2 = parse(buff)
    assert len(vtf2.resources) == 1
    for resource in vtf2.resources:
        assert (
            _squash_xml(resource.mivot_block.content)
            == "<VODMLxmlns=http://www.ivoa.net/xml/mivot><REPORTstatus=KO>NoMivotblock</REPORT></VODML>"
        )
        assert len(resource.tables) == 1


def test_mivot_write_after_resource():
    """
    Build a VOTable, put a MIVOT block in the first resource after another meta resource,
    checks it can be retrieved as well as the following table
    """

    mivot_block = MivotBlock(
        """
    <VODML xmlns="http://www.ivoa.net/xml/mivot" >
      <REPORT status="OK">Unit test mivot block1</REPORT>
      <GLOBALS>
      </GLOBALS>
    </VODML>
    """
    )
    vtf = VOTableFile()
    mivot_resource = Resource()
    mivot_resource.type = "meta"
    mivot_resource.mivot_block = mivot_block
    # pack the meta resource in a top level resource
    r1 = Resource()
    r1.type = "results"
    i1 = tree.Info(name="test_name", value="test_value")
    r1.infos.append(i1)
    meta_resource = Resource()
    meta_resource.type = "meta"
    r1.resources.append(meta_resource)
    r1.resources.append(mivot_resource)
    t1 = tree.TableElement(vtf)
    t1.name = "t1"
    r1.tables.append(t1)
    vtf.resources.append(r1)
    # Push the VOTable in an IOSTream (emulates a disk saving)
    buff = io.BytesIO()
    vtf.to_xml(buff)

    # Read the IOStream (emulates a disk readout)
    buff.seek(0)
    vtf2 = parse(buff)
    assert len(vtf2.resources) == 1
    for resource in vtf2.resources:
        assert _squash_xml(mivot_block.content) == _squash_xml(
            resource.mivot_block.content
        )
        assert len(resource.tables) == 1


def test_mivot_forbidden_write():
    """
    Build a meta resource containing a MIVOT block,
    build the dummy MIVOT block first.
    """
    mivot_block = MivotBlock(
        """
    <VODML xmlns="http://www.ivoa.net/xml/mivot" >
       <REPORT status="KO">Unit test mivot block1</REPORT>
       <GLOBALS/>
    </VODML>
    """
    )
    # package the MIVOT block in the resource
    mivot_resource = Resource()
    mivot_resource.type = "results"

    with pytest.raises(E26):
        # A MIVOT block must be with "type=meta"
        mivot_resource.mivot_block = mivot_block


def test_mivot_order(tmp_path):
    """
    Build a VOTable with 2 resources containing MivotBlock, parse it, and write it in a file.
    Then compare it with another file to see if the order of the elements in a resource is respected,
    in particular the MivotBlock which should be before the tables.
    """

    mivot_block = MivotBlock(
        """
    <VODML xmlns="http://www.ivoa.net/xml/mivot" >
    </VODML>
    """
    )
    vtf = VOTableFile()

    mivot_resource = Resource()
    mivot_resource.type = "meta"
    mivot_resource.mivot_block = mivot_block

    mivot_resource2 = Resource()
    mivot_resource2.type = "meta"
    mivot_resource2.mivot_block = mivot_block

    # R1 : 2 mivot_block, 2 tables, 1 description, 1 info, 1 CooSys
    r1 = Resource()
    r1.type = "results"

    t1 = tree.TableElement(vtf)
    t1.name = "t1"
    t2 = tree.TableElement(vtf)
    t2.name = "t2"

    r1.tables.append(t1)
    r1.tables.append(t2)

    r1.resources.append(mivot_resource)
    r1.resources.append(mivot_resource2)

    cs = tree.CooSys(ID="_XYZ", system="ICRS")
    r1.coordinate_systems.append(cs)
    i1 = tree.Info(name="test_name", value="test_value")
    r1.infos.append(i1)

    vtf.resources.append(r1)

    # R2 : 1 resource "results", 1 mivot_block and 1 table
    r2 = Resource()
    r2.type = "results"

    r3 = Resource()
    r3.type = "results"

    t3 = tree.TableElement(vtf)
    t3.name = "t3"
    r2.tables.append(t3)
    r2.resources.append(mivot_resource)
    r2.resources.append(r3)

    vtf.resources.append(r2)

    # Push the VOTable in an IOSTream (emulates a disk saving)
    buff = io.BytesIO()
    vtf.to_xml(buff)

    # Read the IOStream (emulates a disk readout)
    buff.seek(0)
    vtf2 = parse(buff)

    vpath = get_pkg_data_filename("data/test.order.xml")
    vpath_out = str(tmp_path / "test.order.out.xml")
    vtf2.to_xml(vpath_out)

    # We want to remove the xml header from the VOTable
    with open(vpath_out) as file:
        lines = file.readlines()
    # The xml header is on 2 lines (line 2 and 3)
    del lines[1]
    del lines[1]

    with open(vpath_out, "w") as file:
        file.writelines(lines)

    assert filecmp.cmp(vpath, vpath_out)
