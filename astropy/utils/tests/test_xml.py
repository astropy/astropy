# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

import pytest

from astropy.utils.compat.optional_deps import HAS_BLEACH
from astropy.utils.xml import check, iterparser, unescaper, writer


def test_writer():
    fh = io.StringIO()
    w = writer.XMLWriter(fh)
    with w.tag("html"):
        with w.tag("body"):
            w.data("This is the content")
            w.comment("comment")

    value = "".join(fh.getvalue().split())
    assert value == "<html><body>Thisisthecontent<!--comment--></body></html>"


# The C tokenizer returns one cached string for the element name "TD", which is
# by far the most common name in a VOTable, instead of building a new one per
# cell.  A check that gets this wrong either hands back "TD" for a name that
# only starts like it, such as TDX, or stops matching TD at all and quietly
# builds a new string for every cell.
TD_NEIGHBOUR_DOC = (
    b"<VOTABLE>"
    b"<TR><TD>1</TD><TD>2</TD></TR>"
    b"<T>short</T>"
    b"<TDX>longer</TDX>"
    b"<D>other</D>"
    b"<TD>3</TD>"
    b"</VOTABLE>"
)


def test_iterparser_reports_names_around_the_td_shortcut():
    """Make sure every start event reports the element name as it was
    written, including the names either side of each TD.
    """
    with iterparser.get_xml_iterator(io.BytesIO(TD_NEIGHBOUR_DOC)) as it:
        names = [tag for start, tag, data, pos in it if start]

    assert names == [
        "VOTABLE",
        "TR",
        "TD",
        "TD",
        "T",
        "TDX",
        "D",
        "TD",
    ]


def test_iterparser_returns_one_cached_string_for_every_td():
    """Assert the tag of every TD event is the same string object."""
    with iterparser.get_xml_iterator(io.BytesIO(TD_NEIGHBOUR_DOC)) as it:
        td_tags = [tag for start, tag, data, pos in it if tag == "TD"]

    assert len(td_tags) == 6
    assert all(tag is td_tags[0] for tag in td_tags)


def test_check_id():
    assert check.check_id("Fof32")
    assert check.check_id("_Fof32")
    assert not check.check_id("32Fof")


def test_fix_id():
    assert check.fix_id("Fof32") == "Fof32"
    assert check.fix_id("@#f") == "___f"


def test_check_token():
    assert check.check_token("token")
    assert check.check_token("t o k e n")
    assert not check.check_token("t  oken")
    assert not check.check_token("token\rtoken")
    assert not check.check_token(" foo")
    assert not check.check_token("bar\n")


def test_check_mime_content_type():
    assert check.check_mime_content_type("image/jpeg")
    assert not check.check_mime_content_type("image")


def test_check_anyuri():
    assert check.check_anyuri("https://github.com/astropy/astropy")


def test_unescape_all():
    # str
    url_in = (
        "http://casu.ast.cam.ac.uk/ag/iphas-dsa%2FSubmitCone?"
        "DSACAT=IDR&amp;amp;DSATAB=Emitters&amp;amp;"
    )
    url_out = (
        "http://casu.ast.cam.ac.uk/ag/iphas-dsa/SubmitCone?DSACAT=IDR&DSATAB=Emitters&"
    )
    assert unescaper.unescape_all(url_in) == url_out

    # bytes
    url_in = (
        b"http://casu.ast.cam.ac.uk/ag/iphas-dsa%2FSubmitCone?"
        b"DSACAT=IDR&amp;amp;DSATAB=Emitters&amp;amp;"
    )
    url_out = (
        b"http://casu.ast.cam.ac.uk/ag/iphas-dsa/SubmitCone?DSACAT=IDR&DSATAB=Emitters&"
    )
    assert unescaper.unescape_all(url_in) == url_out


def test_escape_xml():
    s = writer.xml_escape("This & That")
    assert type(s) == str
    assert s == "This &amp; That"

    s = writer.xml_escape(1)
    assert type(s) == str
    assert s == "1"

    s = writer.xml_escape(b"This & That")
    assert type(s) == bytes
    assert s == b"This &amp; That"


@pytest.mark.skipif(HAS_BLEACH, reason="bleach is installed")
def test_escape_xml_without_bleach():
    fh = io.StringIO()
    w = writer.XMLWriter(fh)

    with pytest.raises(
        ValueError, match=r"bleach package is required when HTML escaping is disabled"
    ):
        with w.xml_cleaning_method("bleach_clean"):
            pass


@pytest.mark.skipif(not HAS_BLEACH, reason="requires bleach")
def test_escape_xml_with_bleach():
    fh = io.StringIO()
    w = writer.XMLWriter(fh)

    # Turn off XML escaping, but still sanitize unsafe tags like <script>
    with w.xml_cleaning_method("bleach_clean"):
        w.start("td")
        w.data("<script>x</script> <em>OK</em>")
        w.end(indent=False)
    assert fh.getvalue() == "<td>&lt;script&gt;x&lt;/script&gt; <em>OK</em></td>\n"

    fh = io.StringIO()
    w = writer.XMLWriter(fh)

    # Default is True (all XML tags escaped)
    with w.xml_cleaning_method():
        w.start("td")
        w.data("<script>x</script> <em>OK</em>")
        w.end(indent=False)
    assert (
        fh.getvalue()
        == "<td>&lt;script&gt;x&lt;/script&gt; &lt;em&gt;OK&lt;/em&gt;</td>\n"
    )
