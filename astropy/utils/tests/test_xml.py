# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

import pytest

from astropy.utils.compat.optional_deps import HAS_BLEACH, HAS_LXML
from astropy.utils.xml import check, unescaper, writer


def test_writer():
    fh = io.StringIO()
    w = writer.XMLWriter(fh)
    with w.tag("html"):
        with w.tag("body"):
            w.data("This is the content")
            w.comment("comment")

    value = "".join(fh.getvalue().split())
    assert value == "<html><body>Thisisthecontent<!--comment--></body></html>"


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


# Characters XML does not permit anywhere in a document, See
#     https://www.w3.org/TR/xml/#NT-Char
# The control characters and U+FFFE / U+FFFF are listed in full.  The
# surrogate block U+D800-U+DFFF has 2048 members, so we test a subsample.
FORBIDDEN_XML_CHARS = [
    chr(c)
    for c in [
        *range(0x09),
        0x0B,
        0x0C,
        *range(0x0E, 0x20),
        0xD800,
        0xDBFF,
        0xDC00,
        0xDFFF,
        0xFFFE,
        0xFFFF,
    ]
]

# Characters that look like they might be forbidden but are not.
PERMITTED_XML_CHARS = [
    "\t",
    "\n",
    "\r",
    " ",
    "\x7f",
    "\xa0",
    "\u2028",
    "\ufffd",
    "\U0001f600",
]


@pytest.mark.parametrize("char", FORBIDDEN_XML_CHARS)
def test_escape_xml_rejects_forbidden_chars(char):
    # XML does not allow these, so the writer refuses them.  They cannot be
    # escaped around either: they have no escaped form.
    with pytest.raises(ValueError, match="XML does not permit this character"):
        writer.xml_escape(f"a{char}b")

    with pytest.raises(ValueError, match="XML does not permit this character"):
        writer.xml_escape_cdata(f"a{char}b")


@pytest.mark.parametrize("char", PERMITTED_XML_CHARS)
def test_escape_xml_keeps_permitted_chars(char):
    assert writer.xml_escape(f"a{char}b") == f"a{char}b"
    assert writer.xml_escape_cdata(f"a{char}b") == f"a{char}b"


def test_escape_xml_does_not_drop_nulls():
    with pytest.raises(ValueError, match=r"cannot write U\+0000 at position 3"):
        writer.xml_escape("NGC\x001068")


def test_escape_xml_escapes_but_does_not_check_bytes():
    # These functions are for text.  Bytes are still accepted and escaped,
    # but their encoding is not known here, so they are not checked against
    # the characters XML allows.
    assert writer.xml_escape(b"a<b") == b"a&lt;b"
    assert writer.xml_escape(b"a\x0cb") == b"a\x0cb"
    assert writer.xml_escape(b"a\x00b") == b"a\x00b"


def test_forbidden_char_position_is_in_characters():
    # The position is counted in characters, not bytes.  A byte offset into
    # the UTF-8 encoding would point somewhere else entirely once the text is
    # not ASCII, and the caller is looking at a str.
    e_acute = "\xe9"  # one character, but two bytes in UTF-8
    form_feed = "\x0c"  # the character that should be rejected
    text = e_acute * 3 + form_feed

    assert text.index(form_feed) == 3

    with pytest.raises(ValueError, match=r"cannot write U\+000C at position 3"):
        writer.xml_escape(text)


def test_writer_rejects_forbidden_chars_in_text():
    # Text nodes are wrapped with textwrap.fill() before they are escaped, and
    # wrapping would turn a form feed into a space, so the check has to happen
    # before that.
    fh = io.StringIO()
    w = writer.XMLWriter(fh)
    with pytest.raises(ValueError, match="in text of element 'description'"):
        with w.tag("description"):
            w.data("the\x0ctarget")


def test_writer_rejects_forbidden_chars_in_attribute():
    fh = io.StringIO()
    w = writer.XMLWriter(fh)
    with pytest.raises(ValueError, match="in attribute 'name' of element 'field'"):
        w.start("field", {"name": "the\x0ctarget"})


def test_writer_rejects_forbidden_chars_without_escaping():
    # xml_cleaning_method("none") turns off escaping, but a character XML
    # cannot represent must still be refused.
    fh = io.StringIO()
    w = writer.XMLWriter(fh)
    with pytest.raises(ValueError, match="XML does not permit this character"):
        with w.xml_cleaning_method("none"):
            with w.tag("td"):
                w.data("the\x0ctarget")


@pytest.mark.skipif(not HAS_LXML, reason="requires lxml")
@pytest.mark.parametrize("char", FORBIDDEN_XML_CHARS)
def test_lxml_also_rejects_forbidden_chars(char):
    # Cross-check our XML character implementation against the lxml library's
    # implementation, which enforces the same XML 1.0 Char standard.
    from lxml import etree

    with pytest.raises(ValueError):
        etree.Element("a").set("n", f"a{char}b")


@pytest.mark.skipif(not HAS_LXML, reason="requires lxml")
@pytest.mark.parametrize("char", PERMITTED_XML_CHARS)
def test_lxml_also_accepts_permitted_chars(char):
    from lxml import etree

    element = etree.Element("a")
    element.set("n", f"a{char}b")
    assert element.get("n") == f"a{char}b"


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
