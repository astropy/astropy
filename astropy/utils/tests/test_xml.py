# Licensed under a 3-clause BSD style license - see LICENSE.rst
import io

from ..xml import check, writer


def test_writer():
    fh = io.StringIO()
    w = writer.XMLWriter(fh)
    with w.tag("html"):
        with w.tag("body"):
            w.data("This is the content")
            w.comment("comment")

    value = ''.join(fh.getvalue().split())
    assert value == '<html><body>Thisisthecontent<!--comment--></body></html>'


def test_check_id():
    assert check.check_id("Fof32")
    assert check.check_id("_Fof32")
    assert not check.check_id("32Fof")


def test_fix_id():
    assert check.fix_id("Fof32") == "Fof32"
    assert check.fix_id("@#f") == "___f"


def test_check_token():
    assert check.check_token("token")
    assert not check.check_token("token\rtoken")


def test_check_mime_content_type():
    assert check.check_mime_content_type("image/jpeg")
    assert not check.check_mime_content_type("image")


def test_check_anyuri():
    assert check.check_anyuri("https://github.com/astropy/astropy")
