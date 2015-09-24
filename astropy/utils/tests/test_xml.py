# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six

import io

from ..xml import check, unescaper, writer


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


def test_unescape_all():
    # str
    url_in = 'http://casu.ast.cam.ac.uk/ag/iphas-dsa/SubmitCone?' \
             'DSACAT=IDR&amp;amp;DSATAB=Emitters&amp;amp;'
    url_out = 'http://casu.ast.cam.ac.uk/ag/iphas-dsa/SubmitCone?' \
              'DSACAT=IDR&DSATAB=Emitters&'
    assert unescaper.unescape_all(url_in) == url_out

    # bytes
    url_in = b'http://casu.ast.cam.ac.uk/ag/iphas-dsa/SubmitCone?' \
             b'DSACAT=IDR&amp;amp;DSATAB=Emitters&amp;amp;'
    url_out = b'http://casu.ast.cam.ac.uk/ag/iphas-dsa/SubmitCone?' \
              b'DSACAT=IDR&DSATAB=Emitters&'
    assert unescaper.unescape_all(url_in) == url_out


def test_escape_xml():
    s = writer.xml_escape('This & That')
    assert type(s) == six.text_type
    assert s == 'This &amp; That'

    s = writer.xml_escape(1)
    assert type(s) == str
    assert s == '1'

    s = writer.xml_escape(b'This & That')
    assert type(s) == bytes
    assert s == b'This &amp; That'
