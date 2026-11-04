# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

from astropy.io.votable.validator import html


def test_write_source_line_keeps_xml_forbidden_chars():
    """
    write_source_line() of a line with a form feed at the highlighted
    position writes the form feed unchanged inside the highlight.  The
    report exists to show broken files, so it must not fail on the
    character it is pointing at.
    """
    out = io.StringIO()
    html.write_source_line(out, b"ab\x0ccd", 2)

    assert out.getvalue() == '  ab<span class="highlight">\x0c</span>cd\n\n'
