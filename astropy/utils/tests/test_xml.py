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
