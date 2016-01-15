from ...tests.helper import catch_warnings
from ..exceptions import AstropyDeprecationWarning


def test_import_warning():
    with catch_warnings() as w:
        from ..compat import subprocess

    assert len(w) == 1
    assert w[0].category == AstropyDeprecationWarning
    assert w[0].message.args[0] == ("astropy.utils.compat.subprocess is now deprecated "
                                    "- use the Python subprocess module directly instead")
