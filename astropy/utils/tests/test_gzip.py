from ...tests.helper import catch_warnings
from ..exceptions import AstropyDeprecationWarning
from ...extern.six.moves import reload_module


def test_import_warning():
    # this import *maybe* raises the warning too, but we force the matter below
    from ..compat import gzip  # pylint: disable=W0611

    with catch_warnings() as w:
        # this ensures the warning will be re-issues at reload.  Without this,
        # the warning is missing if you invoke astropy.test twice in the same
        # python session
        if hasattr(gzip, '__warningregistry__'):
            gzip.__warningregistry__.clear()
        reload_module(gzip)

    assert len(w) == 1
    assert w[0].category == AstropyDeprecationWarning
    assert w[0].message.args[0] == ("astropy.utils.compat.gzip is now deprecated "
                                    "- use the gzip module directly instead")
