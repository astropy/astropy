from ...tests.helper import catch_warnings
from ..exceptions import AstropyDeprecationWarning


def test_import_warning():

    with catch_warnings() as w:
        from ..compat import argparse

    assert len(w) == 1
    assert w[0].category == AstropyDeprecationWarning
    assert w[0].message.args[0] == ("astropy.utils.compat.argparse is now "
                                    "deprecated - use the argparse module "
                                    "directly instead")
