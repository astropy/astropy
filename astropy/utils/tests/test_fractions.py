from ...tests.helper import catch_warnings
from ..exceptions import AstropyDeprecationWarning


def test_import_warning():

    with catch_warnings() as w:
        from ..compat import fractions

    assert len(w) == 1
    assert w[0].category == AstropyDeprecationWarning
    assert w[0].message.args[0] == ("astropy.utils.compat.fractions is now "
                                    "deprecated - use the fractions module "
                                    "directly instead")
