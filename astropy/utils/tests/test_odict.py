from ...tests.helper import catch_warnings
from ..exceptions import AstropyDeprecationWarning


def test_init_warning():
    
    from ..compat.odict import OrderedDict
    
    with catch_warnings() as w:
        OrderedDict()
    
    assert len(w) == 1
    assert w[0].category == AstropyDeprecationWarning
    assert w[0].message.args[0] == ("astropy.utils.compat.odict.OrderedDict is "
                                    "now deprecated - import OrderedDict from "
                                    "the collections module instead")
    