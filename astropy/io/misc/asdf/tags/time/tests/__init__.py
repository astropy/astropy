import pytest

from astropy.io.misc.asdf.tests import ASDF_ENTRY_INSTALLED

if not ASDF_ENTRY_INSTALLED:
    pytest.skip(
        "The astropy asdf entry points are not installed", allow_module_level=True
    )
