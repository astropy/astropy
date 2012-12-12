# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains tests for the name resolve convenience module.
"""

# Standard library

# Third party
import numpy as np
import pytest

# Astropy
from ..name_resolve import get_icrs_coordinates

def test_names():

    with pytest.raises(ValueError):
        get_icrs_coordinates("m87h34hhh")

    for name in ["ngc 3642", "m42", "castor", "pollux"]:
        print(get_icrs_coordinates(name))

    assert False
