# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy import constants as const
from astropy.tests.helper import check_pickling_recovery, pickle_protocol  # noqa: F401

originals = [
    const.Constant("h_fake", "Not Planck", 0.0, "J s", 0.0, "fakeref", system="si"),
    const.h,
    const.e.si,
]


@pytest.mark.parametrize("original", originals)
def test_new_constant(pickle_protocol, original):  # noqa: F811
    check_pickling_recovery(original, pickle_protocol)
