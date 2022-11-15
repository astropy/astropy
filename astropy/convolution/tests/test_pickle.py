# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy import convolution as conv
from astropy.tests.helper import check_pickling_recovery, pickle_protocol  # noqa: F401


@pytest.mark.parametrize(
    ("name", "args", "kwargs", "xfail"),
    [
        (conv.CustomKernel, [], {"array": np.random.rand(15)}, False),
        (conv.Gaussian1DKernel, [1.0], {"x_size": 5}, True),
        (conv.Gaussian2DKernel, [1.0], {"x_size": 5, "y_size": 5}, True),
    ],
)
def test_simple_object(pickle_protocol, name, args, kwargs, xfail):  # noqa: F811
    # Tests easily instantiated objects
    if xfail:
        pytest.xfail()
    original = name(*args, **kwargs)
    check_pickling_recovery(original, pickle_protocol)
