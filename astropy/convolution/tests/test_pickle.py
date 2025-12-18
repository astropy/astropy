# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy import convolution as conv
from astropy.tests.helper import check_pickling_recovery, pickle_protocol  # noqa: F401


@pytest.mark.parametrize(
    "original",
    [
        conv.CustomKernel(array=np.random.rand(15)),
        # FIXME: Gaussian1DKernel sometimes fails check_pickling_recovery()
        #        because Gaussian1D (modeling) sometimes fails it.
        #        Fixing the latter will automatically fix the former.
        # conv.Gaussian1DKernel(1.0, x_size=5),
        conv.Gaussian2DKernel(1.0, x_size=5, y_size=5),
    ],
    ids=lambda x: type(x).__name__,
)
def test_simple_object(pickle_protocol, original):  # noqa: F811
    # Tests easily instantiated objects
    check_pickling_recovery(original, pickle_protocol)
