# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from io import StringIO

import numpy as np

from .. import velocities
from ...tests.helper import pytest
from ... import units as u

try:
    import scipy  # pylint: disable=W0611
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


@pytest.mark.skipif('not HAS_SCIPY')
def test_units():
    """ Test if the right units are being returned"""

    z1=2.0
    z2=2.01
    assert velocities.v_from_z(z1,z2).unit == u.Unit('km/s)')


@pytest.mark.skipif('not HAS_SCIPY')
def test_values():
    """ Test the functions 
    """

    # Test values were taken from a match to XIDL
    z1=2.0
    z2=2.01
    assert np.allclose(velocities.v_from_z(z1,z2).value,
                       [-997.8, -997.643, -997.4], rtol=1e-3)


