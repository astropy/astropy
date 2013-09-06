# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools

import numpy as np

from ...tests.helper import pytest
from .. import Time
from ...utils.iers import iers  # used in testing

allclose_jd = functools.partial(np.allclose, rtol=1e-15, atol=0)
allclose_sec = functools.partial(np.allclose, rtol=1e-15, atol=1e-9)
# 1 nanosec atol

try:
    iers.IERS_A.open()  # check if IERS_A is available
except IOError:
    HAS_IERS_A = False
else:
    HAS_IERS_A = True


class TestTimeUT1():
    """Test Time.ut1 using IERS tables"""

    def test_ut1(self):
        t = Time(['2012-06-30 12:00:00', '2012-06-30 23:59:59',
                  '2012-06-30 23:59:60', '2012-07-01 00:00:00',
                  '2012-07-01 12:00:00'], scale='utc')
        t_ut1_jd = t.ut1.jd
        t_comp = np.array([2456108.9999932079,
                           2456109.4999816339,
                           2456109.4999932083,
                           2456109.5000047823,
                           2456110.0000047833])
        assert allclose_jd(t_ut1_jd, t_comp)
        tnow = Time.now()
        with pytest.raises(IndexError):
            tnow.ut1


@pytest.mark.skipif('not HAS_IERS_A')
class TestTimeUT1_IERSA():
    def test_ut1_iers_A(self):
        tnow = Time.now()
        iers_a = iers.IERS_A.open()
        tnow.delta_ut1_utc, status = iers_a.ut1_utc(tnow, return_status=True)
        assert status == iers.FROM_IERS_A_PREDICTION
        tnow_ut1_jd = tnow.ut1.jd
        assert tnow_ut1_jd != tnow.jd
