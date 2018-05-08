# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Some might be indirectly tested already in ``astropy.io.fits.tests``.
"""
import io

import numpy as np

from ..diff import report_diff_values
from ...table import Table


def test_float_comparison():
    """
    Regression test for https://github.com/spacetelescope/PyFITS/issues/21
    """
    f = io.StringIO()

    a = np.float32(0.029751372)
    b = np.float32(0.029751368)

    report_diff_values(f, a, b)
    out = f.getvalue()

    # This test doesn't care about what the exact output is, just that it
    # did show a difference in their text representations
    assert 'a>' in out
    assert 'b>' in out


def test_diff_types():
    """
    Regression test for https://github.com/astropy/astropy/issues/4122
    """
    f = io.StringIO()

    a = 1.0
    b = '1.0'

    report_diff_values(f, a, b)
    out = f.getvalue()

    assert out.lstrip('u') == ("  (float) a> 1.0\n    (str) b> '1.0'\n"
                               "           ? +   +\n")


def test_tablediff():
    """
    Test diff-ing two simple Table objects.
    """
    a = Table.read("""name    obs_date    mag_b  mag_v
M31     2012-01-02  17.0   16.0
M82     2012-10-29  16.2   15.2
M101    2012-10-31  15.1   15.5""", format='ascii')
    b = Table.read("""name    obs_date    mag_b  mag_v
M31     2012-01-02  17.0   16.5
M82     2012-10-29  16.2   15.2
M101    2012-10-30  15.1   15.5
NEW     2018-05-08   nan    9.0""", format='ascii')
    f = io.StringIO()
    report_diff_values(f, a, b)
    out = f.getvalue()
    assert out == ('     name  obs_date  mag_b mag_v\n'
                   '     ---- ---------- ----- -----\n'
                   '  a>  M31 2012-01-02  17.0  16.0\n'
                   '   ?                           ^\n'
                   '  b>  M31 2012-01-02  17.0  16.5\n'
                   '   ?                           ^\n'
                   '      M82 2012-10-29  16.2  15.2\n'
                   '  a> M101 2012-10-31  15.1  15.5\n'
                   '   ?               ^\n'
                   '  b> M101 2012-10-30  15.1  15.5\n'
                   '   ?               ^\n'
                   '  b>  NEW 2018-05-08   nan   9.0\n')
