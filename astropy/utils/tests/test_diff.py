# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Some might be indirectly tested already in ``astropy.io.fits.tests``.
"""
import io

import numpy as np
import pytest

from ..diff import diff_values, report_diff_values, where_not_allclose
from ...table import Table


@pytest.mark.parametrize('a', [np.nan, np.inf, 1.11, 1, 'a'])
def test_diff_values_false(a):
    assert not diff_values(a, a)


@pytest.mark.parametrize(
    ('a', 'b'),
    [(np.inf, np.nan), (1.11, 1.1), (1, 2), (1, 'a'), ('a', 'b')])
def test_diff_values_true(a, b):
    assert diff_values(a, b)


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


@pytest.mark.parametrize('kwargs', [{}, {'atol': 0, 'rtol': 0}])
def test_where_not_allclose(kwargs):
    a = np.array([1, np.nan, np.inf, 4.5])
    b = np.array([1, np.inf, np.nan, 4.6])

    assert where_not_allclose(a, b, **kwargs) == ([3], )
    assert len(where_not_allclose(a, a, **kwargs)[0]) == 0
