# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
from ...extern.six.moves import zip

# THIRD-PARTY
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

# LOCAL
from .. import models
from ..utils import merge_sampleset


def test_no_sampleset():
    m = models.Const1D()
    assert m.sampleset() is None


def test_gaussian1d_sampleset():
    m = models.Gaussian1D(mean=5000, stddev=10)
    ans1 = np.arange(4945, 5056)
    assert_array_equal(m.sampleset(), ans1)

    # n_models > 1
    m = models.Gaussian1D(amplitude=[1, 1], mean=[5000, 6000], stddev=[10, 100],
                          n_models=2)
    ans = [ans1, np.arange(5450, 6560, 10)]
    assert_array_equal(m.sampleset(), ans)


def test_lorentz1d_sampleset():
    fwhm_vals = [11.774100225154747, 117.74100225154747]  # stddev = 5, 50

    m = models.Lorentz1D(x_0=5000, fwhm=fwhm_vals[0])  # stddev = 5
    ans1 = np.arange(4875, 5125.5, 0.5)
    assert_allclose(m.sampleset(), ans1)

    # n_models > 1
    m = models.Lorentz1D(amplitude=[1, 1], x_0=[5000, 6000], fwhm=fwhm_vals,
                         n_models=2)
    ans = [ans1, np.arange(4750, 7255, 5)]
    assert_allclose(m.sampleset(), ans)


def test_box1d_sampleset():
    m = models.Box1D(x_0=5000, width=10)
    ans1 = np.arange(4995 - 0.01, 5005 + 0.02, 0.01)
    assert_array_equal(m.sampleset(), ans1)

    # Minimal points only
    ans = [4995 - 0.01, 4995, 5005, 5005 + 0.01]
    assert_array_equal(m.sampleset(minimal=True), ans)

    # n_models > 1
    m = models.Box1D(amplitude=[1, 1], x_0=[5000, 6000], width=[10, 100],
                     n_models=2)
    x = m.sampleset()
    ans = [ans1, np.arange(5950 - 0.01, 6050 + 0.02, 0.01)]
    # Each answer has different len, so has to be compared separately
    for a, b in zip(x, ans):
        assert_array_equal(a, b)


def test_trapezoid1d_sampleset():
    m = models.Trapezoid1D(amplitude=100, x_0=5000, width=10, slope=1)
    ans1 = [4895, 4995, 5005, 5105]
    assert_array_equal(m.sampleset(), ans1)

    # n_models > 1
    m = models.Trapezoid1D(amplitude=[100, 100], x_0=[5000, 6000],
                           width=[10, 100], slope=[1, 2], n_models=2)
    ans = [ans1, [5900.0, 5950.0, 6050.0, 6100.0]]
    assert_array_equal(m.sampleset(), ans)


def test_mexicanhat1d_sampleset():
    m = models.MexicanHat1D(x_0=5000, sigma=10)
    ans1 = np.arange(4925, 5076, 1)
    assert_array_equal(m.sampleset(), ans1)

    # n_models > 1
    m = models.MexicanHat1D(amplitude=[1, 1], x_0=[5000, 6000], sigma=[10, 100],
                            n_models=2)
    ans = [ans1, np.arange(5250, 6760, 10)]
    assert_array_equal(m.sampleset(), ans)


def test_merge_none():
    m1 = models.PowerLaw1D(x_0=5000, alpha=40)
    m2 = models.Box1D(x_0=5000, width=10)
    ans = m2.sampleset()

    # Left has none
    m = m1 + m2
    assert_array_equal(m.sampleset(), ans)

    # Right has none
    m = m2 + m1
    assert_array_equal(m.sampleset(), ans)

    # Both has none
    m = m1 + m1
    assert m.sampleset() is None


def test_merge_array_only():
    """Test the underlying merging algorithm without models."""
    s1 = [5000, 5000.01, 5000.02, 5000.03, 5000.04, 6000]

    # Make sure "too close together" is taken care of properly
    thres = 1e-12
    s2 = [5000.005, 5000.02 + thres, 5500, 6000]
    ans = [5000, 5000.005, 5000.01, 5000.02, 5000.03, 5000.04, 5500, 6000]
    s = merge_sampleset(s1, s2, threshold=thres)
    ds = s[1:] - s[:-1]
    assert_allclose(s, ans)
    assert np.all(ds > thres)

    # Merging two of the same should not change anything
    s = merge_sampleset(s1, s1)
    assert_array_equal(s, s1)


def test_model_mul_add():
    """Test a normal model operation."""
    m1 = models.Box1D(x_0=1000, width=10)
    m2 = models.Trapezoid1D(x_0=2000, width=100, slope=0.1)
    m3 = models.Const1D(amplitude=5)
    m = (m1 + m2) * m3
    x = m.sampleset()

    # Make sure the result samples important points on Box1D
    s1 = np.hstack([x[:2], x[-6:-4]])  # First and last 2 points of the box
    ans1 = m1.sampleset(minimal=True)
    assert_allclose(s1, ans1)

    # Make sure the result also samples Trapezoid1D properly
    s2 = x[-4:]
    ans2 = m2.sampleset()
    assert_array_equal(s2, ans2)


def test_model_composition():
    """Test | operation on models."""
    g = models.Gaussian1D(mean=5000, stddev=23.548200450309494)  # FWHM = 10
    gx = g.sampleset()
    assert_allclose(g.fwhm, 10)

    # Redshift
    rs = models.RedshiftScaleFactor(z=1.3)
    m = rs.inverse | g  # Shifting wavelength requires the inverse
    s = m.sampleset()
    ans = rs(gx)
    assert_allclose(s, ans)

    # Shift
    sh = models.Shift(offset=100)
    m = sh.inverse | g  # Same idea as redshift
    s = m.sampleset()
    ans = sh(gx)
    assert_allclose(s, ans)

    # Scale
    sc = models.Scale(factor=0.5)
    m = sc.inverse | g  # Same idea as redshift
    s = m.sampleset()
    ans = sc(gx)
    assert_allclose(s, ans)

    # Shift after evaluation should not affect sampleset
    m = g | sh
    assert_array_equal(m.sampleset(), gx)

    # Shift | scale produces nothing
    m = sh | sc
    assert m.sampleset() is None

    # Knows nothing
    g2 = models.Gaussian1D(mean=6000, stddev=10)
    m = g | g2
    assert m.sampleset() is None


def test_multi_n_models():
    """This is currently not supported, so make sure it does not
    generate the wrong results."""
    m1 = models.Box1D(amplitude=[1, 1], x_0=[5000, 6000], width=[10, 100],
                          n_models=2)
    m2 = models.Trapezoid1D(amplitude=[100, 100], x_0=[5000, 6000],
                            width=[10, 100], slope=[1, 2], n_models=2)
    m = m1 + m2
    assert m.sampleset() is None
