# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy import ma

from numpy.testing import assert_allclose

try:
    import matplotlib
    HAS_MATPLOTLIB = True
except:
    HAS_MATPLOTLIB = False

from ..mpl_normalize import ImageNormalize

from ...tests.helper import pytest
from ..stretch import SqrtStretch


@pytest.mark.skipif('HAS_MATPLOTLIB')
def test_error_message():
    with pytest.raises(ImportError) as exc:
        ImageNormalize()
    assert exc.value.args[0] == "matplotlib is required in order to use this class"


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_normalize_scalar():

    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=True)

    assert_allclose(n(6), 0.70710678)


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_normalize_clip():

    data = np.linspace(0., 15., 6)
    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=True)

    output = n(data)

    assert_allclose(output, [0., 0.35355339, 0.70710678, 0.93541435, 1., 1.])

    assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_normalize_noclip():

    data = np.linspace(0., 15., 6)
    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=False)

    output = n(data)

    assert_allclose(output, [np.nan, 0.35355339, 0.70710678, 0.93541435, 1.11803399, 1.27475488])

    assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])

    assert_allclose(n.inverse(n(data))[1:], data[1:])


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_normalize_implicit_autoscale():

    data = np.linspace(0., 15., 6)

    n = ImageNormalize(vmin=None, vmax=10., stretch=SqrtStretch(), clip=False)
    n(data)

    assert n.vmin == np.min(data)
    assert n.vmax == 10.

    n = ImageNormalize(vmin=2., vmax=None, stretch=SqrtStretch(), clip=False)
    n(data)

    assert n.vmin == 2.
    assert n.vmax == np.max(data)


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_masked_normalize_clip():

    data = np.linspace(0., 15., 6)
    mdata = ma.array(data, mask=[0, 0, 1, 0, 0, 0])

    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=True)

    output = n(mdata)

    assert_allclose(output.filled(-10), [0., 0.35355339, 1., 0.93541435, 1., 1.])
    assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_masked_normalize_noclip():

    data = np.linspace(0., 15., 6)
    mdata = ma.array(data, mask=[0, 0, 1, 0, 0, 0])
    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=False)

    output = n(mdata)

    assert_allclose(output.filled(-10), [np.nan, 0.35355339, -10, 0.93541435, 1.11803399, 1.27475488])
    assert_allclose(output.mask, [0, 0, 1, 0, 0, 0])

    assert_allclose(n.inverse(n(data))[1:], data[1:])
