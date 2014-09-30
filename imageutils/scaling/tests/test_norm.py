import numpy as np
from numpy import ma

from numpy.testing import assert_allclose

from ..normalize import ImageNormalize
from ..stretch import SqrtStretch


def test_normalize_scalar():

    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=True)

    assert_allclose(n(6), 0.70710678)


def test_normalize_clip():

    data = np.linspace(0., 15., 6)
    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=True)

    output = n(data)

    assert_allclose(output, [0., 0.35355339, 0.70710678, 0.93541435, 1., 1.])

    assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])


def test_normalize_noclip():

    data = np.linspace(0., 15., 6)
    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=False)

    output = n(data)

    assert_allclose(output, [np.nan, 0.35355339, 0.70710678, 0.93541435, 1.11803399, 1.27475488])

    assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])

    assert_allclose(n.inverse(n(data))[1:], data[1:])


def test_masked_normalize_clip():

    data = np.linspace(0., 15., 6)
    mdata = ma.array(data, mask=[0, 0, 1, 0, 0, 0])

    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=True)

    output = n(mdata)

    assert_allclose(output.filled(-10), [0., 0.35355339, 1., 0.93541435, 1., 1.])
    assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])


def test_masked_normalize_noclip():

    data = np.linspace(0., 15., 6)
    mdata = ma.array(data, mask=[0, 0, 1, 0, 0, 0])
    n = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(), clip=False)

    output = n(mdata)

    assert_allclose(output.filled(-10), [np.nan, 0.35355339, -10, 0.93541435, 1.11803399, 1.27475488])
    assert_allclose(output.mask, [0, 0, 1, 0, 0, 0])

    assert_allclose(n.inverse(n(data))[1:], data[1:])
