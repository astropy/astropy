# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy import ma
from numpy.testing import assert_allclose
from ...tests.helper import pytest
from ..mpl_normalize import ImageNormalize
from ..interval import ManualInterval
from ..stretch import SqrtStretch

try:
    import matplotlib    # pylint: disable=W0611
    HAS_MATPLOTLIB = True
except:
    HAS_MATPLOTLIB = False


DATA = np.linspace(0., 15., 6)


@pytest.mark.skipif('HAS_MATPLOTLIB')
def test_normalize_error_message():
    with pytest.raises(ImportError) as exc:
        ImageNormalize()
    assert (exc.value.args[0] == "matplotlib is required in order to use "
            "this class")


@pytest.mark.skipif('not HAS_MATPLOTLIB')
class TestNormalize(object):

    def test_scalar(self):
        norm = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(),
                              clip=True)
        norm2 = ImageNormalize(data=6, interval=ManualInterval(2, 10),
                               stretch=SqrtStretch(), clip=True)
        assert_allclose(norm(6), 0.70710678)
        assert_allclose(norm(6), norm2(6))

    def test_clip(self):
        norm = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(),
                              clip=True)
        norm2 = ImageNormalize(DATA, interval=ManualInterval(2, 10),
                               stretch=SqrtStretch(), clip=True)
        output = norm(DATA)
        expected = [0., 0.35355339, 0.70710678, 0.93541435, 1., 1.]
        assert_allclose(output, expected)
        assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])
        assert_allclose(output, norm2(DATA))

    def test_noclip(self):
        norm = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(),
                              clip=False)
        norm2 = ImageNormalize(DATA, interval=ManualInterval(2, 10),
                               stretch=SqrtStretch(), clip=False)
        output = norm(DATA)
        expected = [np.nan, 0.35355339, 0.70710678, 0.93541435, 1.11803399,
                    1.27475488]
        assert_allclose(output, expected)
        assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])
        assert_allclose(norm.inverse(norm(DATA))[1:], DATA[1:])
        assert_allclose(output, norm2(DATA))

    def test_implicit_autoscale(self):
        norm = ImageNormalize(vmin=None, vmax=10., stretch=SqrtStretch(),
                              clip=False)
        norm2 = ImageNormalize(DATA, interval=ManualInterval(None, 10),
                               stretch=SqrtStretch(), clip=False)
        output = norm(DATA)
        assert norm.vmin == np.min(DATA)
        assert norm.vmax == 10.
        assert_allclose(output, norm2(DATA))

        norm = ImageNormalize(vmin=2., vmax=None, stretch=SqrtStretch(),
                              clip=False)
        norm2 = ImageNormalize(DATA, interval=ManualInterval(2, None),
                               stretch=SqrtStretch(), clip=False)
        output = norm(DATA)
        assert norm.vmin == 2.
        assert norm.vmax == np.max(DATA)
        assert_allclose(output, norm2(DATA))

    def test_masked_clip(self):
        mdata = ma.array(DATA, mask=[0, 0, 1, 0, 0, 0])
        norm = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(),
                              clip=True)
        norm2 = ImageNormalize(mdata, interval=ManualInterval(2, 10),
                               stretch=SqrtStretch(), clip=True)
        output = norm(mdata)
        expected = [0., 0.35355339, 1., 0.93541435, 1., 1.]
        assert_allclose(output.filled(-10), expected)
        assert_allclose(output.mask, [0, 0, 0, 0, 0, 0])
        assert_allclose(output, norm2(mdata))

    def test_masked_noclip(self):
        mdata = ma.array(DATA, mask=[0, 0, 1, 0, 0, 0])
        norm = ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch(),
                              clip=False)
        norm2 = ImageNormalize(mdata, interval=ManualInterval(2, 10),
                               stretch=SqrtStretch(), clip=False)
        output = norm(mdata)
        expected = [np.nan, 0.35355339, -10, 0.93541435, 1.11803399,
                    1.27475488]
        assert_allclose(output.filled(-10), expected)
        assert_allclose(output.mask, [0, 0, 1, 0, 0, 0])

        assert_allclose(norm.inverse(norm(DATA))[1:], DATA[1:])
        assert_allclose(output, norm2(mdata))
