# Licensed under a 3-clause BSD style license - see LICENSE.rst

from packaging.version import Version
import pytest
import numpy as np
from numpy import ma
from numpy.testing import assert_allclose, assert_equal

from astropy.visualization.mpl_normalize import ImageNormalize, simple_norm, imshow_norm
from astropy.visualization.interval import ManualInterval, PercentileInterval
from astropy.visualization.stretch import LogStretch, PowerStretch, SqrtStretch
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB, HAS_PLT  # noqa

if HAS_MATPLOTLIB:
    import matplotlib
    MATPLOTLIB_LT_32 = Version(matplotlib.__version__) < Version('3.2')

DATA = np.linspace(0., 15., 6)
DATA2 = np.arange(3)
DATA2SCL = 0.5 * DATA2
DATA3 = np.linspace(-3., 3., 7)
STRETCHES = (SqrtStretch(), PowerStretch(0.5), LogStretch())
INVALID = (None, -np.inf, -1)


@pytest.mark.skipif('HAS_MATPLOTLIB')
def test_normalize_error_message():
    with pytest.raises(ImportError) as exc:
        ImageNormalize()
    assert (exc.value.args[0] == "matplotlib is required in order to use "
            "this class.")


@pytest.mark.skipif('not HAS_MATPLOTLIB')
class TestNormalize:
    def test_invalid_interval(self):
        with pytest.raises(TypeError):
            ImageNormalize(vmin=2., vmax=10., interval=ManualInterval,
                           clip=True)

    def test_invalid_stretch(self):
        with pytest.raises(TypeError):
            ImageNormalize(vmin=2., vmax=10., stretch=SqrtStretch,
                           clip=True)

    def test_stretch_none(self):
        with pytest.raises(ValueError):
            ImageNormalize(vmin=2., vmax=10., stretch=None)

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
                              clip=False, invalid=None)
        norm2 = ImageNormalize(DATA, interval=ManualInterval(2, 10),
                               stretch=SqrtStretch(), clip=False,
                               invalid=None)
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

    def test_call_clip(self):
        """Test that the clip keyword is used when calling the object."""
        data = np.arange(5)
        norm = ImageNormalize(vmin=1., vmax=3., clip=False)

        output = norm(data, clip=True)
        assert_equal(output.data, [0, 0, 0.5, 1.0, 1.0])
        assert np.all(~output.mask)

        output = norm(data, clip=False)
        assert_equal(output.data, [-0.5, 0, 0.5, 1.0, 1.5])
        assert np.all(~output.mask)

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
                              clip=False, invalid=None)
        norm2 = ImageNormalize(mdata, interval=ManualInterval(2, 10),
                               stretch=SqrtStretch(), clip=False,
                               invalid=None)
        output = norm(mdata)
        expected = [np.nan, 0.35355339, -10, 0.93541435, 1.11803399,
                    1.27475488]
        assert_allclose(output.filled(-10), expected)
        assert_allclose(output.mask, [0, 0, 1, 0, 0, 0])

        assert_allclose(norm.inverse(norm(DATA))[1:], DATA[1:])
        assert_allclose(output, norm2(mdata))

    def test_invalid_data(self):
        data = np.arange(25.).reshape((5, 5))
        data[2, 2] = np.nan
        data[1, 2] = np.inf
        percent = 85.0
        interval = PercentileInterval(percent)

        # initialized without data
        norm = ImageNormalize(interval=interval)
        norm(data)  # sets vmin/vmax
        assert_equal((norm.vmin, norm.vmax), (1.65, 22.35))

        # initialized with data
        norm2 = ImageNormalize(data, interval=interval)
        assert_equal((norm2.vmin, norm2.vmax), (norm.vmin, norm.vmax))

        norm3 = simple_norm(data, 'linear', percent=percent)
        assert_equal((norm3.vmin, norm3.vmax), (norm.vmin, norm.vmax))

        assert_allclose(norm(data), norm2(data))
        assert_allclose(norm(data), norm3(data))

        norm4 = ImageNormalize()
        norm4(data)  # sets vmin/vmax
        assert_equal((norm4.vmin, norm4.vmax), (0, 24))

        norm5 = ImageNormalize(data)
        assert_equal((norm5.vmin, norm5.vmax), (norm4.vmin, norm4.vmax))

    @pytest.mark.parametrize('stretch', STRETCHES)
    def test_invalid_keyword(self, stretch):
        norm1 = ImageNormalize(stretch=stretch, vmin=-1, vmax=1, clip=False,
                               invalid=None)
        norm2 = ImageNormalize(stretch=stretch, vmin=-1, vmax=1, clip=False)
        norm3 = ImageNormalize(DATA3, stretch=stretch, vmin=-1, vmax=1,
                               clip=False, invalid=-1.)
        result1 = norm1(DATA3)
        result2 = norm2(DATA3)
        result3 = norm3(DATA3)
        assert_equal(result1[0:2], (np.nan, np.nan))
        assert_equal(result2[0:2], (-1., -1.))
        assert_equal(result1[2:], result2[2:])
        assert_equal(result2, result3)


@pytest.mark.skipif('not HAS_MATPLOTLIB')
class TestImageScaling:

    def test_linear(self):
        """Test linear scaling."""
        norm = simple_norm(DATA2, stretch='linear')
        assert_allclose(norm(DATA2), DATA2SCL, atol=0, rtol=1.e-5)

    def test_sqrt(self):
        """Test sqrt scaling."""
        norm1 = simple_norm(DATA2, stretch='sqrt')
        assert_allclose(norm1(DATA2), np.sqrt(DATA2SCL), atol=0, rtol=1.e-5)

    @pytest.mark.parametrize('invalid', INVALID)
    def test_sqrt_invalid_kw(self, invalid):
        stretch = SqrtStretch()
        norm1 = simple_norm(DATA3, stretch='sqrt', min_cut=-1, max_cut=1,
                            clip=False, invalid=invalid)
        norm2 = ImageNormalize(stretch=stretch, vmin=-1, vmax=1, clip=False,
                               invalid=invalid)
        assert_equal(norm1(DATA3), norm2(DATA3))

    def test_power(self):
        """Test power scaling."""
        power = 3.0
        norm = simple_norm(DATA2, stretch='power', power=power)
        assert_allclose(norm(DATA2), DATA2SCL ** power, atol=0, rtol=1.e-5)

    def test_log(self):
        """Test log10 scaling."""
        norm = simple_norm(DATA2, stretch='log')
        ref = np.log10(1000 * DATA2SCL + 1.0) / np.log10(1001.0)
        assert_allclose(norm(DATA2), ref, atol=0, rtol=1.e-5)

    def test_log_with_log_a(self):
        """Test log10 scaling with a custom log_a."""
        log_a = 100
        norm = simple_norm(DATA2, stretch='log', log_a=log_a)
        ref = np.log10(log_a * DATA2SCL + 1.0) / np.log10(log_a + 1)
        assert_allclose(norm(DATA2), ref, atol=0, rtol=1.e-5)

    def test_asinh(self):
        """Test asinh scaling."""
        norm = simple_norm(DATA2, stretch='asinh')
        ref = np.arcsinh(10 * DATA2SCL) / np.arcsinh(10)
        assert_allclose(norm(DATA2), ref, atol=0, rtol=1.e-5)

    def test_asinh_with_asinh_a(self):
        """Test asinh scaling with a custom asinh_a."""
        asinh_a = 0.5
        norm = simple_norm(DATA2, stretch='asinh', asinh_a=asinh_a)
        ref = np.arcsinh(DATA2SCL / asinh_a) / np.arcsinh(1. / asinh_a)
        assert_allclose(norm(DATA2), ref, atol=0, rtol=1.e-5)

    def test_min(self):
        """Test linear scaling."""
        norm = simple_norm(DATA2, stretch='linear', min_cut=1., clip=True)
        assert_allclose(norm(DATA2), [0., 0., 1.], atol=0, rtol=1.e-5)

    def test_percent(self):
        """Test percent keywords."""
        norm = simple_norm(DATA2, stretch='linear', percent=99., clip=True)
        assert_allclose(norm(DATA2), DATA2SCL, atol=0, rtol=1.e-5)

        norm2 = simple_norm(DATA2, stretch='linear', min_percent=0.5,
                            max_percent=99.5, clip=True)
        assert_allclose(norm(DATA2), norm2(DATA2), atol=0, rtol=1.e-5)

    def test_invalid_stretch(self):
        """Test invalid stretch keyword."""
        with pytest.raises(ValueError):
            simple_norm(DATA2, stretch='invalid')


@pytest.mark.skipif('not HAS_PLT')
def test_imshow_norm():
    import matplotlib.pyplot as plt
    image = np.random.randn(10, 10)

    plt.clf()
    ax = plt.subplot(label='test_imshow_norm')
    imshow_norm(image, ax=ax)

    with pytest.raises(ValueError):
        # X and data are the same, can't give both
        imshow_norm(image, X=image, ax=ax)

    with pytest.raises(ValueError):
        # illegal to manually pass in normalization since that defeats the point
        imshow_norm(image, ax=ax, norm=ImageNormalize())

    plt.clf()
    imshow_norm(image, ax=ax, vmin=0, vmax=1)

    # make sure the pyplot version works
    plt.clf()
    imres, norm = imshow_norm(image, ax=None)

    assert isinstance(norm, ImageNormalize)

    plt.close('all')
