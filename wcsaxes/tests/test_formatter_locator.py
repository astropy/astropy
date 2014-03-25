import pytest

import numpy as np
from numpy.testing import assert_almost_equal
from astropy import units as u
from .. import six

from ..formatter_locator import AngleFormatterLocator, ScalarFormatterLocator


class TestAngleFormatterLocator(object):

    def test_no_options(self):

        fl = AngleFormatterLocator()
        assert fl.values is None
        assert fl.number == 5
        assert fl.spacing is None

    def test_too_many_options(self):

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(values=[1.,2.], number=5)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(values=[1.,2.], spacing=5. * u.deg)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(number=5, spacing=5. * u.deg)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(values=[1.,2.], number=5, spacing=5. * u.deg)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

    def test_values(self):

        fl = AngleFormatterLocator(values=[0.1, 1., 14.])
        assert fl.values == [0.1, 1., 14.]
        assert fl.number is None
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values, [0.1, 1., 14.])

    def test_number(self):

        fl = AngleFormatterLocator(number=7)
        assert fl.values is None
        assert fl.number == 7
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values, [35., 40., 45., 50., 55.])

        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, [34.5, 34.75, 35., 35.25, 35.5, 35.75, 36.])

        fl.format = 'dd'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, [35., 36.])

    def test_spacing(self):

        with pytest.raises(TypeError) as exc:
            AngleFormatterLocator(spacing=3.)
        assert exc.value.args[0] == "spacing should be an astropy.units.Quantity instance with units of angle"

        fl = AngleFormatterLocator(spacing=3. * u.degree)
        assert fl.values is None
        assert fl.number is None
        assert fl.spacing == 3. * u.degree

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values, [36., 39., 42., 45., 48., 51., 54.])

        fl.spacing = 30. * u.arcmin
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, [34.5, 35., 35.5, 36.])

        fl.format = 'dd'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, [35., 36.])

    @pytest.mark.parametrize(('format', 'string'), [('dd', six.u('15\xb0')),
                                                    ('dd:mm', six.u('15\xb024\'')),
                                                    ('dd:mm:ss', six.u('15\xb023\'32"')),
                                                    ('dd:mm:ss.s', six.u('15\xb023\'32.0"')),
                                                    ('dd:mm:ss.ssss', six.u('15\xb023\'32.0316"')),
                                                    ('hh', '1h'),
                                                    ('hh:mm', '1h02m'),
                                                    ('hh:mm:ss', '1h01m34s'),
                                                    ('hh:mm:ss.s', '1h01m34.1s'),
                                                    ('hh:mm:ss.ssss', '1h01m34.1354s'),
                                                    ('d', '15'),
                                                    ('d.d', '15.4'),
                                                    ('d.dd', '15.39'),
                                                    ('d.ddd', '15.392'),
                                                    ('m', '924'),
                                                    ('m.m', '923.5'),
                                                    ('m.mm', '923.53'),
                                                    ('s', '55412'),
                                                    ('s.s', '55412.0'),
                                                    ('s.ss', '55412.03'),
                                                    ])
    def test_format(self, format, string):
        fl = AngleFormatterLocator(number=5, format=format)
        assert fl.formatter([15.392231], None)[0] == string

    @pytest.mark.parametrize(('format'), ['x.xxx', 'dd.ss', 'dd:ss', 'mdd:mm:ss'])
    def test_invalid_formats(self, format):
        fl = AngleFormatterLocator(number=5)
        with pytest.raises(ValueError) as exc:
            fl.format = format
        assert exc.value.args[0] == "Invalid format: " + format

    @pytest.mark.parametrize(('format', 'base_spacing'), [('dd', 1. * u.deg),
                                                          ('dd:mm', 1. * u.arcmin),
                                                          ('dd:mm:ss', 1. * u.arcsec),
                                                          ('dd:mm:ss.ss', 0.01 * u.arcsec),
                                                          ('hh', 15. * u.deg),
                                                          ('hh:mm', 15. * u.arcmin),
                                                          ('hh:mm:ss', 15. * u.arcsec),
                                                          ('hh:mm:ss.ss', 0.15 * u.arcsec),
                                                          ('d', 1. * u.deg),
                                                          ('d.d', 0.1 * u.deg),
                                                          ('d.dd', 0.01 * u.deg),
                                                          ('d.ddd', 0.001 * u.deg),
                                                          ('m', 1. * u.arcmin),
                                                          ('m.m', 0.1 * u.arcmin),
                                                          ('m.mm', 0.01 * u.arcmin),
                                                          ('s', 1. * u.arcsec),
                                                          ('s.s', 0.1 * u.arcsec),
                                                          ('s.ss', 0.01 * u.arcsec),
                                                          ])
    def test_base_spacing(self, format, base_spacing):
        fl = AngleFormatterLocator(number=5, format=format)
        assert fl.base_spacing == base_spacing


class TestScalarFormatterLocator(object):

    def test_no_options(self):

        fl = ScalarFormatterLocator()
        assert fl.values is None
        assert fl.number == 5
        assert fl.spacing is None

    def test_too_many_options(self):

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(values=[1.,2.], number=5)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(values=[1.,2.], spacing=5.)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(number=5, spacing=5.)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(values=[1.,2.], number=5, spacing=5.)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

    def test_values(self):

        fl = ScalarFormatterLocator(values=[0.1, 1., 14.])
        assert fl.values == [0.1, 1., 14.]
        assert fl.number is None
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values, [0.1, 1., 14.])

    def test_number(self):

        fl = ScalarFormatterLocator(number=7)
        assert fl.values is None
        assert fl.number == 7
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values, np.linspace(36., 54., 10))

        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, np.linspace(34.4, 36, 9))

        fl.format = 'x'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, [35., 36.])

    def test_spacing(self):

        fl = ScalarFormatterLocator(spacing=3.)
        assert fl.values is None
        assert fl.number is None
        assert fl.spacing == 3.

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values, [36., 39., 42., 45., 48., 51., 54.])

        fl.spacing = 0.5
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, [34.5, 35., 35.5, 36.])

        fl.format = 'x'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values, [35., 36.])

    @pytest.mark.parametrize(('format', 'string'), [('x', '15'),
                                                    ('x.x', '15.4'),
                                                    ('x.xx', '15.39'),
                                                    ('x.xxx', '15.392')])
    def test_format(self, format, string):
        fl = ScalarFormatterLocator(number=5, format=format)
        assert fl.formatter([15.392231], None)[0] == string

    @pytest.mark.parametrize(('format'), ['dd', 'dd:mm', 'xx:mm', 'mx.xxx'])
    def test_invalid_formats(self, format):
        fl = ScalarFormatterLocator(number=5)
        with pytest.raises(ValueError) as exc:
            fl.format = format
        assert exc.value.args[0] == "Invalid format: " + format

    @pytest.mark.parametrize(('format', 'base_spacing'), [('x', 1.),
                                                          ('x.x', 0.1),
                                                          ('x.xxx', 0.001)])
    def test_base_spacing(self, format, base_spacing):
        fl = ScalarFormatterLocator(number=5, format=format)
        assert fl.base_spacing == base_spacing
