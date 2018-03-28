# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_almost_equal

from matplotlib import rc_context

from ....tests.helper import assert_quantity_allclose
from .... import units as u
from ..formatter_locator import AngleFormatterLocator, ScalarFormatterLocator


class TestAngleFormatterLocator:

    def test_no_options(self):

        fl = AngleFormatterLocator()
        assert fl.values is None
        assert fl.number == 5
        assert fl.spacing is None

    def test_too_many_options(self):

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(values=[1., 2.], number=5)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(values=[1., 2.], spacing=5. * u.deg)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(number=5, spacing=5. * u.deg)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            AngleFormatterLocator(values=[1., 2.], number=5, spacing=5. * u.deg)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

    def test_values(self):

        fl = AngleFormatterLocator(values=[0.1, 1., 14.] * u.degree)
        assert fl.values.to_value(u.degree).tolist() == [0.1, 1., 14.]
        assert fl.number is None
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values.to_value(u.degree), [0.1, 1., 14.])

    def test_number(self):

        fl = AngleFormatterLocator(number=7)
        assert fl.values is None
        assert fl.number == 7
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values.to_value(u.degree), [35., 40., 45., 50., 55.])

        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.to_value(u.degree), [34.5, 34.75, 35., 35.25, 35.5, 35.75, 36.])

        fl.format = 'dd'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.to_value(u.degree), [35., 36.])

    def test_spacing(self):

        with pytest.raises(TypeError) as exc:
            AngleFormatterLocator(spacing=3.)
        assert exc.value.args[0] == "spacing should be an astropy.units.Quantity instance with units of angle"

        fl = AngleFormatterLocator(spacing=3. * u.degree)
        assert fl.values is None
        assert fl.number is None
        assert fl.spacing == 3. * u.degree

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values.to_value(u.degree), [36., 39., 42., 45., 48., 51., 54.])

        fl.spacing = 30. * u.arcmin
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.to_value(u.degree), [34.5, 35., 35.5, 36.])

        fl.format = 'dd'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.to_value(u.degree), [35., 36.])

    def test_minor_locator(self):

        fl = AngleFormatterLocator()

        values, spacing = fl.locator(34.3, 55.4)

        minor_values = fl.minor_locator(spacing, 5, 34.3, 55.4)

        assert_almost_equal(minor_values.to_value(u.degree), [36., 37., 38.,
                            39., 41., 42., 43., 44., 46., 47., 48., 49., 51.,
                            52., 53., 54.])

        minor_values = fl.minor_locator(spacing, 2, 34.3, 55.4)

        assert_almost_equal(minor_values.to_value(u.degree), [37.5, 42.5, 47.5, 52.5])

        fl.values = [0.1, 1., 14.] * u.degree

        values, spacing = fl.locator(34.3, 36.1)

        minor_values = fl.minor_locator(spacing, 2, 34.3, 55.4)

        assert_almost_equal(minor_values.to_value(u.degree), [])

    @pytest.mark.parametrize(('format', 'string'), [('dd', '15\xb0'),
                                                    ('dd:mm', '15\xb024\''),
                                                    ('dd:mm:ss', '15\xb023\'32"'),
                                                    ('dd:mm:ss.s', '15\xb023\'32.0"'),
                                                    ('dd:mm:ss.ssss', '15\xb023\'32.0316"'),
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
        assert fl.formatter([15.392231] * u.degree, None)[0] == string

    @pytest.mark.parametrize(('separator', 'format', 'string'), [(('deg', "'", '"'), 'dd', '15deg'),
                                                                 (('deg', "'", '"'), 'dd:mm', '15deg24\''),
                                                                 (('deg', "'", '"'), 'dd:mm:ss', '15deg23\'32"'),
                                                                 ((':', "-", 's'), 'dd:mm:ss.s', '15:23-32.0s'),
                                                                 (':', 'dd:mm:ss.s', '15:23:32.0'),
                                                                 ((':', ":", 's'), 'hh', '1:'),
                                                                 (('-', "-", 's'), 'hh:mm:ss.ssss', '1-01-34.1354s'),
                                                                 (('d', ":", '"'), 'd', '15'),
                                                                 (('d', ":", '"'), 'd.d', '15.4'),
                                                                 ])
    def test_separator(self, separator, format, string):
        fl = AngleFormatterLocator(number=5, format=format)
        fl.sep = separator
        assert fl.formatter([15.392231] * u.degree, None)[0] == string

    def test_latex_format(self):
        fl = AngleFormatterLocator(number=5, format="dd:mm:ss")
        assert fl.formatter([15.392231] * u.degree, None)[0] == '15\xb023\'32"'
        with rc_context(rc={'text.usetex': True}):
            assert fl.formatter([15.392231] * u.degree, None)[0] == "15$^\\circ$23'32\""

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

    def test_incorrect_spacing(self):
        fl = AngleFormatterLocator()
        fl.spacing = 0.032 * u.deg
        fl.format = 'dd:mm:ss'
        assert_almost_equal(fl.spacing.to_value(u.arcsec), 115.)

    def test_decimal_values(self):

        # Regression test for a bug that meant that the spacing was not
        # determined correctly for decimal coordinates

        fl = AngleFormatterLocator()
        fl.format = 'd.dddd'
        assert_quantity_allclose(fl.locator(266.9730, 266.9750)[0],
                                 [266.9735, 266.9740, 266.9745, 266.9750] * u.deg)

    @pytest.mark.parametrize(('spacing', 'string'), [(2 * u.deg, '15\xb0'),
                                                     (2 * u.arcmin, '15\xb024\''),
                                                     (2 * u.arcsec, '15\xb023\'32"'),
                                                     (0.1 * u.arcsec, '15\xb023\'32.0"')])
    def test_formatter_no_format(self, spacing, string):
        fl = AngleFormatterLocator()
        assert fl.formatter([15.392231] * u.degree, spacing)[0] == string


class TestScalarFormatterLocator:

    def test_no_options(self):

        fl = ScalarFormatterLocator(unit=u.m)
        assert fl.values is None
        assert fl.number == 5
        assert fl.spacing is None

    def test_too_many_options(self):

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(values=[1., 2.] * u.m, number=5)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(values=[1., 2.] * u.m, spacing=5. * u.m)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(number=5, spacing=5. * u.m)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

        with pytest.raises(ValueError) as exc:
            ScalarFormatterLocator(values=[1., 2.] * u.m, number=5, spacing=5. * u.m)
        assert exc.value.args[0] == "At most one of values/number/spacing can be specifed"

    def test_values(self):

        fl = ScalarFormatterLocator(values=[0.1, 1., 14.] * u.m, unit=u.m)
        assert fl.values.value.tolist() == [0.1, 1., 14.]
        assert fl.number is None
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values.value, [0.1, 1., 14.])

    def test_number(self):

        fl = ScalarFormatterLocator(number=7, unit=u.m)
        assert fl.values is None
        assert fl.number == 7
        assert fl.spacing is None

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values.value, np.linspace(36., 54., 10))

        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.value, np.linspace(34.4, 36, 9))

        fl.format = 'x'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.value, [35., 36.])

    def test_spacing(self):

        fl = ScalarFormatterLocator(spacing=3. * u.m)
        assert fl.values is None
        assert fl.number is None
        assert fl.spacing == 3. * u.m

        values, spacing = fl.locator(34.3, 55.4)
        assert_almost_equal(values.value, [36., 39., 42., 45., 48., 51., 54.])

        fl.spacing = 0.5 * u.m
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.value, [34.5, 35., 35.5, 36.])

        fl.format = 'x'
        values, spacing = fl.locator(34.3, 36.1)
        assert_almost_equal(values.value, [35., 36.])

    def test_minor_locator(self):

        fl = ScalarFormatterLocator(unit=u.m)

        values, spacing = fl.locator(34.3, 55.4)

        minor_values = fl.minor_locator(spacing, 5, 34.3, 55.4)

        assert_almost_equal(minor_values.value, [36., 37., 38., 39., 41., 42.,
                            43., 44., 46., 47., 48., 49., 51., 52., 53., 54.])
        print('minor_values: ' + str(minor_values))

        minor_values = fl.minor_locator(spacing, 2, 34.3, 55.4)

        assert_almost_equal(minor_values.value, [37.5, 42.5, 47.5, 52.5])

        fl.values = [0.1, 1., 14.] * u.m

        values, spacing = fl.locator(34.3, 36.1)

        minor_values = fl.minor_locator(spacing, 2, 34.3, 55.4)

        assert_almost_equal(minor_values.value, [])

    @pytest.mark.parametrize(('format', 'string'), [('x', '15'),
                                                    ('x.x', '15.4'),
                                                    ('x.xx', '15.39'),
                                                    ('x.xxx', '15.392'),
                                                    ('%g', '15.3922'),
                                                    ('%f', '15.392231'),
                                                    ('%.2f', '15.39'),
                                                    ('%.3f', '15.392')])
    def test_format(self, format, string):
        fl = ScalarFormatterLocator(number=5, format=format, unit=u.m)
        assert fl.formatter([15.392231] * u.m, None)[0] == string

    @pytest.mark.parametrize(('format', 'string'), [('x', '1539'),
                                                    ('x.x', '1539.2'),
                                                    ('x.xx', '1539.22'),
                                                    ('x.xxx', '1539.223')])
    def test_format_unit(self, format, string):
        fl = ScalarFormatterLocator(number=5, format=format, unit=u.m)
        fl.format_unit = u.cm
        assert fl.formatter([15.392231] * u.m, None)[0] == string

    @pytest.mark.parametrize(('format'), ['dd', 'dd:mm', 'xx:mm', 'mx.xxx'])
    def test_invalid_formats(self, format):
        fl = ScalarFormatterLocator(number=5, unit=u.m)
        with pytest.raises(ValueError) as exc:
            fl.format = format
        assert exc.value.args[0] == "Invalid format: " + format

    @pytest.mark.parametrize(('format', 'base_spacing'), [('x', 1. * u.m),
                                                          ('x.x', 0.1 * u.m),
                                                          ('x.xxx', 0.001 * u.m)])
    def test_base_spacing(self, format, base_spacing):
        fl = ScalarFormatterLocator(number=5, format=format, unit=u.m)
        assert fl.base_spacing == base_spacing

    def test_incorrect_spacing(self):
        fl = ScalarFormatterLocator(unit=u.m)
        fl.spacing = 0.032 * u.m
        fl.format = 'x.xx'
        assert_almost_equal(fl.spacing.to_value(u.m), 0.03)
