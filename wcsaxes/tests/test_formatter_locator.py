import pytest

from numpy.testing import assert_almost_equal
from astropy import units as u

from ..formatter_locator import AngleFormatterLocator


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

    expected_formatted = {}
    expected_formatted['dd'] = u'15d'
    expected_formatted['dd:mm'] = u'15d24m'
    expected_formatted['dd:mm:ss'] = u'15d23m32s'
    expected_formatted['dd:mm:ss.s'] = u'15d23m32.0s'
    expected_formatted['dd:mm:ss.ssss'] = u'15d23m32.0316s'
    expected_formatted['hh'] = u'1h'
    expected_formatted['hh:mm'] = u'1h02m'
    expected_formatted['hh:mm:ss'] = u'1h01m34s'
    expected_formatted['hh:mm:ss.s'] = u'1h01m34.1s'
    expected_formatted['hh:mm:ss.ssss'] = u'1h01m34.1354s'
    expected_formatted['d'] = u'15'
    expected_formatted['d.d'] = u'15.4'
    expected_formatted['d.dd'] = u'15.39'
    expected_formatted['d.ddd'] = u'15.392'

    @pytest.mark.parametrize(('format', 'string'),
                             [(x, expected_formatted[x]) for x in expected_formatted])
    def test_format(self, format, string):
        fl = AngleFormatterLocator(number=5, format=format)
        assert fl.formatter([15.392231], None)[0] == string

    @pytest.mark.parametrize(('format'), ['x.xxx', 'dd.ss', 'dd:ss', 'mdd:mm:ss'])
    def test_invalid_formats(self, format):
        fl = AngleFormatterLocator(number=5)
        with pytest.raises(ValueError) as exc:
            fl.format = format
        assert exc.value.args[0] == "Invalid format: " + format

    expected_spacing = {}
    expected_spacing['dd'] = 1. * u.deg
    expected_spacing['dd:mm'] = 1. * u.arcmin
    expected_spacing['dd:mm:ss'] = 1. * u.arcsec
    expected_spacing['dd:mm:ss.ss'] = 0.01 * u.arcsec
    expected_spacing['hh'] = 15. * u.deg
    expected_spacing['hh:mm'] = 15. * u.arcmin
    expected_spacing['hh:mm:ss'] = 15. * u.arcsec
    expected_spacing['hh:mm:ss.ss'] = 0.15 * u.arcsec

    @pytest.mark.parametrize(('format', 'base_spacing'),
                             [(x, expected_spacing[x]) for x in expected_spacing])
    def test_base_spacing(self, format, base_spacing):
        fl = AngleFormatterLocator(number=5, format=format)
        assert fl.base_spacing == base_spacing
