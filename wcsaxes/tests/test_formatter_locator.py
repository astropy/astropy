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

    def test_format(self):

        fl = AngleFormatterLocator(number=5)

        fl.format = 'dd'
        assert fl.formatter([3.], None)[0] == u'3d'

        fl.format = 'dd:mm'
        assert fl.formatter([3.1], None)[0] == u'3d06m'

        fl.format = 'dd:mm:ss'
        assert fl.formatter([3.12], None)[0] == u'3d07m12s'

        fl.format = 'dd:mm:ss.s'
        assert fl.formatter([3.12], None)[0] == u'3d07m12.0s'

        fl.format = 'dd:mm:ss.ssss'
        assert fl.formatter([3.12], None)[0] == u'3d07m12.0000s'

        fl.format = 'hh'
        assert fl.formatter([15.], None)[0] == u'1h'

        fl.format = 'hh:mm'
        assert fl.formatter([15.5], None)[0] == u'1h02m'

        fl.format = 'hh:mm:ss'
        assert fl.formatter([15.55], None)[0] == u'1h02m12s'

        fl.format = 'hh:mm:ss.s'
        assert fl.formatter([15.55], None)[0] == u'1h02m12.0s'

        fl.format = 'hh:mm:ss.ssss'
        assert fl.formatter([15.55], None)[0] == u'1h02m12.0000s'

        fl.format = 'd'
        assert fl.formatter([15.392231], None)[0] == u'15'

        fl.format = 'd.d'
        assert fl.formatter([15.392231], None)[0] == u'15.4'

        fl.format = 'd.dd'
        assert fl.formatter([15.392231], None)[0] == u'15.39'

        fl.format = 'd.ddd'
        assert fl.formatter([15.392231], None)[0] == u'15.392'

    @pytest.mark.parametrize(('format'), ['x.xxx', 'dd.ss', 'dd:ss', 'mdd:mm:ss'])
    def test_invalid_formats(self, format):
        fl = AngleFormatterLocator(number=5)
        with pytest.raises(ValueError) as exc:
            fl.format = format
        assert exc.value.args[0] == "Invalid format: " + format

    def test_base_spacing(self):

        fl = AngleFormatterLocator(number=5)

        fl.format = 'dd'
        assert fl.base_spacing == 1 * u.deg

        fl.format = 'dd:mm'
        assert fl.base_spacing == 1 * u.arcmin

        fl.format = 'dd:mm:ss'
        assert fl.base_spacing == 1 * u.arcsec

        fl.format = 'dd:mm:ss.ss'
        assert fl.base_spacing == 0.01 * u.arcsec

        fl.format = 'hh'
        assert fl.base_spacing == 15. * u.deg

        fl.format = 'hh:mm'
        assert fl.base_spacing == 15. * u.arcmin

        fl.format = 'hh:mm:ss'
        assert fl.base_spacing == 15. * u.arcsec

        fl.format = 'hh:mm:ss.ss'
        assert fl.base_spacing == 0.15 * u.arcsec
