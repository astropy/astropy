from ..angles import Angle
from ... import units as u


def test_format_precision():

    # There are already some tests in test_api.py, but this is a regression
    # test for the bug in issue #1319 which caused incorrect formatting of the
    # seconds for precision=0

    angle = Angle(-1.23456789, unit=u.degree)

    assert angle.format(precision=3) == '-1d14m04.444s'
    assert angle.format(precision=1) == '-1d14m04.4s'
    assert angle.format(precision=0) == '-1d14m04s'

    angle2 = Angle(-1.23456789, unit=u.hourangle)

    assert angle2.format(precision=3, unit=u.hour) == '-1h14m04.444s'
    assert angle2.format(precision=1, unit=u.hour) == '-1h14m04.4s'
    assert angle2.format(precision=0, unit=u.hour) == '-1h14m04s'


def test_format_decimal():

    # There are already some tests in test_api.py, but this is a regression
    # test for the bug in issue #1323 which caused decimal formatting to not
    # work

    angle1 = Angle(2., unit=u.degree)

    assert angle1.format(decimal=True, precision=3) == '2.000'
    assert angle1.format(decimal=True, precision=1) == '2.0'
    assert angle1.format(decimal=True, precision=0) == '2'

    angle2 = Angle(3., unit=u.hourangle)

    # In 0.3 and above .format() should preserve the Angle's input unit, but
    # the change that affects this is not backported to 0.2.x since it depends
    # on deeper changes
    assert angle2.format(decimal=True, precision=3, unit=u.hour) == '3.000'
    assert angle2.format(decimal=True, precision=1, unit=u.hour) == '3.0'
    assert angle2.format(decimal=True, precision=0, unit=u.hour) == '3'

    angle3 = Angle(4., unit=u.radian)

    assert angle3.format(decimal=True, precision=3, unit=u.radian) == '4.000'
    assert angle3.format(decimal=True, precision=1, unit=u.radian) == '4.0'
    assert angle3.format(decimal=True, precision=0, unit=u.radian) == '4'
