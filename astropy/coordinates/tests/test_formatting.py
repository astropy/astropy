from ..angles import Angle
from ... import units as u


def test_to_string_precision():
    # There are already some tests in test_api.py, but this is a regression
    # test for the bug in issue #1319 which caused incorrect formatting of the
    # seconds for precision=0

    angle = Angle(-1.23456789, unit=u.degree)

    assert angle.to_string(precision=3) == '-1d14m04.444s'
    assert angle.to_string(precision=1) == '-1d14m04.4s'
    assert angle.to_string(precision=0) == '-1d14m04s'

    angle2 = Angle(-1.23456789, unit=u.hourangle)

    assert angle2.to_string(precision=3, unit=u.hour) == '-1h14m04.444s'
    assert angle2.to_string(precision=1, unit=u.hour) == '-1h14m04.4s'
    assert angle2.to_string(precision=0, unit=u.hour) == '-1h14m04s'


def test_to_string_decimal():

    # There are already some tests in test_api.py, but this is a regression
    # test for the bug in issue #1323 which caused decimal formatting to not
    # work

    angle1 = Angle(2., unit=u.degree)

    assert angle1.to_string(decimal=True, precision=3) == '2.000'
    assert angle1.to_string(decimal=True, precision=1) == '2.0'
    assert angle1.to_string(decimal=True, precision=0) == '2'

    angle2 = Angle(3., unit=u.hourangle)

    assert angle2.to_string(decimal=True, precision=3) == '3.000'
    assert angle2.to_string(decimal=True, precision=1) == '3.0'
    assert angle2.to_string(decimal=True, precision=0) == '3'

    angle3 = Angle(4., unit=u.radian)

    assert angle3.to_string(decimal=True, precision=3) == '4.000'
    assert angle3.to_string(decimal=True, precision=1) == '4.0'
    assert angle3.to_string(decimal=True, precision=0) == '4'


def test_to_string_formats():
    a = Angle(1.113355, unit=u.deg)
    assert a.to_string(format='latex') == r'$1^\circ06{}^\prime48.07800{}^{\prime\prime}$'
    assert a.to_string(format='unicode') == u'1\xb006\u203248.07800\u2033'

    a = Angle(1.113355, unit=u.hour)
    assert a.to_string(format='latex') == r'$1^\mathrm{h}06^\mathrm{m}48.07800^\mathrm{s}$'
    assert a.to_string(format='unicode') == u'1\u02b006\u1d5048.07800\u02e2'

    a = Angle(1.113355, unit=u.radian)
    assert a.to_string(format='latex') == r'$1.11336\mathrm{rad}$'
    assert a.to_string(format='unicode') == '1.11336rad'


def test_to_string_fields():
    a = Angle(1.113355, unit=u.deg)
    assert a.to_string(fields=1) == r'1d'
    assert a.to_string(fields=2) == r'1d06m'
    assert a.to_string(fields=3) == r'1d06m48.07800s'
