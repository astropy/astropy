from ..angles import Angle
from ... import units as u

def test_format_precision():

    # There are already some tests in test_api.py, but this is a regression
    # test for issue #... which caused incorrect formatting of the seconds for
    # precision=0

    angle = Angle(-1.23456789, unit=u.degree)

    assert angle.format(precision=3) == '-1d14m04.444s'
    assert angle.format(precision=1) == '-1d14m04.4s'
    assert angle.format(precision=0) == '-1d14m04s'

    angle2 = Angle(-1.23456789, unit=u.hour)

    assert angle2.format(precision=3, unit=u.hour) == '-1h14m04.444s'
    assert angle2.format(precision=1, unit=u.hour) == '-1h14m04.4s'
    assert angle2.format(precision=0, unit=u.hour) == '-1h14m04s'
