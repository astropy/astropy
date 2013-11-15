from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six
from ...tests.helper import pytest

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
    assert a.to_string(format='latex') == r'$1^\circ06{}^\prime48.078{}^{\prime\prime}$'
    assert a.to_string(format='unicode') == '1\xb006\u203248.078\u2033'

    a = Angle(1.113355, unit=u.hour)
    assert a.to_string(format='latex') == r'$1^\mathrm{h}06^\mathrm{m}48.078^\mathrm{s}$'
    assert a.to_string(format='unicode') == '1\u02b006\u1d5048.078\u02e2'

    a = Angle(1.113355, unit=u.radian)
    assert a.to_string(format='latex') == r'$1.11336\mathrm{rad}$'
    assert a.to_string(format='unicode') == '1.11336rad'


def test_to_string_fields():
    a = Angle(1.113355, unit=u.deg)
    assert a.to_string(fields=1) == r'1d'
    assert a.to_string(fields=2) == r'1d07m'
    assert a.to_string(fields=3) == r'1d06m48.078s'


def test_sexagesimal_rounding_up():
    a = Angle(359.9999999999, unit=u.deg)

    assert a.to_string(precision=None) == '360d00m00s'
    assert a.to_string(precision=4) == '360d00m00.0000s'
    assert a.to_string(precision=5) == '360d00m00.00000s'
    assert a.to_string(precision=6) == '360d00m00.000000s'
    assert a.to_string(precision=7) == '359d59m59.9999996s'

    a = Angle(3.999999, unit=u.deg)
    assert a.to_string(fields=2, precision=None) == '4d00m'
    assert a.to_string(fields=2, precision=1) == '4d00m'
    assert a.to_string(fields=2, precision=5) == '4d00m'
    assert a.to_string(fields=1, precision=1) == '4d'
    assert a.to_string(fields=1, precision=5) == '4d'


def test_to_string_scalar():
    a = Angle(1.113355, unit=u.deg)
    assert isinstance(a.to_string(), six.text_type)


import numpy as np
from .. import ICRS, FK4, FK4NoETerms, FK5, Galactic, AltAz


@pytest.mark.parametrize('frame', [ICRS, FK4, FK4NoETerms, FK5])
def test_coordinate_to_string_vector_hms(frame):

    C = frame(np.arange(2)*12.05*u.deg, np.arange(2)*13.5*u.deg)
    assert C.to_string(precision=0) == ['-0h00m00s 0d00m00s', '0h48m12s 13d30m00s']
    assert C.to_string(precision=1) == ['-0h00m00.0s 0d00m00.0s', '0h48m12.0s 13d30m00.0s']


@pytest.mark.parametrize('frame', [Galactic, AltAz])
def test_coordinate_to_string_vector_dms(frame):

    C = frame(np.arange(2)*12.05*u.deg, np.arange(2)*13.5*u.deg)
    assert C.to_string(precision=0) == ['-0d00m00s 0d00m00s', '12d03m00s 13d30m00s']
    assert C.to_string(precision=1) == ['-0d00m00.0s 0d00m00.0s', '12d03m00.0s 13d30m00.0s']


@pytest.mark.parametrize('frame', [ICRS, FK4, FK4NoETerms, FK5])
def test_coordinate_to_string_scalar_hms(frame):

    C = frame(12.05*u.deg, 13.5*u.deg)
    assert C.to_string(precision=0) == '0h48m12s 13d30m00s'
    assert C.to_string(precision=1) == '0h48m12.0s 13d30m00.0s'


@pytest.mark.parametrize('frame', [Galactic, AltAz])
def test_coordinate_to_string_scalar_dms(frame):

    C = frame(12.05*u.deg, 13.5*u.deg)
    assert C.to_string(precision=0) == '12d03m00s 13d30m00s'
    assert C.to_string(precision=1) == '12d03m00.0s 13d30m00.0s'


def test_to_string_radian_with_precision():
    """
    Regression test for a bug that caused ``to_string`` to crash for angles in
    radians when specifying the precision.
    """

    # Check that specifying the precision works
    a = Angle(3., unit=u.rad)
    assert a.to_string(precision=3, sep='fromunit') == '3.000rad'
