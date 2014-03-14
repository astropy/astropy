import numpy as np

from ..angles import Longitude, Latitude, Angle
from ...tests.helper import pytest
from ... import units as u

def test_colon():
    a = Angle(1.113355, unit=u.deg)
    assert a.to_string(fields=2, sep=':') == u'1:07'
    assert a.to_string(fields=2, sep='dms') == u'1d07m'
    assert a.to_string(fields=3, sep=':') == u'1:06:48.078'
    assert a.to_string(fields=1, sep=':') == u'1'
