from .. import skyc
from .. import astrom


def test_skyc():
    assert skyc.dow(2451545.0) == 5.0
