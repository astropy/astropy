# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .time import AstroTime
import datetime

def test_jd2gregorian():
    test_date = AstroTime.from_jd(2449443.5)
    assert test_date.datetime == datetime.datetime(1994, 4, 1, 0, 0, 0)

def test_gregorian2jd():
    test_date = AstroTime(datetime.datetime(1994, 4, 1, 0, 0, 0))
    #Should be assert float
    assert test_date.jd == 2449443.5
    