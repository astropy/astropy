# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.time import Time
import datetime

def test_jd2gregorian():
    test_date = Time.from_jd(2449443.5)
    assert test_date.to_calendar_date() == datetime.datetime(1994, 4, 1, 0, 0, 0)

def test_gregorian2jd():
    test_date = Time.from_calendar_date(datetime.datetime(1994, 4, 1, 0, 0, 0))
    #Should be assert float
    assert test_date.to_jd() == 2449443.5