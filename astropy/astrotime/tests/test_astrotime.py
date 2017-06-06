# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import astrotime
import datetime

def test_jd2gregorian():
    test_date = astrotime.AstroTime.from_jd(2449443.5)
    assert test_date.to_date_gregorian() == datetime.datetime(1994, 4, 1, 0, 0, 0)

def test_gregorian2jd():
    test_date = astrotime.AstroTime.from_date_gregorian(datetime.datetime(1994, 4, 1, 0, 0, 0))
    #Should be assert float
    assert test_date.to_jd() == 2449443.5