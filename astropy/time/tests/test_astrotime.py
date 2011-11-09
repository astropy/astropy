# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..astrotime import AstroTime
import datetime

def test_jd2gregorian():
    test_date = AstroTime.from_jd(2449443.5)
    assert test_date.datetime == datetime.datetime(1994, 4, 1, 0, 0, 0)

def test_gregorian2jd():
    test_date = AstroTime(datetime.datetime(1994, 4, 1, 0, 0, 0))
    #Should be assert float
    assert test_date.jd == 2449443.5
    
def test_deltaastrotime():
    a = AstroTime(10)
    b = AstroTime(-5)
    
    d = a - b
    assert d._val == 15000000000L
    assert (d+b)._val == 10000000000L