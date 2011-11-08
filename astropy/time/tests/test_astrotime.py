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
    
def test_addsub_astrotime():
    a = AstroTime(-5)
    b = AstroTime(10)
    
    c = a + b
    assert c._val == 5
    
    c = a - b
    assert c._val == -15
    
    a = AstroTime(1,1e-8)
    b = AstroTime(1,1e-9)
    assert (a+b)._val == 11
    
    
    a = AstroTime(0,1,0)
    b = AstroTime(0,1,1) #jd 1 later
    #uses the jd0 of the first one
    assert (a+b)._val == 86400
    assert (b+a)._val == -86400