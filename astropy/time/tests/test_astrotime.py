# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.time import Time
import datetime
import decimal
def test_jd2utc():
    test_date = Time.from_jd(2449443.5)
    assert test_date.utc == datetime.datetime(1994, 4, 1, 0, 0, 0)

def test_utc2jd():
    test_date = Time.from_utc(datetime.datetime(1994, 4, 1, 0, 0, 0))
    #less than 1 ms

    assert float(test_date.to_jd() - decimal.Decimal(2449443.5)) < 1e-9
