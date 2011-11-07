from .. import astrotime

def test_jd2gregorian():
    test_date = astrotime.AstroTime.from_jd(2449443.5)
    assert test_date.to_gregorian_date() == astrotime.datetime(1994, 4, 1, 0, 0, 0)

def test_gregorian2jd():
    test_date = astrotime.AstroTime.from_gregorian_date(astrotime.datetime(1994, 4, 1, 0, 0, 0))
    #Should be assert float
    assert test_date.to_jd() == 2449443.5