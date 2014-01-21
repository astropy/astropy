from ...extern.six.moves import cPickle as pickle
from ...coordinates import Longitude

def test_basic():
    lon1 = Longitude(1.23, "radian", wrap_angle='180d')
    s = pickle.dumps(lon1)
    lon2 = pickle.loads(s)

def test_pickle_longitude_wrap_angle():
    a = Longitude(1.23, "radian", wrap_angle='180d')
    s = pickle.dumps(a)
    b = pickle.loads(s)

    assert a.rad == b.rad
    assert a.wrap_angle == b.wrap_angle
