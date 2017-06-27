import pytest
import numpy as np

from ...extern.six.moves import zip, cPickle as pickle
from ...coordinates import Longitude
from ... import coordinates as coord
from ...tests.helper import pickle_protocol, check_pickling_recovery  # noqa

# Can't test distances without scipy due to cosmology deps
try:
    import scipy  # pylint: disable=W0611
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


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


_names = [coord.Angle,
          coord.Distance,
          coord.DynamicMatrixTransform,
          coord.ICRS,
          coord.Latitude,
          coord.Longitude,
          coord.StaticMatrixTransform,
          ]

_xfail = [False,
          not HAS_SCIPY,
          True,
          True,
          False,
          True,
          False]

_args = [[0.0],
         [],
         [lambda *args: np.identity(3), coord.ICRS, coord.ICRS],
         [0, 0],
         [0],
         [0],
         [np.identity(3), coord.ICRS, coord.ICRS],
         ]

_kwargs = [{'unit': 'radian'},
           {'z': 0.23},
           {},
           {'unit': ['radian', 'radian']},
           {'unit': 'radian'},
           {'unit': 'radian'},
           {},
           ]


@pytest.mark.parametrize(("name", "args", "kwargs", "xfail"),
                         zip(_names, _args, _kwargs, _xfail))
def test_simple_object(pickle_protocol, name, args, kwargs, xfail):
    # Tests easily instantiated objects
    if xfail:
        pytest.xfail()
    original = name(*args, **kwargs)
    check_pickling_recovery(original, pickle_protocol)
