import pickle
import pytest
import numpy as np

import astropy.units as u
from astropy.coordinates import Longitude, representation
from astropy import coordinates as coord
from astropy.tests.helper import pickle_protocol, check_pickling_recovery  # noqa

# Can't test distances without scipy due to cosmology deps
from astropy.utils.compat.optional_deps import HAS_SCIPY  # noqa


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
                         tuple(zip(_names, _args, _kwargs, _xfail)))
def test_simple_object(pickle_protocol, name, args, kwargs, xfail):
    # Tests easily instantiated objects
    if xfail:
        pytest.xfail()
    original = name(*args, **kwargs)
    check_pickling_recovery(original, pickle_protocol)


class _CustomICRS(coord.ICRS):
    default_representation = coord.PhysicsSphericalRepresentation


@pytest.mark.parametrize(
    "frame",
    [
        coord.SkyOffsetFrame(origin=coord.ICRS(0*u.deg, 0*u.deg)),
        coord.SkyOffsetFrame(5*u.deg, 10*u.deg, origin=coord.Galactic(2*u.deg, -3*u.deg)),
        coord.SkyOffsetFrame(5*u.deg, 10*u.deg, 10*u.pc,
                             origin=coord.Galactic(2*u.deg, -3*u.deg),
                             representation_type=coord.PhysicsSphericalRepresentation),
        coord.SkyOffsetFrame(5*u.deg, 10*u.deg, 0*u.pc,
                             origin=_CustomICRS(2*u.deg, 3*u.deg, 1*u.pc)),
    ]
)
def test_skyoffset_pickle(pickle_protocol, frame):
    """
    This is a regression test for issue #9249:
    https://github.com/astropy/astropy/issues/9249
    """
    check_pickling_recovery(frame, pickle_protocol)
