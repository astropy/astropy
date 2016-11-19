# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to YAML serialization.

Requires `pyyaml <http://pyyaml.org/>`_ to be installed.
"""
import numpy as np

from ...coordinates import SkyCoord, EarthLocation, Angle, Longitude, Latitude
from ... import units as u
from ...time import Time

from ...tests.helper import pytest

try:
    from ..yaml import load, dump
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

pytestmark = pytest.mark.skipif('not HAS_YAML')


@pytest.mark.parametrize('c', [u.m, u.m / u.s, u.hPa, u.dimensionless_unscaled])
def test_unit(c):
    cy = load(dump(c))
    assert c == cy


@pytest.mark.parametrize('c', [Angle('1 2 3', unit='deg'),
                               Longitude('1 2 3', unit='deg'),
                               Latitude('1 2 3', unit='deg'),
                               [[1], [3]] * u.m,
                               np.array([[1], [3]])])
def test_ndarray_subclasses(c):
    cy = load(dump(c))

    assert np.all(c == cy)
    assert c.shape == cy.shape
    assert type(c) is type(cy)
    if hasattr(c, 'unit'):
        assert c.unit == cy.unit


@pytest.mark.parametrize('frame', ['fk4', 'altaz'])
def test_skycoord(frame):
    def compare_coord(c, cy):
        assert c.shape == cy.shape
        assert c.frame.name == cy.frame.name

        assert list(c.get_frame_attr_names()) == list(cy.get_frame_attr_names())
        for attr in c.get_frame_attr_names():
            assert getattr(c, attr) == getattr(cy, attr)

        assert (list(c.representation_component_names) ==
                list(cy.representation_component_names))
        for name in c.representation_component_names:
            assert np.all(getattr(c, attr) == getattr(cy, attr))


    c = SkyCoord([[1,2],[3,4]], [[5,6], [7,8]],
                 unit='deg', frame=frame,
                 obstime=Time.now(),
                 location=EarthLocation(1000, 2000, 3000, unit=u.km))
    cy = load(dump(c))
    compare_coord(c, cy)


def _get_time():
    t = Time([[1],[2]], format='cxcsec',
             location=EarthLocation(1000, 2000, 3000, unit=u.km))
    t.format = 'iso'
    t.precision = 5
    t.delta_ut1_utc = np.array([[3.0], [4.0]])
    t.delta_tdb_tt = np.array([[5.0], [6.0]])
    t.out_subfmt = 'date_hm'

    return t


def test_time():
    t = _get_time()
    ty = load(dump(t))

    assert type(t) is type(ty)
    for attr in ('shape', 'jd1', 'jd2', 'format', 'scale', 'precision', 'in_subfmt',
                 'out_subfmt', 'location', 'delta_ut1_utc', 'delta_tdb_tt'):
        assert np.all(getattr(t, attr) == getattr(ty, attr))


def test_timedelta():
    t = _get_time()
    dt = t - t + 0.1234556 * u.s
    dty = load(dump(dt))

    assert type(dt) is type(dty)
    for attr in ('shape', 'jd1', 'jd2', 'format', 'scale'):
        assert np.all(getattr(dt, attr) == getattr(dty, attr))
