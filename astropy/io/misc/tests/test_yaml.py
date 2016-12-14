# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to YAML serialization.

Requires `pyyaml <http://pyyaml.org/>`_ to be installed.
"""
import numpy as np

from ....coordinates import SkyCoord, EarthLocation, Angle, Longitude, Latitude
from .... import units as u
from ....time import Time
from ....table import QTable
from ....extern.six.moves import StringIO

from ....tests.helper import pytest

try:
    from ..yaml import load, load_all, dump
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

pytestmark = pytest.mark.skipif('not HAS_YAML')


@pytest.mark.parametrize('c', [u.m, u.m / u.s, u.hPa, u.dimensionless_unscaled])
def test_unit(c):
    cy = load(dump(c))
    if isinstance(c, u.CompositeUnit):
        assert c == cy
    else:
        assert c is cy


@pytest.mark.parametrize('c', [Angle('1 2 3', unit='deg'),
                               Longitude('1 2 3', unit='deg'),
                               Latitude('1 2 3', unit='deg'),
                               [[1], [3]] * u.m,
                               np.array([[1, 2], [3, 4]], order='F'),
                               np.array([[1, 2], [3, 4]], order='C'),
                               np.array([1, 2, 3, 4])[::2]])
def test_ndarray_subclasses(c):
    cy = load(dump(c))

    assert np.all(c == cy)
    assert c.shape == cy.shape
    assert type(c) is type(cy)

    cc = 'C_CONTIGUOUS'
    fc = 'F_CONTIGUOUS'
    if c.flags[cc] or c.flags[fc]:
        assert c.flags[cc] == cy.flags[cc]
        assert c.flags[fc] == cy.flags[fc]
    else:
        # Original was not contiguous but round-trip version
        # should be c-contig.
        assert cy.flags[cc]

    if hasattr(c, 'unit'):
        assert c.unit == cy.unit


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


@pytest.mark.parametrize('frame', ['fk4', 'altaz'])
def test_skycoord(frame):

    c = SkyCoord([[1,2],[3,4]], [[5,6], [7,8]],
                 unit='deg', frame=frame,
                 obstime=Time('2016-01-02'),
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


def compare_time(t, ty):
    assert type(t) is type(ty)
    assert np.all(t == ty)
    for attr in ('shape', 'jd1', 'jd2', 'format', 'scale', 'precision', 'in_subfmt',
                 'out_subfmt', 'location', 'delta_ut1_utc', 'delta_tdb_tt'):
        assert np.all(getattr(t, attr) == getattr(ty, attr))


def test_time():
    t = _get_time()
    ty = load(dump(t))
    compare_time(t, ty)


def test_timedelta():
    t = _get_time()
    dt = t - t + 0.1234556 * u.s
    dty = load(dump(dt))

    assert type(dt) is type(dty)
    for attr in ('shape', 'jd1', 'jd2', 'format', 'scale'):
        assert np.all(getattr(dt, attr) == getattr(dty, attr))


def test_load_all():
    t = _get_time()
    unit = u.m / u.s
    c = SkyCoord([[1,2],[3,4]], [[5,6], [7,8]],
                 unit='deg', frame='fk4',
                 obstime=Time('2016-01-02'),
                 location=EarthLocation(1000, 2000, 3000, unit=u.km))

    # Make a multi-document stream
    out = ('---\n' + dump(t)
           + '---\n' + dump(unit)
           + '---\n' + dump(c))

    ty, unity, cy = list(load_all(out))

    compare_time(t, ty)
    compare_coord(c, cy)
    assert unity == unit


@pytest.mark.skipif('not HAS_YAML')
def test_ecsv_astropy_objects_in_meta():
    """
    Test that astropy core objects in ``meta`` are serialized.
    """
    t = QTable([[1, 2] * u.m, [4, 5]], names=['a', 'b'])
    tm = _get_time()
    c = SkyCoord([[1,2],[3,4]], [[5,6], [7,8]],
                 unit='deg', frame='fk4',
                 obstime=Time('2016-01-02'),
                 location=EarthLocation(1000, 2000, 3000, unit=u.km))
    unit = u.m / u.s

    t.meta = {'tm': tm, 'c': c, 'unit': unit}
    out = StringIO()
    t.write(out, format='ascii.ecsv')
    t2 = QTable.read(out.getvalue(), format='ascii.ecsv')

    compare_time(tm, t2.meta['tm'])
    compare_coord(c, t2.meta['c'])
    assert t2.meta['unit'] == unit
