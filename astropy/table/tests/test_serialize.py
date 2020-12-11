# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to table and mixin serialization,
including ECSV, FITS, and HDF5 reader/writers.

Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""
from io import StringIO

import pytest
import numpy as np

from astropy.table import Table, Column, QTable, NdarrayMixin
from astropy.table.table_helpers import simple_table
from astropy.coordinates import (SkyCoord, Latitude, Longitude, Angle, EarthLocation,
                                 SphericalRepresentation, CartesianRepresentation,
                                 SphericalCosLatDifferential)
from astropy.time import Time, TimeDelta
from astropy.units import allclose as quantity_allclose
from astropy.units import QuantityInfo

from astropy.io import ascii
from astropy import units as u

try:
    import yaml  # noqa
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


def assert_objects_equal(obj1, obj2, attrs, compare_class=True):
    if compare_class:
        assert obj1.__class__ is obj2.__class__

    assert obj1.shape == obj2.shape

    info_attrs = ['info.name', 'info.format', 'info.unit', 'info.description']
    for attr in attrs + info_attrs:
        a1 = obj1
        a2 = obj2
        for subattr in attr.split('.'):
            try:
                a1 = getattr(a1, subattr)
                a2 = getattr(a2, subattr)
            except AttributeError:
                a1 = a1[subattr]
                a2 = a2[subattr]

        if isinstance(a1, np.ndarray) and a1.dtype.kind == 'f':
            assert quantity_allclose(a1, a2, rtol=1e-10)
        else:
            assert np.all(a1 == a2)


# TODO: unify with the very similar tests in fits/tests/test_connect.py.
el = EarthLocation(x=[1, 2] * u.km, y=[3, 4] * u.km, z=[5, 6] * u.km)
sr = SphericalRepresentation(
    [0, 1]*u.deg, [2, 3]*u.deg, 1*u.kpc)
cr = CartesianRepresentation(
    [0, 1]*u.pc, [4, 5]*u.pc, [8, 6]*u.pc)
sd = SphericalCosLatDifferential(
    [0, 1]*u.mas/u.yr, [0, 1]*u.mas/u.yr, 10*u.km/u.s)
srd = SphericalRepresentation(sr, differentials=sd)
sc = SkyCoord([1, 2], [3, 4], unit='deg,deg', frame='fk4',
              obstime='J1990.5')
scc = sc.copy()
scc.representation_type = 'cartesian'
tm = Time([51000.5, 51001.5], format='mjd', scale='tai', precision=5, location=el[0])
tm2 = Time(tm, format='iso')
tm3 = Time(tm, location=el)
tm3.info.serialize_method['ecsv'] = 'jd1_jd2'

# NOTE: in the test below the name of the column "x" for the Quantity is
# important since it tests the fix for #10215 (namespace clash, where "x"
# clashes with "el.x").
mixin_cols = {
    'tm': tm,
    'tm2': tm2,
    'tm3': tm3,
    'dt': TimeDelta([1, 2] * u.day),
    'sc': sc,
    'scc': scc,
    'scd': SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,m', frame='fk4',
                    obstime=['J1990.5'] * 2),
    'x': [1, 2] * u.m,
    'qdb': [10, 20] * u.dB(u.mW),
    'qdex': [4.5, 5.5] * u.dex(u.cm / u.s**2),
    'qmag': [21, 22] * u.ABmag,
    'lat': Latitude([1, 2] * u.deg),
    'lon': Longitude([1, 2] * u.deg, wrap_angle=180. * u.deg),
    'ang': Angle([1, 2] * u.deg),
    'el': el,
    'sr': sr,
    'cr': cr,
    'sd': sd,
    'srd': srd,
    'nd': NdarrayMixin([1, 2])
}

time_attrs = ['value', 'shape', 'format', 'scale', 'precision',
              'in_subfmt', 'out_subfmt', 'location']
compare_attrs = {
    'c1': ['data'],
    'c2': ['data'],
    'tm': time_attrs,
    'tm2': time_attrs,
    'tm3': time_attrs,
    'dt': ['shape', 'value', 'format', 'scale'],
    'sc': ['ra', 'dec', 'representation_type', 'frame.name'],
    'scc': ['x', 'y', 'z', 'representation_type', 'frame.name'],
    'scd': ['ra', 'dec', 'distance', 'representation_type', 'frame.name'],
    'x': ['value', 'unit'],
    'qdb': ['value', 'unit'],
    'qdex': ['value', 'unit'],
    'qmag': ['value', 'unit'],
    'lon': ['value', 'unit', 'wrap_angle'],
    'lat': ['value', 'unit'],
    'ang': ['value', 'unit'],
    'el': ['x', 'y', 'z', 'ellipsoid'],
    'nd': ['data'],
    'sr': ['lon', 'lat', 'distance'],
    'cr': ['x', 'y', 'z'],
    'sd': ['d_lon_coslat', 'd_lat', 'd_distance'],
    'srd': ['lon', 'lat', 'distance', 'differentials.s.d_lon_coslat',
            'differentials.s.d_lat', 'differentials.s.d_distance'],
}


@pytest.mark.skipif('not HAS_YAML')
def test_ecsv_mixins_ascii_read_class():
    """Ensure that ascii.read(ecsv_file) returns the correct class
    (QTable if any Quantity subclasses, Table otherwise).
    """
    # Make a table with every mixin type except Quantities
    t = QTable({name: col for name, col in mixin_cols.items()
                if not isinstance(col.info, QuantityInfo)})
    out = StringIO()
    t.write(out, format="ascii.ecsv")
    t2 = ascii.read(out.getvalue(), format='ecsv')
    assert type(t2) is Table

    # Add a single quantity column
    t['lon'] = mixin_cols['lon']

    out = StringIO()
    t.write(out, format="ascii.ecsv")
    t2 = ascii.read(out.getvalue(), format='ecsv')
    assert type(t2) is QTable


@pytest.mark.skipif('not HAS_YAML')
def test_ecsv_mixins_qtable_to_table():
    """Test writing as QTable and reading as Table.  Ensure correct classes
    come out.
    """
    names = sorted(mixin_cols)

    t = QTable([mixin_cols[name] for name in names], names=names)
    out = StringIO()
    t.write(out, format="ascii.ecsv")
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')

    assert t.colnames == t2.colnames

    for name, col in t.columns.items():
        col2 = t2[name]
        attrs = compare_attrs[name]
        compare_class = True

        if isinstance(col.info, QuantityInfo):
            # Downgrade Quantity to Column + unit
            assert type(col2) is Column
            # Class-specific attributes like `value` or `wrap_angle` are lost.
            attrs = ['unit']
            compare_class = False
            # Compare data values here (assert_objects_equal doesn't know how in this case)
            assert np.allclose(col.value, col2, rtol=1e-10)

        assert_objects_equal(col, col2, attrs, compare_class)


@pytest.mark.skipif('not HAS_YAML')
@pytest.mark.parametrize('table_cls', (Table, QTable))
def test_ecsv_mixins_as_one(table_cls):
    """Test write/read all cols at once and validate intermediate column names"""
    names = sorted(mixin_cols)

    serialized_names = ['ang',
                        'cr.x', 'cr.y', 'cr.z',
                        'dt',
                        'el.x', 'el.y', 'el.z',
                        'lat',
                        'lon',
                        'nd',
                        'qdb',
                        'qdex',
                        'qmag',
                        'sc.ra', 'sc.dec',
                        'scc.x', 'scc.y', 'scc.z',
                        'scd.ra', 'scd.dec', 'scd.distance',
                        'scd.obstime',
                        'sd.d_lon_coslat', 'sd.d_lat', 'sd.d_distance',
                        'sr.lon', 'sr.lat', 'sr.distance',
                        'srd.lon', 'srd.lat', 'srd.distance',
                        'srd.differentials.s.d_lon_coslat',
                        'srd.differentials.s.d_lat',
                        'srd.differentials.s.d_distance',
                        'tm',  # serialize_method is formatted_value
                        'tm2',  # serialize_method is formatted_value
                        'tm3.jd1', 'tm3.jd2',    # serialize is jd1_jd2
                        'tm3.location.x', 'tm3.location.y', 'tm3.location.z',
                        'x']

    t = table_cls([mixin_cols[name] for name in names], names=names)

    out = StringIO()
    t.write(out, format="ascii.ecsv")
    t2 = table_cls.read(out.getvalue(), format='ascii.ecsv')

    assert t.colnames == t2.colnames

    # Read as a ascii.basic table (skip all the ECSV junk)
    t3 = table_cls.read(out.getvalue(), format='ascii.basic')
    assert t3.colnames == serialized_names


def make_multidim(col, ndim):
    """Take a col with length=2 and make it N-d by repeating elements.

    For the special case of ndim==1 just return the orignal.

    The output has shape [3] * ndim. By using 3 we can be sure that repeating
    the two input elements gives an output that is sufficiently unique for
    the multidim tests.
    """
    if ndim > 1:
        import itertools
        idxs = [idx for idx, _ in zip(itertools.cycle([0, 1]), range(3 ** ndim))]
        col = col[idxs].reshape([3] * ndim)
    return col


@pytest.mark.skipif('not HAS_YAML')
@pytest.mark.parametrize('name_col', list(mixin_cols.items()))
@pytest.mark.parametrize('table_cls', (Table, QTable))
@pytest.mark.parametrize('ndim', (1, 2, 3))
def test_ecsv_mixins_per_column(table_cls, name_col, ndim):
    """Test write/read one col at a time and do detailed validation"""
    name, col = name_col

    c = make_multidim(np.array([1.0, 2.0]), ndim)
    col = make_multidim(col, ndim)
    t = table_cls([c, col, c], names=['c1', name, 'c2'])
    t[name].info.description = 'description'

    if not t.has_mixin_columns:
        pytest.skip('column is not a mixin (e.g. Quantity subclass in Table)')

    out = StringIO()
    t.write(out, format="ascii.ecsv")
    t_raw = Table.read(out.getvalue(), format='ascii.basic', delimiter=' ', guess=False)
    assert len(t_raw.colnames) >= 3 * 3**(ndim - 1)  # 3 columns, each with 3**(ndim-1) subcols
    t2 = table_cls.read(out.getvalue(), format='ascii.ecsv')

    assert t.colnames == t2.colnames

    for colname in t.colnames:
        assert len(t2[colname].shape) == ndim
        assert_objects_equal(t[colname], t2[colname], compare_attrs[colname])

    # Special case to make sure Column type doesn't leak into Time class data
    if name.startswith('tm'):
        assert t2[name]._time.jd1.__class__ is np.ndarray
        assert t2[name]._time.jd2.__class__ is np.ndarray


@pytest.mark.skipif('not HAS_YAML')
def test_round_trip_masked_table_serialize_mask(tmpdir):
    """Same as prev but set the serialize_method to 'data_mask' so mask is written out"""
    filename = str(tmpdir.join('test.ecsv'))

    t = simple_table(masked=True)  # int, float, and str cols with one masked element
    t['c'][0] = ''  # This would come back as masked for default "" NULL marker

    # MaskedColumn with no masked elements. See table the MaskedColumnInfo class
    # _represent_as_dict() method for info about we test a column with no masked elements.
    t['d'] = [1, 2, 3]

    t.write(filename, serialize_method='data_mask')

    t2 = Table.read(filename)
    assert t2.masked is False
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        assert np.all(t2[name].mask == t[name].mask)
        assert np.all(t2[name] == t[name])

        # Data under the mask round-trips also (unmask data to show this).
        t[name].mask = False
        t2[name].mask = False
        assert np.all(t2[name] == t[name])


@pytest.mark.skipif('not HAS_YAML')
@pytest.mark.parametrize('table_cls', (Table, QTable))
def test_ecsv_round_trip_user_defined_unit(table_cls, tmpdir):
    """Ensure that we can read-back enabled user-defined units."""
    from astropy.utils.compat.context import nullcontext

    # Test adapted from #8897, where it was noted that this works
    # but was not tested.
    filename = str(tmpdir.join('test.ecsv'))
    unit = u.def_unit('bandpass_sol_lum')
    t = table_cls()
    t['l'] = np.arange(5) * unit
    t.write(filename)
    # without the unit enabled, get UnrecognizedUnit
    if table_cls is QTable:
        ctx = pytest.warns(u.UnitsWarning, match=r"'bandpass_sol_lum' did not parse .*")
    else:
        ctx = nullcontext()
    # Note: The read might also generate ResourceWarning, in addition to UnitsWarning
    with ctx:
        t2 = table_cls.read(filename)
    assert isinstance(t2['l'].unit, u.UnrecognizedUnit)
    assert str(t2['l'].unit) == 'bandpass_sol_lum'
    if table_cls is QTable:
        assert np.all(t2['l'].value == t['l'].value)
    else:
        assert np.all(t2['l'] == t['l'])

    # But with it enabled, it works.
    with u.add_enabled_units(unit):
        t3 = table_cls.read(filename)
        assert t3['l'].unit is unit
        assert np.all(t3['l'] == t['l'])

        # Just to be sure, aloso try writing with unit enabled.
        filename2 = str(tmpdir.join('test2.ecsv'))
        t3.write(filename2)
        t4 = table_cls.read(filename)
        assert t4['l'].unit is unit
        assert np.all(t4['l'] == t['l'])
