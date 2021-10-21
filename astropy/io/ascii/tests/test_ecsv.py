# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some of the methods related to the ``ECSV``
reader/writer.
"""
from astropy.table.column import MaskedColumn
import os
import copy
import sys
from io import StringIO
from contextlib import nullcontext

import pytest
import numpy as np
import yaml

from astropy.table import Table, Column, QTable, NdarrayMixin
from astropy.table.table_helpers import simple_table
from astropy.coordinates import (SkyCoord, Latitude, Longitude, Angle, EarthLocation,
                                 SphericalRepresentation, CartesianRepresentation,
                                 SphericalCosLatDifferential)
from astropy.time import Time, TimeDelta
from astropy.units import allclose as quantity_allclose
from astropy.units import QuantityInfo

from astropy.utils.exceptions import AstropyUserWarning

from astropy.io.ascii.ecsv import DELIMITERS
from astropy.io import ascii
from astropy import units as u

from .common import TEST_DIR

DTYPES = ['bool', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32',
          'uint64', 'float16', 'float32', 'float64', 'float128',
          'str']
if os.name == 'nt' or sys.maxsize <= 2**32:
    DTYPES.remove('float128')

T_DTYPES = Table()

for dtype in DTYPES:
    if dtype == 'bool':
        data = np.array([False, True, False])
    elif dtype == 'str':
        data = np.array(['ab 0', 'ab, 1', 'ab2'])
    else:
        data = np.arange(3, dtype=dtype)
    c = Column(data, unit='m / s', description='descr_' + dtype,
               meta={'meta ' + dtype: 1})
    T_DTYPES[dtype] = c

T_DTYPES.meta['comments'] = ['comment1', 'comment2']

# Corresponds to simple_table()
SIMPLE_LINES = ['# %ECSV 1.0',
                '# ---',
                '# datatype:',
                '# - {name: a, datatype: int64}',
                '# - {name: b, datatype: float64}',
                '# - {name: c, datatype: string}',
                '# schema: astropy-2.0',
                'a b c',
                '1 1.0 c',
                '2 2.0 d',
                '3 3.0 e']


def test_write_simple():
    """
    Write a simple table with common types.  This shows the compact version
    of serialization with one line per column.
    """
    t = simple_table()

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    assert out.getvalue().splitlines() == SIMPLE_LINES


def test_write_full():
    """
    Write a full-featured table with common types and explicitly checkout output
    """
    t = T_DTYPES['bool', 'int64', 'float64', 'str']
    lines = ['# %ECSV 1.0',
             '# ---',
             '# datatype:',
             '# - name: bool',
             '#   unit: m / s',
             '#   datatype: bool',
             '#   description: descr_bool',
             '#   meta: {meta bool: 1}',
             '# - name: int64',
             '#   unit: m / s',
             '#   datatype: int64',
             '#   description: descr_int64',
             '#   meta: {meta int64: 1}',
             '# - name: float64',
             '#   unit: m / s',
             '#   datatype: float64',
             '#   description: descr_float64',
             '#   meta: {meta float64: 1}',
             '# - name: str',
             '#   unit: m / s',
             '#   datatype: string',
             '#   description: descr_str',
             '#   meta: {meta str: 1}',
             '# meta: !!omap',
             '# - comments: [comment1, comment2]',
             '# schema: astropy-2.0',
             'bool int64 float64 str',
             'False 0 0.0 "ab 0"',
             'True 1 1.0 "ab, 1"',
             'False 2 2.0 ab2']

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    assert out.getvalue().splitlines() == lines


def test_write_read_roundtrip():
    """
    Write a full-featured table with all types and see that it round-trips on
    readback.  Use both space and comma delimiters.
    """
    t = T_DTYPES
    for delimiter in DELIMITERS:
        out = StringIO()
        t.write(out, format='ascii.ecsv', delimiter=delimiter)

        t2s = [Table.read(out.getvalue(), format='ascii.ecsv'),
               Table.read(out.getvalue(), format='ascii'),
               ascii.read(out.getvalue()),
               ascii.read(out.getvalue(), format='ecsv', guess=False),
               ascii.read(out.getvalue(), format='ecsv')]
        for t2 in t2s:
            assert t.meta == t2.meta
            for name in t.colnames:
                assert t[name].attrs_equal(t2[name])
                assert np.all(t[name] == t2[name])


def test_bad_delimiter():
    """
    Passing a delimiter other than space or comma gives an exception
    """
    out = StringIO()
    with pytest.raises(ValueError) as err:
        T_DTYPES.write(out, format='ascii.ecsv', delimiter='|')
    assert 'only space and comma are allowed' in str(err.value)


def test_bad_header_start():
    """
    Bad header without initial # %ECSV x.x
    """
    lines = copy.copy(SIMPLE_LINES)
    lines[0] = '# %ECV 0.9'
    with pytest.raises(ascii.InconsistentTableError):
        Table.read('\n'.join(lines), format='ascii.ecsv', guess=False)


def test_bad_delimiter_input():
    """
    Illegal delimiter in input
    """
    lines = copy.copy(SIMPLE_LINES)
    lines.insert(2, '# delimiter: |')
    with pytest.raises(ValueError) as err:
        Table.read('\n'.join(lines), format='ascii.ecsv', guess=False)
    assert 'only space and comma are allowed' in str(err.value)


def test_multidim_input():
    """
    Multi-dimensional column in input
    """
    t = Table()
    t['a'] = np.arange(24).reshape(2, 3, 4)
    t['a'].info.description = 'description'
    t['a'].info.meta = {1: 2}
    t['b'] = [1, 2]

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')

    assert np.all(t2['a'] == t['a'])
    assert t2['a'].shape == t['a'].shape
    assert t2['a'].dtype == t['a'].dtype
    assert t2['a'].info.description == t['a'].info.description
    assert t2['a'].info.meta == t['a'].info.meta

    assert np.all(t2['b'] == t['b'])


def test_round_trip_empty_table():
    """Test fix in #5010 for issue #5009 (ECSV fails for empty type with bool type)"""
    t = Table(dtype=[bool, 'i', 'f'], names=['a', 'b', 'c'])
    out = StringIO()
    t.write(out, format='ascii.ecsv')
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')
    assert t.dtype == t2.dtype
    assert len(t2) == 0


def test_csv_ecsv_colnames_mismatch():
    """
    Test that mismatch in column names from normal CSV header vs.
    ECSV YAML header raises the expected exception.
    """
    lines = copy.copy(SIMPLE_LINES)
    header_index = lines.index('a b c')
    lines[header_index] = 'a b d'
    with pytest.raises(ValueError) as err:
        ascii.read(lines, format='ecsv')
    assert "column names from ECSV header ['a', 'b', 'c']" in str(err.value)


def test_regression_5604():
    """
    See https://github.com/astropy/astropy/issues/5604 for more.
    """
    t = Table()
    t.meta = {"foo": 5 * u.km, "foo2": u.s}
    t["bar"] = [7] * u.km

    out = StringIO()
    t.write(out, format="ascii.ecsv")

    assert '!astropy.units.Unit' in out.getvalue()
    assert '!astropy.units.Quantity' in out.getvalue()


def assert_objects_equal(obj1, obj2, attrs, compare_class=True):
    if compare_class:
        assert obj1.__class__ is obj2.__class__

    # For a column that is a native astropy Column, ignore the specified
    # `attrs`. This happens for a mixin like Quantity that is stored in a
    # `Table` (not QTable).
    if isinstance(obj1, Column):
        attrs = []

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

    # For no attrs that means we just compare directly.
    if not attrs:
        if isinstance(obj1, np.ndarray) and obj1.dtype.kind == 'f':
            assert quantity_allclose(obj1, obj2, rtol=1e-15)
        else:
            assert np.all(obj1 == obj2)


# TODO: unify with the very similar tests in fits/tests/test_connect.py
# and misc/tests/test_hd5f.py.
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
scd = SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,m', frame='fk4',
               obstime=['J1990.5'] * 2)
scdc = scd.copy()
scdc.representation_type = 'cartesian'
scpm = SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,pc',
                pm_ra_cosdec=[7, 8]*u.mas/u.yr, pm_dec=[9, 10]*u.mas/u.yr)
scpmrv = SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,pc',
                  pm_ra_cosdec=[7, 8]*u.mas/u.yr, pm_dec=[9, 10]*u.mas/u.yr,
                  radial_velocity=[11, 12]*u.km/u.s)
scrv = SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,pc',
                radial_velocity=[11, 12]*u.km/u.s)
tm = Time([51000.5, 51001.5], format='mjd', scale='tai', precision=5, location=el[0])
tm2 = Time(tm, format='iso')
tm3 = Time(tm, location=el)
tm3.info.serialize_method['ecsv'] = 'jd1_jd2'
obj = Column([{'a': 1}, {'b': [2]}], dtype='object')

# NOTE: in the test below the name of the column "x" for the Quantity is
# important since it tests the fix for #10215 (namespace clash, where "x"
# clashes with "el.x").
mixin_cols = {
    'tm': tm,
    'tm2': tm2,
    'tm3': tm3,
    'dt': TimeDelta([1, 2] * u.day),
    'sc': sc,
    'scd': scd,
    'scdc': scdc,
    'scpm': scpm,
    'scpmrv': scpmrv,
    'scrv': scrv,
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
    'nd': NdarrayMixin([1, 2]),
    'obj': obj
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
    'scd': ['ra', 'dec', 'distance', 'representation_type', 'frame.name'],
    'scdc': ['x', 'y', 'z', 'representation_type', 'frame.name'],
    'scpm': ['ra', 'dec', 'distance', 'pm_ra_cosdec', 'pm_dec',
             'representation_type', 'frame.name'],
    'scpmrv': ['ra', 'dec', 'distance', 'pm_ra_cosdec', 'pm_dec',
               'radial_velocity', 'representation_type', 'frame.name'],
    'scrv': ['ra', 'dec', 'distance', 'radial_velocity', 'representation_type',
             'frame.name'],
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
    'obj': []
}


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
                        'obj',
                        'qdb',
                        'qdex',
                        'qmag',
                        'sc.ra', 'sc.dec',
                        'scd.ra', 'scd.dec', 'scd.distance',
                        'scd.obstime',
                        'scdc.x', 'scdc.y', 'scdc.z',
                        'scdc.obstime',
                        'scpm.ra', 'scpm.dec', 'scpm.distance',
                        'scpm.pm_ra_cosdec', 'scpm.pm_dec',
                        'scpmrv.ra', 'scpmrv.dec', 'scpmrv.distance',
                        'scpmrv.pm_ra_cosdec', 'scpmrv.pm_dec',
                        'scpmrv.radial_velocity',
                        'scrv.ra', 'scrv.dec', 'scrv.distance',
                        'scrv.radial_velocity',
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

    For the special case of ndim==1 just return the original.

    The output has shape [3] * ndim. By using 3 we can be sure that repeating
    the two input elements gives an output that is sufficiently unique for
    the multidim tests.
    """
    if ndim > 1:
        import itertools
        idxs = [idx for idx, _ in zip(itertools.cycle([0, 1]), range(3 ** ndim))]
        col = col[idxs].reshape([3] * ndim)
    return col


@pytest.mark.parametrize('name_col', list(mixin_cols.items()))
@pytest.mark.parametrize('table_cls', (Table, QTable))
@pytest.mark.parametrize('ndim', (1, 2, 3))
def test_ecsv_mixins_per_column(table_cls, name_col, ndim):
    """Test write/read one col at a time and do detailed validation.
    This tests every input column type as 1-d, 2-d and 3-d.
    """
    name, col = name_col

    c = make_multidim(np.array([1.0, 2.0]), ndim)
    col = make_multidim(col, ndim)
    t = table_cls([c, col, c], names=['c1', name, 'c2'])
    t[name].info.description = 'description'

    out = StringIO()
    t.write(out, format="ascii.ecsv")
    t2 = table_cls.read(out.getvalue(), format='ascii.ecsv')

    assert t.colnames == t2.colnames

    for colname in t.colnames:
        assert len(t2[colname].shape) == ndim
        assert_objects_equal(t[colname], t2[colname], compare_attrs[colname])

    # Special case to make sure Column type doesn't leak into Time class data
    if name.startswith('tm'):
        assert t2[name]._time.jd1.__class__ is np.ndarray
        assert t2[name]._time.jd2.__class__ is np.ndarray


def test_round_trip_masked_table_default(tmpdir):
    """Test (mostly) round-trip of MaskedColumn through ECSV using default serialization
    that uses an empty string "" to mark NULL values.  Note:

    >>> simple_table(masked=True)
    <Table masked=True length=3>
      a      b     c
    int64 float64 str1
    ----- ------- ----
       --     1.0    c
        2     2.0   --
        3      --    e
    """
    filename = str(tmpdir.join('test.ecsv'))

    t = simple_table(masked=True)  # int, float, and str cols with one masked element
    t.write(filename)

    t2 = Table.read(filename)
    assert t2.masked is False
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        # From formal perspective the round-trip columns are the "same"
        assert np.all(t2[name].mask == t[name].mask)
        assert np.all(t2[name] == t[name])

        # But peeking under the mask shows that the underlying data are changed
        # because by default ECSV uses "" to represent masked elements.
        t[name].mask = False
        t2[name].mask = False
        assert not np.all(t2[name] == t[name])  # Expected diff


def test_round_trip_masked_table_serialize_mask(tmpdir):
    """Same as prev but set the serialize_method to 'data_mask' so mask is written out"""
    filename = str(tmpdir.join('test.ecsv'))

    t = simple_table(masked=True)  # int, float, and str cols with one masked element
    t['c'][0] = ''  # This would come back as masked for default "" NULL marker

    # MaskedColumn with no masked elements. See table the MaskedColumnInfo class
    # _represent_as_dict() method for info about how we test a column with no masked elements.
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


@pytest.mark.parametrize('table_cls', (Table, QTable))
def test_ecsv_round_trip_user_defined_unit(table_cls, tmpdir):
    """Ensure that we can read-back enabled user-defined units."""

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


def test_read_masked_bool():
    txt = """\
# %ECSV 1.0
# ---
# datatype:
# - {name: col0, datatype: bool}
# schema: astropy-2.0
col0
1
0
True
""
False
"""
    dat = ascii.read(txt, format='ecsv')
    col = dat['col0']
    assert isinstance(col, MaskedColumn)
    assert np.all(col.mask == [False, False, False, True, False])
    assert np.all(col == [True, False, True, False, False])


@pytest.mark.parametrize('serialize_method', ['null_value', 'data_mask'])
@pytest.mark.parametrize('dtype', [np.int64, np.float64, bool, str])
@pytest.mark.parametrize('delimiter', [',', ' '])
def test_roundtrip_multidim_masked_array(serialize_method, dtype, delimiter):
    # TODO also test empty string with null value
    t = Table()
    col = MaskedColumn(np.arange(12).reshape(2, 3, 2), dtype=dtype)
    if dtype is str:
        # np does something funny and gives a dtype of U21.
        col = col.astype('U2')
    col.mask[0, 0, 0] = True
    col.mask[1, 1, 1] = True
    t['a'] = col
    t['b'] = ['x', 'y']  # Add another column for kicks
    out = StringIO()
    t.write(out, format='ascii.ecsv', serialize_method=serialize_method)
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')

    assert t2.masked is False
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        assert t2[name].dtype == t[name].dtype
        if hasattr(t[name], 'mask'):
            assert np.all(t2[name].mask == t[name].mask)
        assert np.all(t2[name] == t[name])


@pytest.mark.parametrize('subtype', ['some-user-type', 'complex'])
def test_multidim_unknown_subtype(subtype):
    """Test an ECSV file with a string type but unknown subtype"""
    txt = f"""\
# %ECSV 1.0
# ---
# datatype:
# - name: a
#   datatype: string
#   subtype: {subtype}
# schema: astropy-2.0
a
[1,2]
[3,4]"""
    with pytest.warns(AstropyUserWarning,
                      match=rf"unexpected subtype '{subtype}' set for column 'a'"):
        t = ascii.read(txt, format='ecsv')

    assert t['a'].dtype.kind == 'U'
    assert t['a'][0] == '[1,2]'


def test_multidim_bad_shape():
    """Test a malformed ECSV file"""
    txt = """\
# %ECSV 1.0
# ---
# datatype:
# - name: a
#   datatype: string
#   subtype: int64[3]
# schema: astropy-2.0
a
[1,2]
[3,4]"""
    with pytest.raises(ValueError, match="column 'a' failed to convert: shape mismatch"):
        Table.read(txt, format='ascii.ecsv')


def test_write_not_json_serializable():
    t = Table()
    t['a'] = np.array([set([1, 2]), 1], dtype=object)
    match = "could not convert column 'a' to string: Object of type set is not JSON serializable"
    out = StringIO()
    with pytest.raises(TypeError, match=match):
        t.write(out, format='ascii.ecsv')


def test_read_not_json_serializable():
    """Test a malformed ECSV file"""
    txt = """\
# %ECSV 1.0
# ---
# datatype:
# - {name: a, datatype: string, subtype: json}
# schema: astropy-2.0
a
fail
[3,4]"""
    match = "column 'a' failed to convert: column value is not valid JSON"
    with pytest.raises(ValueError, match=match):
        Table.read(txt, format='ascii.ecsv')


def test_read_complex():
    """Test an ECSV file with a complex column"""
    txt = """\
# %ECSV 1.0
# ---
# datatype:
# - {name: a, datatype: complex}
# schema: astropy-2.0
a
1+1j
2+2j"""
    match = "datatype 'complex' of column 'a' is not in allowed values"
    with pytest.raises(ValueError, match=match):
        Table.read(txt, format='ascii.ecsv')


def test_read_bad_datatype_for_object_subtype():
    """Test a malformed ECSV file"""
    txt = """\
# %ECSV 1.0
# ---
# datatype:
# - {name: a, datatype: int64, subtype: json}
# schema: astropy-2.0
a
fail
[3,4]"""
    match = "column 'a' failed to convert: datatype of column 'a' must be \"string\""
    with pytest.raises(ValueError, match=match):
        Table.read(txt, format='ascii.ecsv')


def test_read_bad_datatype():
    """Test a malformed ECSV file"""
    txt = """\
# %ECSV 1.0
# ---
# datatype:
# - {name: a, datatype: object}
# schema: astropy-2.0
a
fail
[3,4]"""
    match = r"column 'a' is not in allowed values \('bool', 'int8', 'int16', 'int32'"
    with pytest.raises(ValueError, match=match):
        Table.read(txt, format='ascii.ecsv')


def test_full_repr_roundtrip():
    """Test round-trip of float values to full precision even with format
    specified"""
    t = Table()
    t['a'] = np.array([np.pi, 1/7], dtype=np.float64)
    t['a'].info.format = '.2f'
    out = StringIO()
    t.write(out, format='ascii.ecsv')
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')
    assert np.all(t['a'] == t2['a'])
    assert t2['a'].info.format == '.2f'


#############################################################################
# Define a number of specialized columns for testing and the expected values
# of `datatype` for each column.
#############################################################################

# First here is some helper code used to make the expected outputs code.
def _get_ecsv_header_dict(text):
    lines = [line.strip() for line in text.splitlines()]
    lines = [line[2:] for line in lines if line.startswith('#')]
    lines = lines[2:]  # Get rid of the header
    out = yaml.safe_load('\n'.join(lines))
    return out


def _make_expected_values(cols):
    from pprint import pformat
    for name, col in cols.items():
        t = Table()
        t[name] = col
        out = StringIO()
        t.write(out, format='ascii.ecsv')
        hdr = _get_ecsv_header_dict(out.getvalue())
        fmt_hdr = pformat(hdr['datatype'])
        print(f'exps[{name!r}] =', fmt_hdr[:1])
        print(fmt_hdr[1:])
        print()


# Expected values of `datatype` for each column
exps = {}
cols = {}

# Run of the mill scalar for completeness
cols['scalar'] = np.array([1, 2], dtype=np.int16)
exps['scalar'] = [
    {'datatype': 'int16', 'name': 'scalar'}]

# Array of lists that works as a 2-d variable array. This is just treated
# as an object.
cols['2-d variable array lists'] = c = np.empty(shape=(2,), dtype=object)
c[0] = [[1, 2], ["a", 4]]
c[1] = [[1, 2, 3], [4, 5.25, 6]]
exps['2-d variable array lists'] = [
    {'datatype': 'string',
     'name': '2-d variable array lists',
     'subtype': 'json'}]

# Array of numpy arrays that is a 2-d variable array
cols['2-d variable array numpy'] = c = np.empty(shape=(2,), dtype=object)
c[0] = np.array([[1, 2], [3, 4]], dtype=np.float32)
c[1] = np.array([[1, 2, 3], [4, 5.5, 6]], dtype=np.float32)
exps['2-d variable array numpy'] = [
    {'datatype': 'string',
     'name': '2-d variable array numpy',
     'subtype': 'float32[2,null]'}]

cols['1-d variable array lists'] = np.array([[1, 2], [3, 4, 5]], dtype=object)
exps['1-d variable array lists'] = [
    {'datatype': 'string',
     'name': '1-d variable array lists',
     'subtype': 'json'}]

# Variable-length array
cols['1-d variable array numpy'] = np.array(
    [np.array([1, 2], dtype=np.uint8),
     np.array([3, 4, 5], dtype=np.uint8)], dtype=object)
exps['1-d variable array numpy'] = [
    {'datatype': 'string',
     'name': '1-d variable array numpy',
     'subtype': 'uint8[null]'}]

cols['1-d variable array numpy str'] = np.array(
    [np.array(['a', 'b']),
     np.array(['c', 'd', 'e'])], dtype=object)
exps['1-d variable array numpy str'] = [
    {'datatype': 'string',
     'name': '1-d variable array numpy str',
     'subtype': 'string[null]'}]

cols['1-d variable array numpy bool'] = np.array(
    [np.array([True, False]),
     np.array([True, False, True])], dtype=object)
exps['1-d variable array numpy bool'] = [
    {'datatype': 'string',
     'name': '1-d variable array numpy bool',
     'subtype': 'bool[null]'}]

cols['1-d regular array'] = np.array([[1, 2], [3, 4]], dtype=np.int8)
exps['1-d regular array'] = [
    {'datatype': 'string',
     'name': '1-d regular array',
     'subtype': 'int8[2]'}]

cols['2-d regular array'] = np.arange(8, dtype=np.float16).reshape(2, 2, 2)
exps['2-d regular array'] = [
    {'datatype': 'string',
     'name': '2-d regular array',
     'subtype': 'float16[2,2]'}]

cols['scalar object'] = np.array([{'a': 1}, {'b':2}], dtype=object)
exps['scalar object'] = [
    {'datatype': 'string', 'name': 'scalar object', 'subtype': 'json'}]

cols['1-d object'] = np.array(
    [[{'a': 1}, {'b':2}],
     [{'a': 1}, {'b':2}]], dtype=object)
exps['1-d object'] = [
    {'datatype': 'string',
     'name': '1-d object',
     'subtype': 'json[2]'}]


@pytest.mark.parametrize('name,col,exp',
                         list(zip(cols, cols.values(), exps.values())))
def test_specialized_columns(name, col, exp):
    """Test variable length lists, multidim columns, object columns.
    """
    t = Table()
    t[name] = col
    out = StringIO()
    t.write(out, format='ascii.ecsv')
    hdr = _get_ecsv_header_dict(out.getvalue())
    assert hdr['datatype'] == exp
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        assert t2[name].dtype == t[name].dtype
        for val1, val2 in zip(t2[name], t[name]):
            if isinstance(val1, np.ndarray):
                assert val1.dtype == val2.dtype
            assert np.all(val1 == val2)


def test_full_subtypes():
    """Read ECSV file created by M. Taylor that includes scalar, fixed array,
    variable array for all datatypes. This file has missing values for all
    columns as both per-value null and blank entries for the entire column
    value.

    Note: original file was modified to include blank values in f_float and
    f_double columns.
    """
    t = Table.read(os.path.join(TEST_DIR, 'data', 'subtypes.ecsv'))
    colnames = ('i_index,'
                's_byte,s_short,s_int,s_long,s_float,s_double,s_string,s_boolean,'
                'f_byte,f_short,f_int,f_long,f_float,f_double,f_string,f_boolean,'
                'v_byte,v_short,v_int,v_long,v_float,v_double,v_string,v_boolean,'
                'm_int,m_double').split(',')
    assert t.colnames == colnames

    type_map = {'byte': 'int8',
                'short': 'int16',
                'int': 'int32',
                'long': 'int64',
                'float': 'float32',
                'double': 'float64',
                'string': 'str',
                'boolean': 'bool'}

    for col in t.itercols():
        info = col.info
        if info.name == 'i_index':
            continue

        assert isinstance(col, MaskedColumn)

        type_name = info.name[2:]  # short, int, etc
        subtype = info.name[:1]

        if subtype == 's':  # Scalar
            assert col.shape == (16,)

        if subtype == 'f':  # Fixed array
            assert col.shape == (16, 3)

        if subtype == 'v':  # Variable array
            assert col.shape == (16,)
            assert info.dtype.name == 'object'
            for val in col:
                assert isinstance(val, np.ndarray)
                assert val.dtype.name.startswith(type_map[type_name])
                assert len(val) in [0, 1, 2, 3]
        else:
            assert info.dtype.name.startswith(type_map[type_name])


def test_masked_empty_subtypes():
    """Test blank field in subtypes. Similar to previous test but with explicit
    checks of values"""
    txt = """
    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: o, datatype: string, subtype: json}
    # - {name: f, datatype: string, subtype: 'int64[2]'}
    # - {name: v, datatype: string, subtype: 'int64[null]'}
    # schema: astropy-2.0
    o f v
    null [0,1] [1]
    "" "" ""
    [1,2] [2,3] [2,3]
    """
    t = Table.read(txt, format='ascii.ecsv')
    assert np.all(t['o'] == np.array([None, -1, [1, 2]], dtype=object))
    assert np.all(t['o'].mask == [False, True, False])

    exp = np.ma.array([[0, 1], [-1, -1], [2, 3]], mask=[[0, 0], [1, 1], [0, 0]])
    assert np.all(t['f'] == exp)
    assert np.all(t['f'].mask == exp.mask)

    assert np.all(t['v'][0] == [1])
    assert np.all(t['v'][2] == [2, 3])
    assert np.all(t['v'].mask == [False, True, False])


def test_masked_vals_in_array_subtypes():
    """Test null values in fixed and variable array subtypes."""
    t = Table()
    t['f'] = np.ma.array([[1, 2], [3, 4]], mask=[[0, 1], [1, 0]], dtype=np.int64)
    t['v'] = np.empty(2, dtype=object)
    t['v'][0] = np.ma.array([1, 2], mask=[0, 1], dtype=np.int64)
    t['v'][1] = np.ma.array([3, 4, 5], mask=[1, 0, 0], dtype=np.int64)

    out = StringIO()
    t.write(out, format='ascii.ecsv')
    txt = """
    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: f, datatype: string, subtype: 'int64[2]'}
    # - {name: v, datatype: string, subtype: 'int64[null]'}
    # schema: astropy-2.0
    f v
    [1,null] [1,null]
    [null,4] [null,4,5]
    """
    hdr = _get_ecsv_header_dict(out.getvalue())
    hdr_exp = _get_ecsv_header_dict(txt)
    assert hdr == hdr_exp
    t2 = Table.read(out.getvalue(), format='ascii.ecsv')
    assert t2.colnames == t.colnames
    for name in t2.colnames:
        assert t2[name].dtype == t[name].dtype
        assert type(t2[name]) is type(t[name])
        for val1, val2 in zip(t2[name], t[name]):
            if isinstance(val1, np.ndarray):
                assert val1.dtype == val2.dtype
            if isinstance(val1, np.ma.MaskedArray):
                assert np.all(val1.mask == val2.mask)
            assert np.all(val1 == val2)


def test_guess_ecsv_with_one_column():
    """Except for ECSV, guessing always requires at least 2 columns"""
    txt = """
    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: col, datatype: string, description: hello}
    # schema: astropy-2.0
    col
    1
    2
    """
    t = ascii.read(txt)
    assert t['col'].dtype.kind == 'U'  # would be int with basic format
    assert t['col'].description == 'hello'
