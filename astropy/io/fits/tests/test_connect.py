import os
import gc
import sys
import pathlib
import warnings

import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.io.fits.column import (_parse_tdisp_format, _fortran_to_python_format,
                                    python_to_tdisp)

from astropy.io.fits import HDUList, PrimaryHDU, BinTableHDU, ImageHDU

from astropy.io import fits

from astropy import units as u
from astropy.table import Table, QTable, NdarrayMixin, Column
from astropy.table.table_helpers import simple_table
from astropy.tests.helper import catch_warnings
from astropy.units.format.fits import UnitScaleError
from astropy.utils.exceptions import (AstropyUserWarning,
                                      AstropyDeprecationWarning)

from astropy.coordinates import SkyCoord, Latitude, Longitude, Angle, EarthLocation
from astropy.time import Time, TimeDelta
from astropy.units.quantity import QuantityInfo

try:
    import yaml  # pylint: disable=W0611  # noqa
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

DATA = os.path.join(os.path.dirname(__file__), 'data')


def equal_data(a, b):
    for name in a.dtype.names:
        if not np.all(a[name] == b[name]):
            return False
    return True


class TestSingleTable:

    def setup_class(self):
        self.data = np.array(list(zip([1, 2, 3, 4],
                                      ['a', 'b', 'c', 'd'],
                                      [2.3, 4.5, 6.7, 8.9])),
                             dtype=[('a', int), ('b', 'U1'), ('c', float)])

    def test_simple(self, tmpdir):
        filename = str(tmpdir.join('test_simple.fts'))
        t1 = Table(self.data)
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename)
        assert equal_data(t1, t2)

    def test_simple_pathlib(self, tmpdir):
        filename = pathlib.Path(str(tmpdir.join('test_simple.fit')))
        t1 = Table(self.data)
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename)
        assert equal_data(t1, t2)

    def test_simple_meta(self, tmpdir):
        filename = str(tmpdir.join('test_simple.fits'))
        t1 = Table(self.data)
        t1.meta['A'] = 1
        t1.meta['B'] = 2.3
        t1.meta['C'] = 'spam'
        t1.meta['comments'] = ['this', 'is', 'a', 'long', 'comment']
        t1.meta['HISTORY'] = ['first', 'second', 'third']
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename)
        assert equal_data(t1, t2)
        for key in t1.meta:
            if isinstance(t1.meta, list):
                for i in range(len(t1.meta[key])):
                    assert t1.meta[key][i] == t2.meta[key][i]
            else:
                assert t1.meta[key] == t2.meta[key]

    def test_simple_meta_conflicting(self, tmpdir):
        filename = str(tmpdir.join('test_simple.fits'))
        t1 = Table(self.data)
        t1.meta['ttype1'] = 'spam'
        with catch_warnings() as l:
            t1.write(filename, overwrite=True)
        assert len(l) == 1
        assert str(l[0].message).startswith(
            'Meta-data keyword ttype1 will be ignored since it conflicts with a FITS reserved keyword')

    def test_simple_noextension(self, tmpdir):
        """
        Test that file type is recognized without extension
        """
        filename = str(tmpdir.join('test_simple'))
        t1 = Table(self.data)
        t1.write(filename, overwrite=True, format='fits')
        t2 = Table.read(filename)
        assert equal_data(t1, t2)

    @pytest.mark.parametrize('table_type', (Table, QTable))
    def test_with_units(self, table_type, tmpdir):
        filename = str(tmpdir.join('test_with_units.fits'))
        t1 = table_type(self.data)
        t1['a'].unit = u.m
        t1['c'].unit = u.km / u.s
        t1.write(filename, overwrite=True)
        t2 = table_type.read(filename)
        assert equal_data(t1, t2)
        assert t2['a'].unit == u.m
        assert t2['c'].unit == u.km / u.s

    @pytest.mark.skipif('not HAS_YAML')
    def test_with_custom_units_qtable(self, tmpdir):
        # Test only for QTable - for Table's Column, new units are dropped
        # (as is checked in test_write_drop_nonstandard_units).
        filename = str(tmpdir.join('test_with_units.fits'))
        unit = u.def_unit('bandpass_sol_lum')
        t = QTable()
        t['l'] = np.ones(5) * unit
        with catch_warnings(AstropyUserWarning) as w:
            t.write(filename, overwrite=True)
        assert len(w) == 1
        assert 'bandpass_sol_lum' in str(w[0].message)
        # Just reading back, the data is fine but the unit is not recognized.
        with catch_warnings() as w:
            t2 = QTable.read(filename)
        assert isinstance(t2['l'].unit, u.UnrecognizedUnit)
        assert str(t2['l'].unit) == 'bandpass_sol_lum'
        assert len(w) == 1
        assert "'bandpass_sol_lum' did not parse" in str(w[0].message)
        assert np.all(t2['l'].value == t['l'].value)

        # But if we enable the unit, it should be recognized.
        with u.add_enabled_units(unit):
            t3 = QTable.read(filename)
            assert t3['l'].unit is unit
            assert equal_data(t3, t)

            # Regression check for #8897; write used to fail when a custom
            # unit was enabled.
            with catch_warnings(u.UnitsWarning) as w:
                t3.write(filename, overwrite=True)
            assert len(w) == 0

    @pytest.mark.parametrize('table_type', (Table, QTable))
    def test_with_format(self, table_type, tmpdir):
        filename = str(tmpdir.join('test_with_format.fits'))
        t1 = table_type(self.data)
        t1['a'].format = '{:5d}'
        t1['b'].format = '{:>20}'
        t1['c'].format = '{:6.2f}'
        t1.write(filename, overwrite=True)
        t2 = table_type.read(filename)
        assert equal_data(t1, t2)
        assert t2['a'].format == '{:5d}'
        assert t2['b'].format == '{:>20}'
        assert t2['c'].format == '{:6.2f}'

    def test_masked(self, tmpdir):
        filename = str(tmpdir.join('test_masked.fits'))
        t1 = Table(self.data, masked=True)
        t1.mask['a'] = [1, 0, 1, 0]
        t1.mask['b'] = [1, 0, 0, 1]
        t1.mask['c'] = [0, 1, 1, 0]
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename)
        assert t2.masked
        assert equal_data(t1, t2)
        assert np.all(t1['a'].mask == t2['a'].mask)
        # Disabled for now, as there is no obvious way to handle masking of
        # non-integer columns in FITS
        # TODO: Re-enable these tests if some workaround for this can be found
        # assert np.all(t1['b'].mask == t2['b'].mask)
        # assert np.all(t1['c'].mask == t2['c'].mask)

    def test_masked_nan(self, tmpdir):
        filename = str(tmpdir.join('test_masked_nan.fits'))
        data = np.array(list(zip([5.2, 8.4, 3.9, 6.3],
                                 [2.3, 4.5, 6.7, 8.9])),
                                dtype=[('a', np.float64), ('b', np.float32)])
        t1 = Table(data, masked=True)
        t1.mask['a'] = [1, 0, 1, 0]
        t1.mask['b'] = [1, 0, 0, 1]
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename)
        np.testing.assert_array_almost_equal(t2['a'], [np.nan, 8.4, np.nan, 6.3])
        np.testing.assert_array_almost_equal(t2['b'], [np.nan, 4.5, 6.7, np.nan])
        # assert t2.masked
        # t2.masked = false currently, as the only way to determine whether a table is masked
        # while reading is to check whether col.null is present. For float columns, col.null
        # is not initialized

    def test_read_from_fileobj(self, tmpdir):
        filename = str(tmpdir.join('test_read_from_fileobj.fits'))
        hdu = BinTableHDU(self.data)
        hdu.writeto(filename, overwrite=True)
        with open(filename, 'rb') as f:
            t = Table.read(f)
        assert equal_data(t, self.data)

    def test_read_with_nonstandard_units(self):
        hdu = BinTableHDU(self.data)
        hdu.columns[0].unit = 'RADIANS'
        hdu.columns[1].unit = 'spam'
        hdu.columns[2].unit = 'millieggs'
        t = Table.read(hdu)
        assert equal_data(t, self.data)

    @pytest.mark.parametrize('table_type', (Table, QTable))
    def test_write_drop_nonstandard_units(self, table_type, tmpdir):
        # While we are generous on input (see above), we are strict on
        # output, dropping units not recognized by the fits standard.
        filename = str(tmpdir.join('test_nonstandard_units.fits'))
        spam = u.def_unit('spam')
        t = table_type()
        t['a'] = [1., 2., 3.] * spam
        with catch_warnings() as w:
            t.write(filename)
        assert len(w) == 1
        assert 'spam' in str(w[0].message)
        if table_type is Table or not HAS_YAML:
            assert ('cannot be recovered in reading. '
                    'If pyyaml is installed') in str(w[0].message)
        else:
            assert 'lost to non-astropy fits readers' in str(w[0].message)

        with fits.open(filename) as ff:
            hdu = ff[1]
            assert 'TUNIT1' not in hdu.header

    def test_memmap(self, tmpdir):
        filename = str(tmpdir.join('test_simple.fts'))
        t1 = Table(self.data)
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename, memmap=False)
        t3 = Table.read(filename, memmap=True)
        assert equal_data(t2, t3)
        # To avoid issues with --open-files, we need to remove references to
        # data that uses memory mapping and force the garbage collection
        del t1, t2, t3
        gc.collect()

    @pytest.mark.parametrize('memmap', (False, True))
    def test_character_as_bytes(self, tmpdir, memmap):
        filename = str(tmpdir.join('test_simple.fts'))
        t1 = Table(self.data)
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename, character_as_bytes=False, memmap=memmap)
        t3 = Table.read(filename, character_as_bytes=True, memmap=memmap)
        assert t2['b'].dtype.kind == 'U'
        assert t3['b'].dtype.kind == 'S'
        assert equal_data(t2, t3)
        # To avoid issues with --open-files, we need to remove references to
        # data that uses memory mapping and force the garbage collection
        del t1, t2, t3
        gc.collect()

    def test_oned_single_element(self, tmpdir):
        filename = str(tmpdir.join('test_oned_single_element.fits'))
        table = Table({'x': [[1], [2]]})
        table.write(filename, overwrite=True)

        read = Table.read(filename)
        assert read['x'].shape == (2, 1)
        assert len(read['x'][0]) == 1


class TestMultipleHDU:

    def setup_class(self):
        self.data1 = np.array(list(zip([1, 2, 3, 4],
                                       ['a', 'b', 'c', 'd'],
                                       [2.3, 4.5, 6.7, 8.9])),
                              dtype=[('a', int), ('b', 'U1'), ('c', float)])
        self.data2 = np.array(list(zip([1.4, 2.3, 3.2, 4.7],
                                       [2.3, 4.5, 6.7, 8.9])),
                              dtype=[('p', float), ('q', float)])
        self.data3 = np.array(list(zip([1, 2, 3, 4],
                                       [2.3, 4.5, 6.7, 8.9])),
                              dtype=[('A', int), ('B', float)])
        hdu0 = PrimaryHDU()
        hdu1 = BinTableHDU(self.data1, name='first')
        hdu2 = BinTableHDU(self.data2, name='second')
        hdu3 = ImageHDU(np.ones((3, 3)), name='third')
        hdu4 = BinTableHDU(self.data3)

        self.hdus = HDUList([hdu0, hdu1, hdu2, hdu3, hdu4])
        self.hdusb = HDUList([hdu0, hdu3, hdu2, hdu1])
        self.hdus3 = HDUList([hdu0, hdu3, hdu2])
        self.hdus2 = HDUList([hdu0, hdu1, hdu3])
        self.hdus1 = HDUList([hdu0, hdu1])

    def teardown_class(self):
        del self.hdus

    def setup_method(self, method):
        warnings.filterwarnings('always')

    def test_read(self, tmpdir):
        filename = str(tmpdir.join('test_read.fits'))
        self.hdus.writeto(filename)
        with pytest.warns(AstropyUserWarning,
                          match=r"hdu= was not specified but multiple tables "
                                r"are present, reading in first available "
                                r"table \(hdu=1\)"):
            t = Table.read(filename)
        assert equal_data(t, self.data1)

        filename = str(tmpdir.join('test_read_2.fits'))
        self.hdusb.writeto(filename)
        with pytest.warns(AstropyUserWarning,
                          match=r"hdu= was not specified but multiple tables "
                                r"are present, reading in first available "
                                r"table \(hdu=2\)"):
            t3 = Table.read(filename)
        assert equal_data(t3, self.data2)

    def test_read_with_hdu_0(self, tmpdir):
        filename = str(tmpdir.join('test_read_with_hdu_0.fits'))
        self.hdus.writeto(filename)
        with pytest.raises(ValueError) as exc:
            Table.read(filename, hdu=0)
        assert exc.value.args[0] == 'No table found in hdu=0'

    @pytest.mark.parametrize('hdu', [1, 'first'])
    def test_read_with_hdu_1(self, tmpdir, hdu):
        filename = str(tmpdir.join('test_read_with_hdu_1.fits'))
        self.hdus.writeto(filename)
        with catch_warnings() as l:
            t = Table.read(filename, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data1)

    @pytest.mark.parametrize('hdu', [2, 'second'])
    def test_read_with_hdu_2(self, tmpdir, hdu):
        filename = str(tmpdir.join('test_read_with_hdu_2.fits'))
        self.hdus.writeto(filename)
        with catch_warnings() as l:
            t = Table.read(filename, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data2)

    @pytest.mark.parametrize('hdu', [3, 'third'])
    def test_read_with_hdu_3(self, tmpdir, hdu):
        filename = str(tmpdir.join('test_read_with_hdu_3.fits'))
        self.hdus.writeto(filename)
        with pytest.raises(ValueError, match='No table found in hdu=3'):
            Table.read(filename, hdu=hdu)

    def test_read_with_hdu_4(self, tmpdir):
        filename = str(tmpdir.join('test_read_with_hdu_4.fits'))
        self.hdus.writeto(filename)
        with catch_warnings() as l:
            t = Table.read(filename, hdu=4)
        assert len(l) == 0
        assert equal_data(t, self.data3)

    @pytest.mark.parametrize('hdu', [2, 3, '1', 'second', ''])
    def test_read_with_hdu_missing(self, tmpdir, hdu):
        filename = str(tmpdir.join('test_warn_with_hdu_1.fits'))
        self.hdus1.writeto(filename)
        with pytest.warns(AstropyDeprecationWarning,
                          match=rf"Specified hdu={hdu} not found, "
                                r"reading in first available table \(hdu=1\)"):
            t1 = Table.read(filename, hdu=hdu)
        assert equal_data(t1, self.data1)

    @pytest.mark.parametrize('hdu', [0, 2, 'third'])
    def test_read_with_hdu_warning(self, tmpdir, hdu):
        filename = str(tmpdir.join('test_warn_with_hdu_2.fits'))
        self.hdus2.writeto(filename)
        with pytest.warns(AstropyDeprecationWarning,
                          match=rf"No table found in specified hdu={hdu}, "
                                r"reading in first available table \(hdu=1\)"):
            t2 = Table.read(filename, hdu=hdu)
        assert equal_data(t2, self.data1)

    @pytest.mark.parametrize('hdu', [0, 1, 'third'])
    def test_read_in_last_hdu(self, tmpdir, hdu):
        filename = str(tmpdir.join('test_warn_with_hdu_3.fits'))
        self.hdus3.writeto(filename)
        with pytest.warns(AstropyDeprecationWarning,
                          match=rf"No table found in specified hdu={hdu}, "
                                r"reading in first available table \(hdu=2\)"):
            t3 = Table.read(filename, hdu=hdu)
        assert equal_data(t3, self.data2)

    def test_read_from_hdulist(self):
        with pytest.warns(AstropyUserWarning,
                          match=r"hdu= was not specified but multiple tables "
                                r"are present, reading in first available "
                                r"table \(hdu=1\)"):
            t = Table.read(self.hdus)
        assert equal_data(t, self.data1)

        with pytest.warns(AstropyUserWarning,
                          match=r"hdu= was not specified but multiple tables "
                                r"are present, reading in first available "
                                r"table \(hdu=2\)"):
            t3 = Table.read(self.hdusb)
        assert equal_data(t3, self.data2)

    def test_read_from_hdulist_with_hdu_0(self):
        with pytest.raises(ValueError) as exc:
            Table.read(self.hdus, hdu=0)
        assert exc.value.args[0] == 'No table found in hdu=0'

    @pytest.mark.parametrize('hdu', [1, 'first', None])
    def test_read_from_hdulist_with_single_table(self, hdu):
        with catch_warnings() as l:
            t = Table.read(self.hdus1, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data1)

    @pytest.mark.parametrize('hdu', [1, 'first'])
    def test_read_from_hdulist_with_hdu_1(self, hdu):
        with catch_warnings() as l:
            t = Table.read(self.hdus, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data1)

    @pytest.mark.parametrize('hdu', [2, 'second'])
    def test_read_from_hdulist_with_hdu_2(self, hdu):
        with catch_warnings() as l:
            t = Table.read(self.hdus, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data2)

    @pytest.mark.parametrize('hdu', [3, 'third'])
    def test_read_from_hdulist_with_hdu_3(self, hdu):
        with pytest.raises(ValueError, match='No table found in hdu=3'):
            Table.read(self.hdus, hdu=hdu)

    @pytest.mark.parametrize('hdu', [0, 2, 'third'])
    def test_read_from_hdulist_with_hdu_warning(self, hdu):
        with pytest.warns(AstropyDeprecationWarning,
                          match=rf"No table found in specified hdu={hdu}, "
                                r"reading in first available table \(hdu=1\)"):\
            t2 = Table.read(self.hdus2, hdu=hdu)
        assert equal_data(t2, self.data1)

    @pytest.mark.parametrize('hdu', [2, 3, '1', 'second', ''])
    def test_read_from_hdulist_with_hdu_missing(self, hdu):
        with pytest.warns(AstropyDeprecationWarning,
                          match=rf"Specified hdu={hdu} not found, "
                                r"reading in first available table \(hdu=1\)"):
            t1 = Table.read(self.hdus1, hdu=hdu)
        assert equal_data(t1, self.data1)

    @pytest.mark.parametrize('hdu', [0, 1, 'third'])
    def test_read_from_hdulist_in_last_hdu(self, hdu):
        with pytest.warns(AstropyDeprecationWarning,
                          match=rf"No table found in specified hdu={hdu}, "
                                r"reading in first available table \(hdu=2\)"):
            t3 = Table.read(self.hdus3, hdu=hdu)
        assert equal_data(t3, self.data2)

    @pytest.mark.parametrize('hdu', [None, 1, 'first'])
    def test_read_from_single_hdu(self, hdu):
        with catch_warnings() as l:
            t = Table.read(self.hdus[1])
        assert len(l) == 0
        assert equal_data(t, self.data1)


def test_masking_regression_1795():
    """
    Regression test for #1795 - this bug originally caused columns where TNULL
    was not defined to have their first element masked.
    """
    t = Table.read(os.path.join(DATA, 'tb.fits'))
    assert np.all(t['c1'].mask == np.array([False, False]))
    assert np.all(t['c2'].mask == np.array([False, False]))
    assert np.all(t['c3'].mask == np.array([False, False]))
    assert np.all(t['c4'].mask == np.array([False, False]))
    assert np.all(t['c1'].data == np.array([1, 2]))
    assert np.all(t['c2'].data == np.array([b'abc', b'xy ']))
    assert_allclose(t['c3'].data, np.array([3.70000007153, 6.6999997139]))
    assert np.all(t['c4'].data == np.array([False, True]))


def test_scale_error():
    a = [1, 4, 5]
    b = [2.0, 5.0, 8.2]
    c = ['x', 'y', 'z']
    t = Table([a, b, c], names=('a', 'b', 'c'), meta={'name': 'first table'})
    t['a'].unit = '1.2'
    with pytest.raises(UnitScaleError) as exc:
        t.write('t.fits', format='fits', overwrite=True)
    assert exc.value.args[0] == "The column 'a' could not be stored in FITS format because it has a scale '(1.2)' that is not recognized by the FITS standard. Either scale the data or change the units."


@pytest.mark.parametrize('tdisp_str, format_return',
                         [('EN10.5', ('EN', '10', '5', None)),
                          ('F6.2', ('F', '6', '2', None)),
                          ('B5.10', ('B', '5', '10', None)),
                          ('E10.5E3', ('E', '10', '5', '3')),
                          ('A21', ('A', '21', None, None))])
def test_parse_tdisp_format(tdisp_str, format_return):
    assert _parse_tdisp_format(tdisp_str) == format_return


@pytest.mark.parametrize('tdisp_str, format_str_return',
                         [('G15.4E2', '{:15.4g}'),
                          ('Z5.10', '{:5x}'),
                          ('I6.5', '{:6d}'),
                          ('L8', '{:>8}'),
                          ('E20.7', '{:20.7e}')])
def test_fortran_to_python_format(tdisp_str, format_str_return):
    assert _fortran_to_python_format(tdisp_str) == format_str_return


@pytest.mark.parametrize('fmt_str, tdisp_str',
                         [('{:3d}', 'I3'),
                          ('3d', 'I3'),
                          ('7.3f', 'F7.3'),
                          ('{:>4}', 'A4'),
                          ('{:7.4f}', 'F7.4'),
                          ('%5.3g', 'G5.3'),
                          ('%10s', 'A10'),
                          ('%.4f', 'F13.4')])
def test_python_to_tdisp(fmt_str, tdisp_str):
    assert python_to_tdisp(fmt_str) == tdisp_str


def test_logical_python_to_tdisp():
    assert python_to_tdisp('{:>7}', logical_dtype=True) == 'L7'


def test_bool_column(tmpdir):
    """
    Regression test for https://github.com/astropy/astropy/issues/1953

    Ensures that Table columns of bools are properly written to a FITS table.
    """

    arr = np.ones(5, dtype=bool)
    arr[::2] == np.False_

    t = Table([arr])
    t.write(str(tmpdir.join('test.fits')), overwrite=True)

    with fits.open(str(tmpdir.join('test.fits'))) as hdul:
        assert hdul[1].data['col0'].dtype == np.dtype('bool')
        assert np.all(hdul[1].data['col0'] == arr)


def test_unicode_column(tmpdir):
    """
    Test that a column of unicode strings is still written as one
    byte-per-character in the FITS table (so long as the column can be ASCII
    encoded).

    Regression test for one of the issues fixed in
    https://github.com/astropy/astropy/pull/4228
    """

    t = Table([np.array(['a', 'b', 'cd'])])
    t.write(str(tmpdir.join('test.fits')), overwrite=True)

    with fits.open(str(tmpdir.join('test.fits'))) as hdul:
        assert np.all(hdul[1].data['col0'] == ['a', 'b', 'cd'])
        assert hdul[1].header['TFORM1'] == '2A'

    t2 = Table([np.array(['\N{SNOWMAN}'])])

    with pytest.raises(UnicodeEncodeError):
        t2.write(str(tmpdir.join('test.fits')), overwrite=True)


def test_unit_warnings_read_write(tmpdir):
    filename = str(tmpdir.join('test_unit.fits'))
    t1 = Table([[1, 2], [3, 4]], names=['a', 'b'])
    t1['a'].unit = 'm/s'
    t1['b'].unit = 'not-a-unit'

    with catch_warnings() as l:
        t1.write(filename, overwrite=True)
        assert len(l) == 1
        assert str(l[0].message).startswith("'not-a-unit' did not parse as fits unit")

    with catch_warnings() as l:
        Table.read(filename, hdu=1)
    assert len(l) == 0


def test_convert_comment_convention(tmpdir):
    """
    Regression test for https://github.com/astropy/astropy/issues/6079
    """
    filename = os.path.join(DATA, 'stddata.fits')
    with pytest.warns(AstropyUserWarning, match=r'hdu= was not specified but '
                      r'multiple tables are present'):
        t = Table.read(filename)

    assert t.meta['comments'] == [
        '',
        ' *** End of mandatory fields ***',
        '',
        '',
        ' *** Column names ***',
        '',
        '',
        ' *** Column formats ***',
        ''
    ]


def assert_objects_equal(obj1, obj2, attrs, compare_class=True):
    if compare_class:
        assert obj1.__class__ is obj2.__class__

    info_attrs = ['info.name', 'info.format', 'info.unit', 'info.description', 'info.meta']
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

        # Mixin info.meta can None instead of empty OrderedDict(), #6720 would
        # fix this.
        if attr == 'info.meta':
            if a1 is None:
                a1 = {}
            if a2 is None:
                a2 = {}

        assert np.all(a1 == a2)

# Testing FITS table read/write with mixins.  This is mostly
# copied from ECSV mixin testing.


el = EarthLocation(x=1 * u.km, y=3 * u.km, z=5 * u.km)
el2 = EarthLocation(x=[1, 2] * u.km, y=[3, 4] * u.km, z=[5, 6] * u.km)
sc = SkyCoord([1, 2], [3, 4], unit='deg,deg', frame='fk4',
              obstime='J1990.5')
scc = sc.copy()
scc.representation_type = 'cartesian'
tm = Time([2450814.5, 2450815.5], format='jd', scale='tai', location=el)

# NOTE: in the test below the name of the column "x" for the Quantity is
# important since it tests the fix for #10215 (namespace clash, where "x"
# clashes with "el2.x").
mixin_cols = {
    'tm': tm,
    'dt': TimeDelta([1, 2] * u.day),
    'sc': sc,
    'scc': scc,
    'scd': SkyCoord([1, 2], [3, 4], [5, 6], unit='deg,deg,m', frame='fk4',
                    obstime=['J1990.5', 'J1991.5']),
    'x': [1, 2] * u.m,
    'lat': Latitude([1, 2] * u.deg),
    'lon': Longitude([1, 2] * u.deg, wrap_angle=180. * u.deg),
    'ang': Angle([1, 2] * u.deg),
    'el2': el2,
}

time_attrs = ['value', 'shape', 'format', 'scale', 'location']
compare_attrs = {
    'c1': ['data'],
    'c2': ['data'],
    'tm': time_attrs,
    'dt': ['shape', 'value', 'format', 'scale'],
    'sc': ['ra', 'dec', 'representation_type', 'frame.name'],
    'scc': ['x', 'y', 'z', 'representation_type', 'frame.name'],
    'scd': ['ra', 'dec', 'distance', 'representation_type', 'frame.name'],
    'x': ['value', 'unit'],
    'lon': ['value', 'unit', 'wrap_angle'],
    'lat': ['value', 'unit'],
    'ang': ['value', 'unit'],
    'el2': ['x', 'y', 'z', 'ellipsoid'],
    'nd': ['x', 'y', 'z'],
}


@pytest.mark.skipif('not HAS_YAML')
def test_fits_mixins_qtable_to_table(tmpdir):
    """Test writing as QTable and reading as Table.  Ensure correct classes
    come out.
    """
    filename = str(tmpdir.join('test_simple.fits'))

    names = sorted(mixin_cols)

    t = QTable([mixin_cols[name] for name in names], names=names)
    t.write(filename, format='fits')
    t2 = Table.read(filename, format='fits', astropy_native=True)

    assert t.colnames == t2.colnames

    for name, col in t.columns.items():
        col2 = t2[name]

        # Special-case Time, which does not yet support round-tripping
        # the format.
        if isinstance(col2, Time):
            col2.format = col.format

        attrs = compare_attrs[name]
        compare_class = True

        if isinstance(col.info, QuantityInfo):
            # Downgrade Quantity to Column + unit
            assert type(col2) is Column
            # Class-specific attributes like `value` or `wrap_angle` are lost.
            attrs = ['unit']
            compare_class = False
            # Compare data values here (assert_objects_equal doesn't know how in this case)
            assert np.all(col.value == col2)

        assert_objects_equal(col, col2, attrs, compare_class)


@pytest.mark.skipif('not HAS_YAML')
@pytest.mark.parametrize('table_cls', (Table, QTable))
def test_fits_mixins_as_one(table_cls, tmpdir):
    """Test write/read all cols at once and validate intermediate column names"""
    filename = str(tmpdir.join('test_simple.fits'))
    names = sorted(mixin_cols)

    serialized_names = ['ang',
                        'dt.jd1', 'dt.jd2',
                        'el2.x', 'el2.y', 'el2.z',
                        'lat',
                        'lon',
                        'sc.ra', 'sc.dec',
                        'scc.x', 'scc.y', 'scc.z',
                        'scd.ra', 'scd.dec', 'scd.distance',
                        'scd.obstime.jd1', 'scd.obstime.jd2',
                        'tm',  # serialize_method is formatted_value
                        'x',
                        ]

    t = table_cls([mixin_cols[name] for name in names], names=names)
    t.meta['C'] = 'spam'
    t.meta['comments'] = ['this', 'is', 'a', 'comment']
    t.meta['history'] = ['first', 'second', 'third']

    t.write(filename, format="fits")

    t2 = table_cls.read(filename, format='fits', astropy_native=True)
    assert t2.meta['C'] == 'spam'
    assert t2.meta['comments'] == ['this', 'is', 'a', 'comment']
    assert t2.meta['HISTORY'] == ['first', 'second', 'third']

    assert t.colnames == t2.colnames

    # Read directly via fits and confirm column names
    with fits.open(filename) as hdus:
        assert hdus[1].columns.names == serialized_names


@pytest.mark.skipif('not HAS_YAML')
@pytest.mark.parametrize('name_col', list(mixin_cols.items()))
@pytest.mark.parametrize('table_cls', (Table, QTable))
def test_fits_mixins_per_column(table_cls, name_col, tmpdir):
    """Test write/read one col at a time and do detailed validation"""
    filename = str(tmpdir.join('test_simple.fits'))
    name, col = name_col

    c = [1.0, 2.0]
    t = table_cls([c, col, c], names=['c1', name, 'c2'])
    t[name].info.description = 'my \n\n\n description'
    t[name].info.meta = {'list': list(range(50)), 'dict': {'a': 'b' * 200}}

    if not t.has_mixin_columns:
        pytest.skip('column is not a mixin (e.g. Quantity subclass in Table)')

    if isinstance(t[name], NdarrayMixin):
        pytest.xfail('NdarrayMixin not supported')

    t.write(filename, format="fits")
    t2 = table_cls.read(filename, format='fits', astropy_native=True)

    assert t.colnames == t2.colnames

    for colname in t.colnames:
        assert_objects_equal(t[colname], t2[colname], compare_attrs[colname])

    # Special case to make sure Column type doesn't leak into Time class data
    if name.startswith('tm'):
        assert t2[name]._time.jd1.__class__ is np.ndarray
        assert t2[name]._time.jd2.__class__ is np.ndarray


def test_warn_for_dropped_info_attributes(tmpdir, monkeypatch):
    # make sure that yaml cannot be imported if it is available
    monkeypatch.setitem(sys.modules, 'yaml', None)

    filename = str(tmpdir.join('test.fits'))
    t = Table([[1, 2]])
    t['col0'].info.description = 'hello'
    with catch_warnings() as warns:
        t.write(filename, overwrite=True)
    assert len(warns) == 1
    assert str(warns[0].message).startswith(
        "table contains column(s) with defined 'format'")


def test_error_for_mixins_but_no_yaml(tmpdir, monkeypatch):
    # make sure that yaml cannot be imported if it is available
    monkeypatch.setitem(sys.modules, 'yaml', None)

    filename = str(tmpdir.join('test.fits'))
    t = Table([mixin_cols['sc']])
    with pytest.raises(TypeError) as err:
        t.write(filename)
    assert "cannot write type SkyCoord column 'col0' to FITS without PyYAML" in str(err.value)


@pytest.mark.skipif('not HAS_YAML')
def test_info_attributes_with_no_mixins(tmpdir):
    """Even if there are no mixin columns, if there is metadata that would be lost it still
    gets serialized
    """
    filename = str(tmpdir.join('test.fits'))
    t = Table([[1.0, 2.0]])
    t['col0'].description = 'hello' * 40
    t['col0'].format = '{:8.4f}'
    t['col0'].meta['a'] = {'b': 'c'}
    t.write(filename, overwrite=True)

    t2 = Table.read(filename)
    assert t2['col0'].description == 'hello' * 40
    assert t2['col0'].format == '{:8.4f}'
    assert t2['col0'].meta['a'] == {'b': 'c'}


@pytest.mark.skipif('not HAS_YAML')
@pytest.mark.parametrize('method', ['set_cols', 'names', 'class'])
def test_round_trip_masked_table_serialize_mask(tmpdir, method):
    """
    Same as previous test but set the serialize_method to 'data_mask' so mask is
    written out and the behavior is all correct.
    """
    filename = str(tmpdir.join('test.fits'))

    t = simple_table(masked=True)  # int, float, and str cols with one masked element

    # MaskedColumn but no masked elements.  See table the MaskedColumnInfo class
    # _represent_as_dict() method for info about we test a column with no masked elements.
    t['d'] = [1, 2, 3]

    if method == 'set_cols':
        for col in t.itercols():
            col.info.serialize_method['fits'] = 'data_mask'
        t.write(filename)
    elif method == 'names':
        t.write(filename, serialize_method={'a': 'data_mask', 'b': 'data_mask',
                                            'c': 'data_mask', 'd': 'data_mask'})
    elif method == 'class':
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
def test_read_serialized_without_yaml(tmpdir, monkeypatch):
    filename = str(tmpdir.join('test.fits'))
    t = Table([mixin_cols['sc']])
    t.write(filename)

    monkeypatch.setitem(sys.modules, 'yaml', None)
    with pytest.warns(AstropyUserWarning):
        t2 = Table.read(filename)

    assert t2.colnames == ['col0.ra', 'col0.dec']
    assert len(t2) == 2


@pytest.mark.skipif('not HAS_YAML')
def test_meta_not_modified(tmpdir):
    filename = str(tmpdir.join('test.fits'))
    t = Table(data=[Column([1, 2], 'a', description='spam')])
    t.meta['comments'] = ['a', 'b']
    assert len(t.meta) == 1
    t.write(filename)
    assert len(t.meta) == 1
    assert t.meta['comments'] == ['a', 'b']
