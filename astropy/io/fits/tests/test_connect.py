import os
import sys
import warnings

import numpy as np
from numpy.testing import assert_allclose

from ... import fits
from .. import HDUList, PrimaryHDU, BinTableHDU
from ....table import Table
from .... import units as u
from .... import log
from ....tests.helper import pytest, catch_warnings
from astropy.units.format.fits import UnitScaleError

DATA = os.path.join(os.path.dirname(__file__), 'data')


def equal_data(a, b):
    for name in a.dtype.names:
        if not np.all(a[name] == b[name]):
            return False
    return True


class TestSingleTable(object):

    def setup_class(self):
        self.data = np.array(list(zip([1, 2, 3, 4],
                                      ['a', 'b', 'c', 'd'],
                                      [2.3, 4.5, 6.7, 8.9])),
                             dtype=[(str('a'), int), (str('b'), str('U1')), (str('c'), float)])

    def test_simple(self, tmpdir):
        filename = str(tmpdir.join('test_simple.fits'))
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
        t1.meta['COMMENT'] = ['this', 'is', 'a', 'long', 'comment']
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

    def test_with_units(self, tmpdir):
        filename = str(tmpdir.join('test_with_units.fits'))
        t1 = Table(self.data)
        t1['a'].unit = u.m
        t1['c'].unit = u.km / u.s
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename)
        assert equal_data(t1, t2)
        assert t2['a'].unit == u.m
        assert t2['c'].unit == u.km / u.s

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
                                dtype=[(str('a'), np.float64), (str('b'), np.float32)])
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
        hdu.writeto(filename)
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


class TestMultipleHDU(object):

    def setup_class(self):
        self.data1 = np.array(list(zip([1, 2, 3, 4],
                                       ['a', 'b', 'c', 'd'],
                                       [2.3, 4.5, 6.7, 8.9])),
                              dtype=[(str('a'), int), (str('b'), str('U1')), (str('c'), float)])
        self.data2 = np.array(list(zip([1.4, 2.3, 3.2, 4.7],
                                       [2.3, 4.5, 6.7, 8.9])),
                              dtype=[(str('p'), float), (str('q'), float)])
        hdu1 = PrimaryHDU()
        hdu2 = BinTableHDU(self.data1, name='first')
        hdu3 = BinTableHDU(self.data2, name='second')

        self.hdus = HDUList([hdu1, hdu2, hdu3])

    def setup_method(self, method):
        warnings.filterwarnings('always')

    def test_read(self, tmpdir):
        filename = str(tmpdir.join('test_read.fits'))
        self.hdus.writeto(filename)
        with catch_warnings() as l:
            t = Table.read(filename)
        assert len(l) == 1
        assert str(l[0].message).startswith(
            'hdu= was not specified but multiple tables are present, reading in first available table (hdu=1)')
        assert equal_data(t, self.data1)

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

    def test_read_from_hdulist(self):
        with catch_warnings() as l:
            t = Table.read(self.hdus)
        assert len(l) == 1
        assert str(l[0].message).startswith(
            'hdu= was not specified but multiple tables are present, reading in first available table (hdu=1)')
        assert equal_data(t, self.data1)

    def test_read_from_hdulist_with_hdu_0(self, tmpdir):
        with pytest.raises(ValueError) as exc:
            Table.read(self.hdus, hdu=0)
        assert exc.value.args[0] == 'No table found in hdu=0'

    @pytest.mark.parametrize('hdu', [1, 'first'])
    def test_read_from_hdulist_with_hdu_1(self, tmpdir, hdu):
        with catch_warnings() as l:
            t = Table.read(self.hdus, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data1)

    @pytest.mark.parametrize('hdu', [2, 'second'])
    def test_read_from_hdulist_with_hdu_2(self, tmpdir, hdu):
        with catch_warnings() as l:
            t = Table.read(self.hdus, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data2)

    def test_read_from_single_hdu(self):
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
    assert np.all(t['c1'].data == np.array([1,2]))
    assert np.all(t['c2'].data == np.array(['abc', 'xy ']))
    assert_allclose(t['c3'].data, np.array([3.70000007153, 6.6999997139]))
    assert np.all(t['c4'].data == np.array([False, True]))


def test_scale_error():
    a = [1, 4, 5]
    b = [2.0, 5.0, 8.2]
    c = ['x', 'y', 'z']
    t = Table([a, b, c], names=('a', 'b', 'c'), meta={'name': 'first table'})
    t['a'].unit='1.2'
    with pytest.raises(UnitScaleError) as exc:
        t.write('t.fits',format='fits', overwrite=True)
    assert exc.value.args[0]=="The column 'a' could not be stored in FITS format because it has a scale '(1.2)' that is not recognized by the FITS standard. Either scale the data or change the units."


def test_bool_column(tmpdir):
    """
    Regression test for https://github.com/astropy/astropy/issues/1953

    Ensures that Table columns of bools are properly written to a FITS table.
    """

    arr = np.ones(5, dtype=bool)
    arr[::2] == False

    t = Table([arr])
    t.write(str(tmpdir.join('test.fits')), overwrite=True)

    with fits.open(str(tmpdir.join('test.fits'))) as hdul:
        assert hdul[1].data['col0'].dtype == np.dtype('bool')
        assert np.all(hdul[1].data['col0'] == arr)
