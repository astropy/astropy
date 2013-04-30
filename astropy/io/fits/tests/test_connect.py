import pytest
import numpy as np
import warnings

from astropy.io.fits import HDUList, PrimaryHDU, BinTableHDU
from astropy.table import Table
from astropy import units as u
from astropy import log


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
                             dtype=[('a', int), ('b', 'U1'), ('c', float)])

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
        t1['a'].units = u.m
        t1['c'].units = u.km / u.s
        t1.write(filename, overwrite=True)
        t2 = Table.read(filename)
        assert equal_data(t1, t2)
        assert t2['a'].units == u.m
        assert t2['c'].units == u.km / u.s

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
        assert np.all(t1['b'].mask == t2['b'].mask)
        assert np.all(t1['c'].mask == t2['c'].mask)


class TestMultipleHDU(object):

    def setup_class(self):
        self.data1 = np.array(list(zip([1, 2, 3, 4],
                                       ['a', 'b', 'c', 'd'],
                                       [2.3, 4.5, 6.7, 8.9])),
                              dtype=[('a', int), ('b', 'U1'), ('c', float)])
        self.data2 = np.array(list(zip([1.4, 2.3, 3.2, 4.7],
                                       [2.3, 4.5, 6.7, 8.9])),
                              dtype=[('p', float), ('q', float)])
        hdu1 = PrimaryHDU()
        hdu2 = BinTableHDU(self.data1, name='first')
        hdu3 = BinTableHDU(self.data2, name='second')

        self.hdus = HDUList([hdu1, hdu2, hdu3])

    def setup_method(self, method):
        warnings.filterwarnings('always')

    def test_read(self, tmpdir):
        filename = str(tmpdir.join('test_read.fits'))
        self.hdus.writeto(filename)
        with log.log_to_list() as l:
            t = Table.read(filename)
        assert len(l) == 1
        assert l[0].message == 'hdu= was not specified but multiple tables are present, reading in first available table (hdu=1)'
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
        with log.log_to_list() as l:
            t = Table.read(filename, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data1)

    @pytest.mark.parametrize('hdu', [2, 'second'])
    def test_read_with_hdu_2(self, tmpdir, hdu):
        filename = str(tmpdir.join('test_read_with_hdu_2.fits'))
        self.hdus.writeto(filename)
        with log.log_to_list() as l:
            t = Table.read(filename, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data2)

    def test_read_from_hdulist(self):
        with log.log_to_list() as l:
            t = Table.read(self.hdus)
        assert len(l) == 1
        assert l[0].message == 'hdu= was not specified but multiple tables are present, reading in first available table (hdu=1)'
        assert equal_data(t, self.data1)

    def test_read_from_hdulist_with_hdu_0(self, tmpdir):
        with pytest.raises(ValueError) as exc:
            Table.read(self.hdus, hdu=0)
        assert exc.value.args[0] == 'No table found in hdu=0'

    @pytest.mark.parametrize('hdu', [1, 'first'])
    def test_read_from_hdulist_with_hdu_1(self, tmpdir, hdu):
        with log.log_to_list() as l:
            t = Table.read(self.hdus, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data1)

    @pytest.mark.parametrize('hdu', [2, 'second'])
    def test_read_from_hdulist_with_hdu_2(self, tmpdir, hdu):
        with log.log_to_list() as l:
            t = Table.read(self.hdus, hdu=hdu)
        assert len(l) == 0
        assert equal_data(t, self.data2)

    def test_read_from_single_hdu(self):
        with log.log_to_list() as l:
            t = Table.read(self.hdus[1])
        assert len(l) == 0
        assert equal_data(t, self.data1)
