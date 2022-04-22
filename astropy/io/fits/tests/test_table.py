# Licensed under a 3-clause BSD style license - see PYFITS.rst

import contextlib
import copy
import gc
import pickle
import re
import sys
import warnings

import pytest
import numpy as np
from numpy import char as chararray

try:
    import objgraph
    HAVE_OBJGRAPH = True
except ImportError:
    HAVE_OBJGRAPH = False

from astropy.io import fits
from astropy.table import Table
from astropy.units import UnitsWarning, Unit, UnrecognizedUnit
from astropy.utils.compat import NUMPY_LT_1_22, NUMPY_LT_1_22_1
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning

from astropy.io.fits.column import ColumnAttribute, Delayed, NUMPY2FITS
from astropy.io.fits.util import decode_ascii
from astropy.io.fits.verify import VerifyError
from . import FitsTestCase


def comparefloats(a, b):
    """
    Compare two float scalars or arrays and see if they are consistent

    Consistency is determined ensuring the difference is less than the
    expected amount. Return True if consistent, False if any differences.
    """

    aa = a
    bb = b
    # compute expected precision
    if aa.dtype.name == 'float32' or bb.dtype.name == 'float32':
        precision = 0.000001
    else:
        precision = 0.0000000000000001
    precision = 0.00001  # until precision problem is fixed in astropy.io.fits
    diff = np.absolute(aa - bb)
    mask0 = aa == 0
    masknz = aa != 0.
    if np.any(mask0):
        if diff[mask0].max() != 0.:
            return False
    if np.any(masknz):
        if (diff[masknz] / np.absolute(aa[masknz])).max() > precision:
            return False
    return True


def comparerecords(a, b):
    """
    Compare two record arrays

    Does this field by field, using approximation testing for float columns
    (Complex not yet handled.)
    Column names not compared, but column types and sizes are.
    """

    nfieldsa = len(a.dtype.names)
    nfieldsb = len(b.dtype.names)
    if nfieldsa != nfieldsb:
        print("number of fields don't match")
        return False
    for i in range(nfieldsa):
        fielda = a.field(i)
        fieldb = b.field(i)
        if fielda.dtype.char == 'S':
            fielda = decode_ascii(fielda)
        if fieldb.dtype.char == 'S':
            fieldb = decode_ascii(fieldb)
        if (not isinstance(fielda, type(fieldb)) and not
                isinstance(fieldb, type(fielda))):
            print("type(fielda): ", type(fielda), " fielda: ", fielda)
            print("type(fieldb): ", type(fieldb), " fieldb: ", fieldb)
            print(f'field {i} type differs')
            return False
        if len(fielda) and isinstance(fielda[0], np.floating):
            if not comparefloats(fielda, fieldb):
                print("fielda: ", fielda)
                print("fieldb: ", fieldb)
                print(f'field {i} differs')
                return False
        elif (isinstance(fielda, fits.column._VLF) or
              isinstance(fieldb, fits.column._VLF)):
            for row in range(len(fielda)):
                if np.any(fielda[row] != fieldb[row]):
                    print(f'fielda[{row}]: {fielda[row]}')
                    print(f'fieldb[{row}]: {fieldb[row]}')
                    print(f'field {i} differs in row {row}')
        else:
            if np.any(fielda != fieldb):
                print("fielda: ", fielda)
                print("fieldb: ", fieldb)
                print(f'field {i} differs')
                return False
    return True


def _assert_attr_col(new_tbhdu, tbhdu):
    """
    Helper function to compare column attributes
    """
    # Double check that the headers are equivalent
    assert tbhdu.columns.names == new_tbhdu.columns.names
    attrs = [k for k, v in fits.Column.__dict__.items()
             if isinstance(v, ColumnAttribute)]
    for name in tbhdu.columns.names:
        col = tbhdu.columns[name]
        new_col = new_tbhdu.columns[name]
        for attr in attrs:
            if getattr(col, attr) and getattr(new_col, attr):
                assert getattr(col, attr) == getattr(new_col, attr)


class TestTableFunctions(FitsTestCase):
    def test_constructor_copies_header(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/153

        Ensure that a header from one HDU is copied when used to initialize new
        HDU.

        This is like the test of the same name in test_image, but tests this
        for tables as well.
        """

        ifd = fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU()])
        thdr = ifd[1].header
        thdr['FILENAME'] = 'labq01i3q_rawtag.fits'

        thdu = fits.BinTableHDU(header=thdr)
        ofd = fits.HDUList(thdu)
        ofd[0].header['FILENAME'] = 'labq01i3q_flt.fits'

        # Original header should be unchanged
        assert thdr['FILENAME'] == 'labq01i3q_rawtag.fits'

    def test_open(self):
        # open some existing FITS files:
        tt = fits.open(self.data('tb.fits'))
        fd = fits.open(self.data('test0.fits'))

        # create some local arrays
        a1 = chararray.array(['abc', 'def', 'xx'])
        r1 = np.array([11., 12., 13.], dtype=np.float32)

        # create a table from scratch, using a mixture of columns from existing
        # tables and locally created arrays:

        # first, create individual column definitions

        c1 = fits.Column(name='abc', format='3A', array=a1)
        c2 = fits.Column(name='def', format='E', array=r1)
        a3 = np.array([3, 4, 5], dtype='i2')
        c3 = fits.Column(name='xyz', format='I', array=a3)
        a4 = np.array([1, 2, 3], dtype='i2')
        c4 = fits.Column(name='t1', format='I', array=a4)
        a5 = np.array([3 + 3j, 4 + 4j, 5 + 5j], dtype='c8')
        c5 = fits.Column(name='t2', format='C', array=a5)

        # Note that X format must be two-D array
        a6 = np.array([[0], [1], [0]], dtype=np.uint8)
        c6 = fits.Column(name='t3', format='X', array=a6)
        a7 = np.array([101, 102, 103], dtype='i4')
        c7 = fits.Column(name='t4', format='J', array=a7)
        a8 = np.array([[1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1],
                       [0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0],
                       [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1]], dtype=np.uint8)
        c8 = fits.Column(name='t5', format='11X', array=a8)

        # second, create a column-definitions object for all columns in a table

        x = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8])

        tbhdu = fits.BinTableHDU.from_columns(x)

        # another way to create a table is by using existing table's
        # information:

        x2 = fits.ColDefs(tt[1])
        t2 = fits.BinTableHDU.from_columns(x2, nrows=2)
        ra = np.rec.array([
            (1, 'abc', 3.7000002861022949, 0),
            (2, 'xy ', 6.6999998092651367, 1)], names='c1, c2, c3, c4')

        assert comparerecords(t2.data, ra)

        # the table HDU's data is a subclass of a record array, so we can
        # access one row like this:

        assert tbhdu.data[1][0] == a1[1]
        assert tbhdu.data[1][1] == r1[1]
        assert tbhdu.data[1][2] == a3[1]
        assert tbhdu.data[1][3] == a4[1]
        assert tbhdu.data[1][4] == a5[1]
        assert (tbhdu.data[1][5] == a6[1].view('bool')).all()
        assert tbhdu.data[1][6] == a7[1]
        assert (tbhdu.data[1][7] == a8[1]).all()

        # and a column like this:
        assert str(tbhdu.data.field('abc')) == "['abc' 'def' 'xx']"

        # An alternative way to create a column-definitions object is from an
        # existing table.
        _ = fits.ColDefs(tt[1])

        # now we write out the newly created table HDU to a FITS file:
        fout = fits.HDUList(fits.PrimaryHDU())
        fout.append(tbhdu)
        fout.writeto(self.temp('tableout1.fits'), overwrite=True)

        with fits.open(self.temp('tableout1.fits')) as f2:
            temp = f2[1].data.field(7)
            assert (temp[0] == [True, True, False, True, False, True,
                                True, True, False, False, True]).all()

        # An alternative way to create an output table FITS file:
        fout2 = fits.open(self.temp('tableout2.fits'), 'append')
        fout2.append(fd[0])
        fout2.append(tbhdu)
        fout2.close()
        tt.close()
        fd.close()

    def test_binary_table(self):
        # binary table:
        t = fits.open(self.data('tb.fits'))
        assert t[1].header['tform1'] == '1J'

        info = {'name': ['c1', 'c2', 'c3', 'c4'],
                'format': ['1J', '3A', '1E', '1L'],
                'unit': ['', '', '', ''],
                'null': [-2147483647, '', '', ''],
                'bscale': ['', '', 3, ''],
                'bzero': ['', '', 0.4, ''],
                'disp': ['I11', 'A3', 'G15.7', 'L6'],
                'start': ['', '', '', ''],
                'dim': ['', '', '', ''],
                'coord_inc': ['', '', '', ''],
                'coord_type': ['', '', '', ''],
                'coord_unit': ['', '', '', ''],
                'coord_ref_point': ['', '', '', ''],
                'coord_ref_value': ['', '', '', ''],
                'time_ref_pos': ['', '', '', '']}

        assert t[1].columns.info(output=False) == info

        ra = np.rec.array([
            (1, 'abc', 3.7000002861022949, 0),
            (2, 'xy ', 6.6999998092651367, 1)], names='c1, c2, c3, c4')

        assert comparerecords(t[1].data, ra[:2])

        # Change scaled field and scale back to the original array
        t[1].data.field('c4')[0] = 1
        t[1].data._scale_back()
        assert str(np.rec.recarray.field(t[1].data, 'c4')) == '[84 84]'

        # look at data column-wise
        assert (t[1].data.field(0) == np.array([1, 2])).all()

        # When there are scaled columns, the raw data are in data._parent

        t.close()

    def test_ascii_table(self):
        # ASCII table
        a = fits.open(self.data('ascii.fits'))
        ra1 = np.rec.array([
            (10.123000144958496, 37),
            (5.1999998092651367, 23),
            (15.609999656677246, 17),
            (0.0, 0),
            (345.0, 345)], names='c1, c2')
        assert comparerecords(a[1].data, ra1)

        # Test slicing
        a2 = a[1].data[2:][2:]
        ra2 = np.rec.array([(345.0, 345)], names='c1, c2')

        assert comparerecords(a2, ra2)

        assert (a2.field(1) == np.array([345])).all()

        ra3 = np.rec.array([
            (10.123000144958496, 37),
            (15.609999656677246, 17),
            (345.0, 345)
        ], names='c1, c2')

        assert comparerecords(a[1].data[::2], ra3)

        # Test Start Column

        a1 = chararray.array(['abcd', 'def'])
        r1 = np.array([11., 12.])
        c1 = fits.Column(name='abc', format='A3', start=19, array=a1)
        c2 = fits.Column(name='def', format='E', start=3, array=r1)
        c3 = fits.Column(name='t1', format='I', array=[91, 92, 93])
        hdu = fits.TableHDU.from_columns([c2, c1, c3])

        assert (dict(hdu.data.dtype.fields) ==
                {'abc': (np.dtype('|S3'), 18),
                 'def': (np.dtype('|S15'), 2),
                 't1': (np.dtype('|S10'), 21)})
        hdu.writeto(self.temp('toto.fits'), overwrite=True)
        hdul = fits.open(self.temp('toto.fits'))
        assert comparerecords(hdu.data, hdul[1].data)
        hdul.close()

        # Test Scaling

        r1 = np.array([11., 12.])
        c2 = fits.Column(name='def', format='D', array=r1, bscale=2.3,
                         bzero=0.6)
        hdu = fits.TableHDU.from_columns([c2])
        hdu.writeto(self.temp('toto.fits'), overwrite=True)
        with open(self.temp('toto.fits')) as f:
            assert '4.95652173913043548D+00' in f.read()
        with fits.open(self.temp('toto.fits')) as hdul:
            assert comparerecords(hdu.data, hdul[1].data)

        # Test Integer precision according to width

        c1 = fits.Column(name='t2', format='I2', array=[91, 92, 93])
        c2 = fits.Column(name='t4', format='I5', array=[91, 92, 93])
        c3 = fits.Column(name='t8', format='I10', array=[91, 92, 93])
        hdu = fits.TableHDU.from_columns([c1, c2, c3])

        assert c1.array.dtype == np.int16
        assert c2.array.dtype == np.int32
        assert c3.array.dtype == np.int64
        hdu.writeto(self.temp('toto.fits'), overwrite=True)
        with fits.open(self.temp('toto.fits')) as hdul:
            assert comparerecords(hdu.data, hdul[1].data)

        a.close()

    def test_endianness(self):
        x = np.ndarray((1,), dtype=object)
        channelsIn = np.array([3], dtype='uint8')
        x[0] = channelsIn
        col = fits.Column(name="Channels", format="PB()", array=x)
        cols = fits.ColDefs([col])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.name = "RFI"
        tbhdu.writeto(self.temp('testendian.fits'), overwrite=True)
        hduL = fits.open(self.temp('testendian.fits'))
        rfiHDU = hduL['RFI']
        data = rfiHDU.data
        channelsOut = data.field('Channels')[0]
        assert (channelsIn == channelsOut).all()
        hduL.close()

    def test_column_endianness(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/77
        (Astropy doesn't preserve byte order of non-native order column arrays)
        """

        a = [1., 2., 3., 4.]
        a1 = np.array(a, dtype='<f8')
        a2 = np.array(a, dtype='>f8')

        col1 = fits.Column(name='a', format='D', array=a1)
        col2 = fits.Column(name='b', format='D', array=a2)
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.BinTableHDU.from_columns(cols)

        assert (tbhdu.data['a'] == a1).all()
        assert (tbhdu.data['b'] == a2).all()

        # Double check that the array is converted to the correct byte-order
        # for FITS (big-endian).
        tbhdu.writeto(self.temp('testendian.fits'), overwrite=True)
        with fits.open(self.temp('testendian.fits')) as hdul:
            assert (hdul[1].data['a'] == a2).all()
            assert (hdul[1].data['b'] == a2).all()

    def test_recarray_to_bintablehdu(self):
        bright = np.rec.array(
            [(1, 'Serius', -1.45, 'A1V'),
             (2, 'Canopys', -0.73, 'F0Ib'),
             (3, 'Rigil Kent', -0.1, 'G2V')],
            formats='int16,a20,float32,a10',
            names='order,name,mag,Sp')
        hdu = fits.BinTableHDU(bright)
        assert comparerecords(hdu.data, bright)
        hdu.writeto(self.temp('toto.fits'), overwrite=True)
        hdul = fits.open(self.temp('toto.fits'))
        assert comparerecords(hdu.data, hdul[1].data)
        assert comparerecords(bright, hdul[1].data)
        hdul.close()

    def test_numpy_ndarray_to_bintablehdu(self):
        desc = np.dtype({'names': ['order', 'name', 'mag', 'Sp'],
                         'formats': ['int', 'S20', 'float32', 'S10']})
        a = np.array([(1, 'Serius', -1.45, 'A1V'),
                      (2, 'Canopys', -0.73, 'F0Ib'),
                      (3, 'Rigil Kent', -0.1, 'G2V')], dtype=desc)
        hdu = fits.BinTableHDU(a)
        assert comparerecords(hdu.data, a.view(fits.FITS_rec))
        hdu.writeto(self.temp('toto.fits'), overwrite=True)
        hdul = fits.open(self.temp('toto.fits'))
        assert comparerecords(hdu.data, hdul[1].data)
        hdul.close()

    def test_numpy_ndarray_to_bintablehdu_with_unicode(self):
        desc = np.dtype({'names': ['order', 'name', 'mag', 'Sp'],
                         'formats': ['int', 'U20', 'float32', 'U10']})
        a = np.array([(1, 'Serius', -1.45, 'A1V'),
                      (2, 'Canopys', -0.73, 'F0Ib'),
                      (3, 'Rigil Kent', -0.1, 'G2V')], dtype=desc)
        hdu = fits.BinTableHDU(a)
        assert comparerecords(hdu.data, a.view(fits.FITS_rec))
        hdu.writeto(self.temp('toto.fits'), overwrite=True)
        hdul = fits.open(self.temp('toto.fits'))
        assert comparerecords(hdu.data, hdul[1].data)
        hdul.close()

    def test_new_table_from_recarray(self):
        bright = np.rec.array([(1, 'Serius', -1.45, 'A1V'),
                               (2, 'Canopys', -0.73, 'F0Ib'),
                               (3, 'Rigil Kent', -0.1, 'G2V')],
                              formats='int16,a20,float64,a10',
                              names='order,name,mag,Sp')
        hdu = fits.TableHDU.from_columns(bright, nrows=2)

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert (id(hdu.data._coldefs.columns[0].array) ==
                id(hdu.data._coldefs._arrays[0]))
        assert (id(hdu.data._coldefs.columns[0].array) ==
                id(hdu.columns.columns[0].array))
        assert (id(hdu.data._coldefs.columns[0].array) ==
                id(hdu.columns._arrays[0]))

        # Ensure I can change the value of one data element and it effects
        # all of the others.
        hdu.data[0][0] = 213

        assert hdu.data[0][0] == 213
        assert hdu.data._coldefs._arrays[0][0] == 213
        assert hdu.data._coldefs.columns[0].array[0] == 213
        assert hdu.columns._arrays[0][0] == 213
        assert hdu.columns.columns[0].array[0] == 213

        hdu.data._coldefs._arrays[0][0] = 100

        assert hdu.data[0][0] == 100
        assert hdu.data._coldefs._arrays[0][0] == 100
        assert hdu.data._coldefs.columns[0].array[0] == 100
        assert hdu.columns._arrays[0][0] == 100
        assert hdu.columns.columns[0].array[0] == 100

        hdu.data._coldefs.columns[0].array[0] = 500
        assert hdu.data[0][0] == 500
        assert hdu.data._coldefs._arrays[0][0] == 500
        assert hdu.data._coldefs.columns[0].array[0] == 500
        assert hdu.columns._arrays[0][0] == 500
        assert hdu.columns.columns[0].array[0] == 500

        hdu.columns._arrays[0][0] = 600
        assert hdu.data[0][0] == 600
        assert hdu.data._coldefs._arrays[0][0] == 600
        assert hdu.data._coldefs.columns[0].array[0] == 600
        assert hdu.columns._arrays[0][0] == 600
        assert hdu.columns.columns[0].array[0] == 600

        hdu.columns.columns[0].array[0] = 800
        assert hdu.data[0][0] == 800
        assert hdu.data._coldefs._arrays[0][0] == 800
        assert hdu.data._coldefs.columns[0].array[0] == 800
        assert hdu.columns._arrays[0][0] == 800
        assert hdu.columns.columns[0].array[0] == 800

        assert (hdu.data.field(0) ==
                np.array([800, 2], dtype=np.int16)).all()
        assert hdu.data[0][1] == 'Serius'
        assert hdu.data[1][1] == 'Canopys'
        assert (hdu.data.field(2) ==
                np.array([-1.45, -0.73], dtype=np.float64)).all()
        assert hdu.data[0][3] == 'A1V'
        assert hdu.data[1][3] == 'F0Ib'

        hdu.writeto(self.temp('toto.fits'), overwrite=True)

        with fits.open(self.temp('toto.fits')) as hdul:
            assert (hdul[1].data.field(0) ==
                    np.array([800, 2], dtype=np.int16)).all()
            assert hdul[1].data[0][1] == 'Serius'
            assert hdul[1].data[1][1] == 'Canopys'
            assert (hdul[1].data.field(2) ==
                    np.array([-1.45, -0.73], dtype=np.float64)).all()
            assert hdul[1].data[0][3] == 'A1V'
            assert hdul[1].data[1][3] == 'F0Ib'
        del hdul

        hdu = fits.BinTableHDU.from_columns(bright, nrows=2)
        tmp = np.rec.array([(1, 'Serius', -1.45, 'A1V'),
                            (2, 'Canopys', -0.73, 'F0Ib')],
                           formats='int16,a20,float64,a10',
                           names='order,name,mag,Sp')
        assert comparerecords(hdu.data, tmp)
        hdu.writeto(self.temp('toto.fits'), overwrite=True)
        with fits.open(self.temp('toto.fits')) as hdul:
            assert comparerecords(hdu.data, hdul[1].data)

    def test_new_fitsrec(self):
        """
        Tests creating a new FITS_rec object from a multi-field ndarray.
        """

        with fits.open(self.data('tb.fits')) as h:
            data = h[1].data
        new_data = np.array([(3, 'qwe', 4.5, False)], dtype=data.dtype)
        appended = np.append(data, new_data).view(fits.FITS_rec)
        assert repr(appended).startswith('FITS_rec(')
        # This test used to check the entire string representation of FITS_rec,
        # but that has problems between different numpy versions.  Instead just
        # check that the FITS_rec was created, and we'll let subsequent tests
        # worry about checking values and such

    def test_appending_a_column(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        counts = np.array([412, 434, 408, 417])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[0, 1, 0, 0])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.writeto(self.temp('table2.fits'))

        # Append the rows of table 2 after the rows of table 1
        # The column definitions are assumed to be the same

        # Open the two files we want to append
        t1 = fits.open(self.temp('table1.fits'))
        t2 = fits.open(self.temp('table2.fits'))

        # Get the number of rows in the table from the first file
        nrows1 = t1[1].data.shape[0]

        # Get the total number of rows in the resulting appended table
        nrows = t1[1].data.shape[0] + t2[1].data.shape[0]

        assert (t1[1].columns._arrays[1] is t1[1].columns.columns[1].array)

        # Create a new table that consists of the data from the first table
        # but has enough space in the ndarray to hold the data from both tables
        hdu = fits.BinTableHDU.from_columns(t1[1].columns, nrows=nrows)

        # For each column in the tables append the data from table 2 after the
        # data from table 1.
        for i in range(len(t1[1].columns)):
            hdu.data.field(i)[nrows1:] = t2[1].data.field(i)

        hdu.writeto(self.temp('newtable.fits'))

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 4, (), '', ''),
                (1, '', 1, 'BinTableHDU', 19, '8R x 5C', '[10A, J, 10A, 5E, L]',
                 '')]

        assert fits.info(self.temp('newtable.fits'), output=False) == info

        z = np.array([0., 0., 0., 0., 0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '', z, True),
             ('NGC2', 334, '', z, False),
             ('NGC3', 308, '', z, True),
             ('NCG4', 317, '', z, True),
             ('NGC5', 412, '', z, False),
             ('NGC6', 434, '', z, True),
             ('NGC7', 408, '', z, False),
             ('NCG8', 417, '', z, False)],
            formats='a10,u4,a10,5f4,l')

        assert comparerecords(hdu.data, array)

        # Verify that all of the references to the data point to the same
        # numarray
        hdu.data[0][1] = 300
        assert hdu.data._coldefs._arrays[1][0] == 300
        assert hdu.data._coldefs.columns[1].array[0] == 300
        assert hdu.columns._arrays[1][0] == 300
        assert hdu.columns.columns[1].array[0] == 300
        assert hdu.data[0][1] == 300

        hdu.data._coldefs._arrays[1][0] = 200
        assert hdu.data._coldefs._arrays[1][0] == 200
        assert hdu.data._coldefs.columns[1].array[0] == 200
        assert hdu.columns._arrays[1][0] == 200
        assert hdu.columns.columns[1].array[0] == 200
        assert hdu.data[0][1] == 200

        hdu.data._coldefs.columns[1].array[0] = 100
        assert hdu.data._coldefs._arrays[1][0] == 100
        assert hdu.data._coldefs.columns[1].array[0] == 100
        assert hdu.columns._arrays[1][0] == 100
        assert hdu.columns.columns[1].array[0] == 100
        assert hdu.data[0][1] == 100

        hdu.columns._arrays[1][0] = 90
        assert hdu.data._coldefs._arrays[1][0] == 90
        assert hdu.data._coldefs.columns[1].array[0] == 90
        assert hdu.columns._arrays[1][0] == 90
        assert hdu.columns.columns[1].array[0] == 90
        assert hdu.data[0][1] == 90

        hdu.columns.columns[1].array[0] = 80
        assert hdu.data._coldefs._arrays[1][0] == 80
        assert hdu.data._coldefs.columns[1].array[0] == 80
        assert hdu.columns._arrays[1][0] == 80
        assert hdu.columns.columns[1].array[0] == 80
        assert hdu.data[0][1] == 80

        # Same verification from the file
        hdul = fits.open(self.temp('newtable.fits'))
        hdu = hdul[1]
        hdu.data[0][1] = 300
        assert hdu.data._coldefs._arrays[1][0] == 300
        assert hdu.data._coldefs.columns[1].array[0] == 300
        assert hdu.columns._arrays[1][0] == 300
        assert hdu.columns.columns[1].array[0] == 300
        assert hdu.data[0][1] == 300

        hdu.data._coldefs._arrays[1][0] = 200
        assert hdu.data._coldefs._arrays[1][0] == 200
        assert hdu.data._coldefs.columns[1].array[0] == 200
        assert hdu.columns._arrays[1][0] == 200
        assert hdu.columns.columns[1].array[0] == 200
        assert hdu.data[0][1] == 200

        hdu.data._coldefs.columns[1].array[0] = 100
        assert hdu.data._coldefs._arrays[1][0] == 100
        assert hdu.data._coldefs.columns[1].array[0] == 100
        assert hdu.columns._arrays[1][0] == 100
        assert hdu.columns.columns[1].array[0] == 100
        assert hdu.data[0][1] == 100

        hdu.columns._arrays[1][0] = 90
        assert hdu.data._coldefs._arrays[1][0] == 90
        assert hdu.data._coldefs.columns[1].array[0] == 90
        assert hdu.columns._arrays[1][0] == 90
        assert hdu.columns.columns[1].array[0] == 90
        assert hdu.data[0][1] == 90

        hdu.columns.columns[1].array[0] = 80
        assert hdu.data._coldefs._arrays[1][0] == 80
        assert hdu.data._coldefs.columns[1].array[0] == 80
        assert hdu.columns._arrays[1][0] == 80
        assert hdu.columns.columns[1].array[0] == 80
        assert hdu.data[0][1] == 80

        t1.close()
        t2.close()
        hdul.close()

    def test_adding_a_column(self):
        # Tests adding a column to a table.
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        assert tbhdu.columns.names == ['target', 'counts', 'notes', 'spectrum']
        coldefs1 = coldefs + c5

        tbhdu1 = fits.BinTableHDU.from_columns(coldefs1)
        assert tbhdu1.columns.names == ['target', 'counts', 'notes',
                                        'spectrum', 'flag']

        z = np.array([0., 0., 0., 0., 0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '', z, True),
             ('NGC2', 334, '', z, False),
             ('NGC3', 308, '', z, True),
             ('NCG4', 317, '', z, True)],
            formats='a10,u4,a10,5f4,l')
        assert comparerecords(tbhdu1.data, array)

    def test_adding_a_column_inplace(self):
        # Tests adding a column to a table.
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        assert tbhdu.columns.names == ['target', 'counts', 'notes', 'spectrum']

        tbhdu.columns.add_col(c5)
        assert tbhdu.columns.names == ['target', 'counts', 'notes',
                                       'spectrum', 'flag']

        z = np.array([0., 0., 0., 0., 0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '', z, True),
             ('NGC2', 334, '', z, False),
             ('NGC3', 308, '', z, True),
             ('NCG4', 317, '', z, True)],
            formats='a10,u4,a10,5f4,l')
        assert comparerecords(tbhdu.data, array)

    def test_adding_a_column_to_file(self):
        hdul = fits.open(self.data('table.fits'))
        tbhdu = hdul[1]
        col = fits.Column(name='a', array=np.array([1, 2]), format='K')
        tbhdu.columns.add_col(col)
        assert tbhdu.columns.names == ['target', 'V_mag', 'a']
        array = np.rec.array(
            [('NGC1001', 11.1, 1),
             ('NGC1002', 12.3, 2),
             ('NGC1003', 15.2, 0)],
            formats='a20,f4,i8')
        assert comparerecords(tbhdu.data, array)
        hdul.close()

    def test_removing_a_column_inplace(self):
        # Tests adding a column to a table.
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        assert tbhdu.columns.names == ['target', 'counts', 'notes',
                                       'spectrum', 'flag']

        tbhdu.columns.del_col('flag')

        assert tbhdu.columns.names == ['target', 'counts', 'notes', 'spectrum']
        z = np.array([0., 0., 0., 0., 0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '', z),
             ('NGC2', 334, '', z),
             ('NGC3', 308, '', z),
             ('NCG4', 317, '', z)],
            formats='a10,u4,a10,5f4')
        assert comparerecords(tbhdu.data, array)

        tbhdu.columns.del_col('counts')
        tbhdu.columns.del_col('notes')

        assert tbhdu.columns.names == ['target', 'spectrum']
        array = np.rec.array(
            [('NGC1', z),
             ('NGC2', z),
             ('NGC3', z),
             ('NCG4', z)],
            formats='a10,5f4')
        assert comparerecords(tbhdu.data, array)

    def test_removing_a_column_from_file(self):
        hdul = fits.open(self.data('table.fits'))
        tbhdu = hdul[1]
        tbhdu.columns.del_col('V_mag')
        assert tbhdu.columns.names == ['target']
        array = np.rec.array(
            [('NGC1001', ),
             ('NGC1002', ),
             ('NGC1003', )],
            formats='a20')
        assert comparerecords(tbhdu.data, array)
        hdul.close()

    def test_merge_tables(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        counts = np.array([412, 434, 408, 417])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = fits.Column(name='target1', format='10A', array=names)
        c2 = fits.Column(name='counts1', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes1', format='A10')
        c4 = fits.Column(name='spectrum1', format='5E')
        c5 = fits.Column(name='flag1', format='L', array=[0, 1, 0, 0])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.writeto(self.temp('table2.fits'))

        # Merge the columns of table 2 after the columns of table 1
        # The column names are assumed to be different

        # Open the two files we want to append
        t1 = fits.open(self.temp('table1.fits'))
        t2 = fits.open(self.temp('table2.fits'))

        hdu = fits.BinTableHDU.from_columns(t1[1].columns + t2[1].columns)

        z = np.array([0., 0., 0., 0., 0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '', z, True, 'NGC5', 412, '', z, False),
             ('NGC2', 334, '', z, False, 'NGC6', 434, '', z, True),
             ('NGC3', 308, '', z, True, 'NGC7', 408, '', z, False),
             ('NCG4', 317, '', z, True, 'NCG8', 417, '', z, False)],
            formats='a10,u4,a10,5f4,l,a10,u4,a10,5f4,l')
        assert comparerecords(hdu.data, array)

        hdu.writeto(self.temp('newtable.fits'))

        # Verify that all of the references to the data point to the same
        # numarray
        hdu.data[0][1] = 300
        assert hdu.data._coldefs._arrays[1][0] == 300
        assert hdu.data._coldefs.columns[1].array[0] == 300
        assert hdu.columns._arrays[1][0] == 300
        assert hdu.columns.columns[1].array[0] == 300
        assert hdu.data[0][1] == 300

        hdu.data._coldefs._arrays[1][0] = 200
        assert hdu.data._coldefs._arrays[1][0] == 200
        assert hdu.data._coldefs.columns[1].array[0] == 200
        assert hdu.columns._arrays[1][0] == 200
        assert hdu.columns.columns[1].array[0] == 200
        assert hdu.data[0][1] == 200

        hdu.data._coldefs.columns[1].array[0] = 100
        assert hdu.data._coldefs._arrays[1][0] == 100
        assert hdu.data._coldefs.columns[1].array[0] == 100
        assert hdu.columns._arrays[1][0] == 100
        assert hdu.columns.columns[1].array[0] == 100
        assert hdu.data[0][1] == 100

        hdu.columns._arrays[1][0] = 90
        assert hdu.data._coldefs._arrays[1][0] == 90
        assert hdu.data._coldefs.columns[1].array[0] == 90
        assert hdu.columns._arrays[1][0] == 90
        assert hdu.columns.columns[1].array[0] == 90
        assert hdu.data[0][1] == 90

        hdu.columns.columns[1].array[0] = 80
        assert hdu.data._coldefs._arrays[1][0] == 80
        assert hdu.data._coldefs.columns[1].array[0] == 80
        assert hdu.columns._arrays[1][0] == 80
        assert hdu.columns.columns[1].array[0] == 80
        assert hdu.data[0][1] == 80

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 4, (), '', ''),
                (1, '', 1, 'BinTableHDU', 30, '4R x 10C',
                 '[10A, J, 10A, 5E, L, 10A, J, 10A, 5E, L]', '')]

        assert fits.info(self.temp('newtable.fits'), output=False) == info

        hdul = fits.open(self.temp('newtable.fits'))
        hdu = hdul[1]

        assert (hdu.columns.names ==
                ['target', 'counts', 'notes', 'spectrum', 'flag', 'target1',
                 'counts1', 'notes1', 'spectrum1', 'flag1'])

        z = np.array([0., 0., 0., 0., 0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '', z, True, 'NGC5', 412, '', z, False),
             ('NGC2', 334, '', z, False, 'NGC6', 434, '', z, True),
             ('NGC3', 308, '', z, True, 'NGC7', 408, '', z, False),
             ('NCG4', 317, '', z, True, 'NCG8', 417, '', z, False)],
            formats='a10,u4,a10,5f4,l,a10,u4,a10,5f4,l')
        assert comparerecords(hdu.data, array)

        # Same verification from the file
        hdu.data[0][1] = 300
        assert hdu.data._coldefs._arrays[1][0] == 300
        assert hdu.data._coldefs.columns[1].array[0] == 300
        assert hdu.columns._arrays[1][0] == 300
        assert hdu.columns.columns[1].array[0] == 300
        assert hdu.data[0][1] == 300

        hdu.data._coldefs._arrays[1][0] = 200
        assert hdu.data._coldefs._arrays[1][0] == 200
        assert hdu.data._coldefs.columns[1].array[0] == 200
        assert hdu.columns._arrays[1][0] == 200
        assert hdu.columns.columns[1].array[0] == 200
        assert hdu.data[0][1] == 200

        hdu.data._coldefs.columns[1].array[0] = 100
        assert hdu.data._coldefs._arrays[1][0] == 100
        assert hdu.data._coldefs.columns[1].array[0] == 100
        assert hdu.columns._arrays[1][0] == 100
        assert hdu.columns.columns[1].array[0] == 100
        assert hdu.data[0][1] == 100

        hdu.columns._arrays[1][0] = 90
        assert hdu.data._coldefs._arrays[1][0] == 90
        assert hdu.data._coldefs.columns[1].array[0] == 90
        assert hdu.columns._arrays[1][0] == 90
        assert hdu.columns.columns[1].array[0] == 90
        assert hdu.data[0][1] == 90

        hdu.columns.columns[1].array[0] = 80
        assert hdu.data._coldefs._arrays[1][0] == 80
        assert hdu.data._coldefs.columns[1].array[0] == 80
        assert hdu.columns._arrays[1][0] == 80
        assert hdu.columns.columns[1].array[0] == 80
        assert hdu.data[0][1] == 80

        t1.close()
        t2.close()
        hdul.close()

    def test_modify_column_attributes(self):
        """Regression test for https://github.com/astropy/astropy/issues/996

        This just tests one particular use case, but it should apply pretty
        well to other similar cases.
        """

        NULLS = {'a': 2, 'b': 'b', 'c': 2.3}

        data = np.array(list(zip([1, 2, 3, 4],
                                 ['a', 'b', 'c', 'd'],
                                 [2.3, 4.5, 6.7, 8.9])),
                        dtype=[('a', int), ('b', 'S1'), ('c', float)])

        b = fits.BinTableHDU(data=data)
        for col in b.columns:
            col.null = NULLS[col.name]

        b.writeto(self.temp('test.fits'), overwrite=True)

        with fits.open(self.temp('test.fits')) as hdul:
            header = hdul[1].header
            assert header['TNULL1'] == 2
            assert header['TNULL2'] == 'b'
            assert header['TNULL3'] == 2.3

    def test_multidimension_table_from_numpy_rec_columns(self):
        """Regression test for https://github.com/astropy/astropy/issues/5280
        and https://github.com/astropy/astropy/issues/5287

        multidimentional tables can now be written with the correct TDIM.
        Author: Stephen Bailey.
        """

        dtype = [
            (str('x'), (str, 5)),        # 1D column of 5-character strings
            (str('y'), (str, 3), (4,)),  # 2D column; each row is four 3-char strings
        ]
        data = np.zeros(2, dtype=dtype)
        data['x'] = ['abcde', 'xyz']
        data['y'][0] = ['A', 'BC', 'DEF', '123']
        data['y'][1] = ['X', 'YZ', 'PQR', '999']
        table = Table(data)

        # Test convenience functions io.fits.writeto / getdata
        fits.writeto(self.temp('test.fits'), data)
        dx = fits.getdata(self.temp('test.fits'))
        assert data['x'].dtype == dx['x'].dtype
        assert data['y'].dtype == dx['y'].dtype
        assert np.all(data['x'] == dx['x']), 'x: {} != {}'.format(data['x'], dx['x'])
        assert np.all(data['y'] == dx['y']), 'y: {} != {}'.format(data['y'], dx['y'])

        # Test fits.BinTableHDU(data) and avoid convenience functions
        hdu0 = fits.PrimaryHDU()
        hdu1 = fits.BinTableHDU(data)
        hx = fits.HDUList([hdu0, hdu1])
        hx.writeto(self.temp('test2.fits'))
        fx = fits.open(self.temp('test2.fits'))
        dx = fx[1].data
        fx.close()
        assert data['x'].dtype == dx['x'].dtype
        assert data['y'].dtype == dx['y'].dtype
        assert np.all(data['x'] == dx['x']), 'x: {} != {}'.format(data['x'], dx['x'])
        assert np.all(data['y'] == dx['y']), 'y: {} != {}'.format(data['y'], dx['y'])

        # Test Table write and read
        table.write(self.temp('test3.fits'))
        tx = Table.read(self.temp('test3.fits'), character_as_bytes=False)
        assert table['x'].dtype == tx['x'].dtype
        assert table['y'].dtype == tx['y'].dtype
        assert np.all(table['x'] == tx['x']), 'x: {} != {}'.format(table['x'], tx['x'])
        assert np.all(table['y'] == tx['y']), 'y: {} != {}'.format(table['y'], tx['y'])

    def test_mask_array(self):
        t = fits.open(self.data('table.fits'))
        tbdata = t[1].data
        mask = tbdata.field('V_mag') > 12
        newtbdata = tbdata[mask]
        hdu = fits.BinTableHDU(newtbdata)
        hdu.writeto(self.temp('newtable.fits'))

        hdul = fits.open(self.temp('newtable.fits'))

        # match to a regex rather than a specific string.
        expect = r"\[\('NGC1002',\s+12.3[0-9]*\) \(\'NGC1003\',\s+15.[0-9]+\)\]"
        assert re.match(expect, str(hdu.data))
        assert re.match(expect, str(hdul[1].data))

        t.close()
        hdul.close()

    def test_slice_a_row(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        t1 = fits.open(self.temp('table1.fits'))
        row = t1[1].data[2]
        assert row['counts'] == 308
        a, b, c = row[1:4]
        assert a == counts[2]
        assert b == ''
        assert (c == np.array([0., 0., 0., 0., 0.], dtype=np.float32)).all()
        row['counts'] = 310
        assert row['counts'] == 310

        row[1] = 315
        assert row['counts'] == 315

        assert row[1:4]['counts'] == 315

        pytest.raises(KeyError, lambda r: r[1:4]['flag'], row)

        row[1:4]['counts'] = 300
        assert row[1:4]['counts'] == 300
        assert row['counts'] == 300

        row[1:4][0] = 400
        assert row[1:4]['counts'] == 400
        row[1:4]['counts'] = 300
        assert row[1:4]['counts'] == 300

        # Test stepping for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/59
        row[1:4][::-1][-1] = 500
        assert row[1:4]['counts'] == 500
        row[1:4:2][0] = 300
        assert row[1:4]['counts'] == 300

        pytest.raises(KeyError, lambda r: r[1:4]['flag'], row)

        assert row[1:4].field(0) == 300
        assert row[1:4].field('counts') == 300

        pytest.raises(KeyError, row[1:4].field, 'flag')

        row[1:4].setfield('counts', 500)
        assert row[1:4].field(0) == 500

        pytest.raises(KeyError, row[1:4].setfield, 'flag', False)

        assert t1[1].data._coldefs._arrays[1][2] == 500
        assert t1[1].data._coldefs.columns[1].array[2] == 500
        assert t1[1].columns._arrays[1][2] == 500
        assert t1[1].columns.columns[1].array[2] == 500
        assert t1[1].data[2][1] == 500

        t1.close()

    def test_fits_record_len(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        t1 = fits.open(self.temp('table1.fits'))

        assert len(t1[1].data[0]) == 5
        assert len(t1[1].data[0][0:4]) == 4
        assert len(t1[1].data[0][0:5]) == 5
        assert len(t1[1].data[0][0:6]) == 5
        assert len(t1[1].data[0][0:7]) == 5
        assert len(t1[1].data[0][1:4]) == 3
        assert len(t1[1].data[0][1:5]) == 4
        assert len(t1[1].data[0][1:6]) == 4
        assert len(t1[1].data[0][1:7]) == 4

        t1.close()

    def test_add_data_by_rows(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu1 = fits.BinTableHDU.from_columns(coldefs)

        c1 = fits.Column(name='target', format='10A')
        c2 = fits.Column(name='counts', format='J', unit='DN')
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L')
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.BinTableHDU.from_columns(coldefs, nrows=5)

        # Test assigning data to a tables row using a FITS_record
        tbhdu.data[0] = tbhdu1.data[0]
        tbhdu.data[4] = tbhdu1.data[3]

        # Test assigning data to a tables row using a tuple
        tbhdu.data[2] = ('NGC1', 312, 'A Note',
                         np.array([1.1, 2.2, 3.3, 4.4, 5.5], dtype=np.float32),
                         True)

        # Test assigning data to a tables row using a list
        tbhdu.data[3] = ['JIM1', '33', 'A Note',
                         np.array([1., 2., 3., 4., 5.], dtype=np.float32),
                         True]

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert (id(tbhdu.data._coldefs.columns[0].array) ==
                id(tbhdu.data._coldefs._arrays[0]))
        assert (id(tbhdu.data._coldefs.columns[0].array) ==
                id(tbhdu.columns.columns[0].array))
        assert (id(tbhdu.data._coldefs.columns[0].array) ==
                id(tbhdu.columns._arrays[0]))

        assert tbhdu.data[0][1] == 312
        assert tbhdu.data._coldefs._arrays[1][0] == 312
        assert tbhdu.data._coldefs.columns[1].array[0] == 312
        assert tbhdu.columns._arrays[1][0] == 312
        assert tbhdu.columns.columns[1].array[0] == 312
        assert tbhdu.columns.columns[0].array[0] == 'NGC1'
        assert tbhdu.columns.columns[2].array[0] == ''
        assert (tbhdu.columns.columns[3].array[0] ==
                np.array([0., 0., 0., 0., 0.], dtype=np.float32)).all()
        assert tbhdu.columns.columns[4].array[0] == True  # noqa

        assert tbhdu.data[3][1] == 33
        assert tbhdu.data._coldefs._arrays[1][3] == 33
        assert tbhdu.data._coldefs.columns[1].array[3] == 33
        assert tbhdu.columns._arrays[1][3] == 33
        assert tbhdu.columns.columns[1].array[3] == 33
        assert tbhdu.columns.columns[0].array[3] == 'JIM1'
        assert tbhdu.columns.columns[2].array[3] == 'A Note'
        assert (tbhdu.columns.columns[3].array[3] ==
                np.array([1., 2., 3., 4., 5.], dtype=np.float32)).all()
        assert tbhdu.columns.columns[4].array[3] == True  # noqa

    def test_assign_multiple_rows_to_table(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu1 = fits.BinTableHDU.from_columns(coldefs)

        counts = np.array([112, 134, 108, 117])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[0, 1, 0, 0])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)
        tbhdu.data[0][3] = np.array([1., 2., 3., 4., 5.], dtype=np.float32)

        tbhdu2 = fits.BinTableHDU.from_columns(tbhdu1.data, nrows=9)

        # Assign the 4 rows from the second table to rows 5 thru 8 of the
        # new table.  Note that the last row of the new table will still be
        # initialized to the default values.
        tbhdu2.data[4:] = tbhdu.data

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert (id(tbhdu2.data._coldefs.columns[0].array) ==
                id(tbhdu2.data._coldefs._arrays[0]))
        assert (id(tbhdu2.data._coldefs.columns[0].array) ==
                id(tbhdu2.columns.columns[0].array))
        assert (id(tbhdu2.data._coldefs.columns[0].array) ==
                id(tbhdu2.columns._arrays[0]))

        assert tbhdu2.data[0][1] == 312
        assert tbhdu2.data._coldefs._arrays[1][0] == 312
        assert tbhdu2.data._coldefs.columns[1].array[0] == 312
        assert tbhdu2.columns._arrays[1][0] == 312
        assert tbhdu2.columns.columns[1].array[0] == 312
        assert tbhdu2.columns.columns[0].array[0] == 'NGC1'
        assert tbhdu2.columns.columns[2].array[0] == ''
        assert (tbhdu2.columns.columns[3].array[0] ==
                np.array([0., 0., 0., 0., 0.], dtype=np.float32)).all()
        assert tbhdu2.columns.columns[4].array[0] == True  # noqa

        assert tbhdu2.data[4][1] == 112
        assert tbhdu2.data._coldefs._arrays[1][4] == 112
        assert tbhdu2.data._coldefs.columns[1].array[4] == 112
        assert tbhdu2.columns._arrays[1][4] == 112
        assert tbhdu2.columns.columns[1].array[4] == 112
        assert tbhdu2.columns.columns[0].array[4] == 'NGC5'
        assert tbhdu2.columns.columns[2].array[4] == ''
        assert (tbhdu2.columns.columns[3].array[4] ==
                np.array([1., 2., 3., 4., 5.], dtype=np.float32)).all()
        assert tbhdu2.columns.columns[4].array[4] == False  # noqa
        assert tbhdu2.columns.columns[1].array[8] == 0
        assert tbhdu2.columns.columns[0].array[8] == ''
        assert tbhdu2.columns.columns[2].array[8] == ''
        assert (tbhdu2.columns.columns[3].array[8] ==
                np.array([0., 0., 0., 0., 0.], dtype=np.float32)).all()
        assert tbhdu2.columns.columns[4].array[8] == False  # noqa

    def test_verify_data_references(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        # Verify that original ColDefs object has independent Column
        # objects.
        assert id(coldefs.columns[0]) != id(c1)

        # Verify that original ColDefs object has independent ndarray
        # objects.
        assert id(coldefs.columns[0].array) != id(names)

        # Verify that original ColDefs object references the same data
        # object as the original Column object.
        assert id(coldefs.columns[0].array) == id(c1.array)
        assert id(coldefs.columns[0].array) == id(coldefs._arrays[0])

        # Verify new HDU has an independent ColDefs object.
        assert id(coldefs) != id(tbhdu.columns)

        # Verify new HDU has independent Column objects.
        assert id(coldefs.columns[0]) != id(tbhdu.columns.columns[0])

        # Verify new HDU has independent ndarray objects.
        assert (id(coldefs.columns[0].array) !=
                id(tbhdu.columns.columns[0].array))

        # Verify that both ColDefs objects in the HDU reference the same
        # Coldefs object.
        assert id(tbhdu.columns) == id(tbhdu.data._coldefs)

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert (id(tbhdu.data._coldefs.columns[0].array) ==
                id(tbhdu.data._coldefs._arrays[0]))
        assert (id(tbhdu.data._coldefs.columns[0].array) ==
                id(tbhdu.columns.columns[0].array))
        assert (id(tbhdu.data._coldefs.columns[0].array) ==
                id(tbhdu.columns._arrays[0]))

        tbhdu.writeto(self.temp('table1.fits'))

        t1 = fits.open(self.temp('table1.fits'))

        t1[1].data[0][1] = 213

        assert t1[1].data[0][1] == 213
        assert t1[1].data._coldefs._arrays[1][0] == 213
        assert t1[1].data._coldefs.columns[1].array[0] == 213
        assert t1[1].columns._arrays[1][0] == 213
        assert t1[1].columns.columns[1].array[0] == 213

        t1[1].data._coldefs._arrays[1][0] = 100

        assert t1[1].data[0][1] == 100
        assert t1[1].data._coldefs._arrays[1][0] == 100
        assert t1[1].data._coldefs.columns[1].array[0] == 100
        assert t1[1].columns._arrays[1][0] == 100
        assert t1[1].columns.columns[1].array[0] == 100

        t1[1].data._coldefs.columns[1].array[0] = 500
        assert t1[1].data[0][1] == 500
        assert t1[1].data._coldefs._arrays[1][0] == 500
        assert t1[1].data._coldefs.columns[1].array[0] == 500
        assert t1[1].columns._arrays[1][0] == 500
        assert t1[1].columns.columns[1].array[0] == 500

        t1[1].columns._arrays[1][0] = 600
        assert t1[1].data[0][1] == 600
        assert t1[1].data._coldefs._arrays[1][0] == 600
        assert t1[1].data._coldefs.columns[1].array[0] == 600
        assert t1[1].columns._arrays[1][0] == 600
        assert t1[1].columns.columns[1].array[0] == 600

        t1[1].columns.columns[1].array[0] = 800
        assert t1[1].data[0][1] == 800
        assert t1[1].data._coldefs._arrays[1][0] == 800
        assert t1[1].data._coldefs.columns[1].array[0] == 800
        assert t1[1].columns._arrays[1][0] == 800
        assert t1[1].columns.columns[1].array[0] == 800

        t1.close()

    def test_new_table_with_ndarray(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        tbhdu1 = fits.BinTableHDU.from_columns(tbhdu.data.view(np.ndarray))

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert (id(tbhdu1.data._coldefs.columns[0].array) ==
                id(tbhdu1.data._coldefs._arrays[0]))
        assert (id(tbhdu1.data._coldefs.columns[0].array) ==
                id(tbhdu1.columns.columns[0].array))
        assert (id(tbhdu1.data._coldefs.columns[0].array) ==
                id(tbhdu1.columns._arrays[0]))

        # Ensure I can change the value of one data element and it effects
        # all of the others.
        tbhdu1.data[0][1] = 213

        assert tbhdu1.data[0][1] == 213
        assert tbhdu1.data._coldefs._arrays[1][0] == 213
        assert tbhdu1.data._coldefs.columns[1].array[0] == 213
        assert tbhdu1.columns._arrays[1][0] == 213
        assert tbhdu1.columns.columns[1].array[0] == 213

        tbhdu1.data._coldefs._arrays[1][0] = 100

        assert tbhdu1.data[0][1] == 100
        assert tbhdu1.data._coldefs._arrays[1][0] == 100
        assert tbhdu1.data._coldefs.columns[1].array[0] == 100
        assert tbhdu1.columns._arrays[1][0] == 100
        assert tbhdu1.columns.columns[1].array[0] == 100

        tbhdu1.data._coldefs.columns[1].array[0] = 500
        assert tbhdu1.data[0][1] == 500
        assert tbhdu1.data._coldefs._arrays[1][0] == 500
        assert tbhdu1.data._coldefs.columns[1].array[0] == 500
        assert tbhdu1.columns._arrays[1][0] == 500
        assert tbhdu1.columns.columns[1].array[0] == 500

        tbhdu1.columns._arrays[1][0] = 600
        assert tbhdu1.data[0][1] == 600
        assert tbhdu1.data._coldefs._arrays[1][0] == 600
        assert tbhdu1.data._coldefs.columns[1].array[0] == 600
        assert tbhdu1.columns._arrays[1][0] == 600
        assert tbhdu1.columns.columns[1].array[0] == 600

        tbhdu1.columns.columns[1].array[0] = 800
        assert tbhdu1.data[0][1] == 800
        assert tbhdu1.data._coldefs._arrays[1][0] == 800
        assert tbhdu1.data._coldefs.columns[1].array[0] == 800
        assert tbhdu1.columns._arrays[1][0] == 800
        assert tbhdu1.columns.columns[1].array[0] == 800

        tbhdu1.writeto(self.temp('table1.fits'))

        t1 = fits.open(self.temp('table1.fits'))

        t1[1].data[0][1] = 213

        assert t1[1].data[0][1] == 213
        assert t1[1].data._coldefs._arrays[1][0] == 213
        assert t1[1].data._coldefs.columns[1].array[0] == 213
        assert t1[1].columns._arrays[1][0] == 213
        assert t1[1].columns.columns[1].array[0] == 213

        t1[1].data._coldefs._arrays[1][0] = 100

        assert t1[1].data[0][1] == 100
        assert t1[1].data._coldefs._arrays[1][0] == 100
        assert t1[1].data._coldefs.columns[1].array[0] == 100
        assert t1[1].columns._arrays[1][0] == 100
        assert t1[1].columns.columns[1].array[0] == 100

        t1[1].data._coldefs.columns[1].array[0] = 500
        assert t1[1].data[0][1] == 500
        assert t1[1].data._coldefs._arrays[1][0] == 500
        assert t1[1].data._coldefs.columns[1].array[0] == 500
        assert t1[1].columns._arrays[1][0] == 500
        assert t1[1].columns.columns[1].array[0] == 500

        t1[1].columns._arrays[1][0] = 600
        assert t1[1].data[0][1] == 600
        assert t1[1].data._coldefs._arrays[1][0] == 600
        assert t1[1].data._coldefs.columns[1].array[0] == 600
        assert t1[1].columns._arrays[1][0] == 600
        assert t1[1].columns.columns[1].array[0] == 600

        t1[1].columns.columns[1].array[0] = 800
        assert t1[1].data[0][1] == 800
        assert t1[1].data._coldefs._arrays[1][0] == 800
        assert t1[1].data._coldefs.columns[1].array[0] == 800
        assert t1[1].columns._arrays[1][0] == 800
        assert t1[1].columns.columns[1].array[0] == 800

        t1.close()

    def test_new_table_with_fits_rec(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.BinTableHDU.from_columns(coldefs)

        tbhdu.data[0][1] = 213

        assert tbhdu.data[0][1] == 213
        assert tbhdu.data._coldefs._arrays[1][0] == 213
        assert tbhdu.data._coldefs.columns[1].array[0] == 213
        assert tbhdu.columns._arrays[1][0] == 213
        assert tbhdu.columns.columns[1].array[0] == 213

        tbhdu.data._coldefs._arrays[1][0] = 100

        assert tbhdu.data[0][1] == 100
        assert tbhdu.data._coldefs._arrays[1][0] == 100
        assert tbhdu.data._coldefs.columns[1].array[0] == 100
        assert tbhdu.columns._arrays[1][0] == 100
        assert tbhdu.columns.columns[1].array[0] == 100

        tbhdu.data._coldefs.columns[1].array[0] = 500
        assert tbhdu.data[0][1] == 500
        assert tbhdu.data._coldefs._arrays[1][0] == 500
        assert tbhdu.data._coldefs.columns[1].array[0] == 500
        assert tbhdu.columns._arrays[1][0] == 500
        assert tbhdu.columns.columns[1].array[0] == 500

        tbhdu.columns._arrays[1][0] = 600
        assert tbhdu.data[0][1] == 600
        assert tbhdu.data._coldefs._arrays[1][0] == 600
        assert tbhdu.data._coldefs.columns[1].array[0] == 600
        assert tbhdu.columns._arrays[1][0] == 600
        assert tbhdu.columns.columns[1].array[0] == 600

        tbhdu.columns.columns[1].array[0] = 800
        assert tbhdu.data[0][1] == 800
        assert tbhdu.data._coldefs._arrays[1][0] == 800
        assert tbhdu.data._coldefs.columns[1].array[0] == 800
        assert tbhdu.columns._arrays[1][0] == 800
        assert tbhdu.columns.columns[1].array[0] == 800

        tbhdu.columns.columns[1].array[0] = 312

        tbhdu.writeto(self.temp('table1.fits'))

        t1 = fits.open(self.temp('table1.fits'))

        t1[1].data[0][1] = 1
        fr = t1[1].data
        assert t1[1].data[0][1] == 1
        assert t1[1].data._coldefs._arrays[1][0] == 1
        assert t1[1].data._coldefs.columns[1].array[0] == 1
        assert t1[1].columns._arrays[1][0] == 1
        assert t1[1].columns.columns[1].array[0] == 1
        assert fr[0][1] == 1
        assert fr._coldefs._arrays[1][0] == 1
        assert fr._coldefs.columns[1].array[0] == 1

        fr._coldefs.columns[1].array[0] = 312

        tbhdu1 = fits.BinTableHDU.from_columns(fr)

        i = 0
        for row in tbhdu1.data:
            for j in range(len(row)):
                if isinstance(row[j], np.ndarray):
                    assert (row[j] == tbhdu.data[i][j]).all()
                else:
                    assert row[j] == tbhdu.data[i][j]
            i = i + 1

        tbhdu1.data[0][1] = 213

        assert t1[1].data[0][1] == 312
        assert t1[1].data._coldefs._arrays[1][0] == 312
        assert t1[1].data._coldefs.columns[1].array[0] == 312
        assert t1[1].columns._arrays[1][0] == 312
        assert t1[1].columns.columns[1].array[0] == 312
        assert fr[0][1] == 312
        assert fr._coldefs._arrays[1][0] == 312
        assert fr._coldefs.columns[1].array[0] == 312
        assert tbhdu1.data[0][1] == 213
        assert tbhdu1.data._coldefs._arrays[1][0] == 213
        assert tbhdu1.data._coldefs.columns[1].array[0] == 213
        assert tbhdu1.columns._arrays[1][0] == 213
        assert tbhdu1.columns.columns[1].array[0] == 213

        t1[1].data[0][1] = 10

        assert t1[1].data[0][1] == 10
        assert t1[1].data._coldefs._arrays[1][0] == 10
        assert t1[1].data._coldefs.columns[1].array[0] == 10
        assert t1[1].columns._arrays[1][0] == 10
        assert t1[1].columns.columns[1].array[0] == 10
        assert fr[0][1] == 10
        assert fr._coldefs._arrays[1][0] == 10
        assert fr._coldefs.columns[1].array[0] == 10
        assert tbhdu1.data[0][1] == 213
        assert tbhdu1.data._coldefs._arrays[1][0] == 213
        assert tbhdu1.data._coldefs.columns[1].array[0] == 213
        assert tbhdu1.columns._arrays[1][0] == 213
        assert tbhdu1.columns.columns[1].array[0] == 213

        tbhdu1.data._coldefs._arrays[1][0] = 666

        assert t1[1].data[0][1] == 10
        assert t1[1].data._coldefs._arrays[1][0] == 10
        assert t1[1].data._coldefs.columns[1].array[0] == 10
        assert t1[1].columns._arrays[1][0] == 10
        assert t1[1].columns.columns[1].array[0] == 10
        assert fr[0][1] == 10
        assert fr._coldefs._arrays[1][0] == 10
        assert fr._coldefs.columns[1].array[0] == 10
        assert tbhdu1.data[0][1] == 666
        assert tbhdu1.data._coldefs._arrays[1][0] == 666
        assert tbhdu1.data._coldefs.columns[1].array[0] == 666
        assert tbhdu1.columns._arrays[1][0] == 666
        assert tbhdu1.columns.columns[1].array[0] == 666

        t1.close()

    def test_bin_table_hdu_constructor(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu1 = fits.BinTableHDU.from_columns(coldefs)

        hdu = fits.BinTableHDU(tbhdu1.data)

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert (id(hdu.data._coldefs.columns[0].array) ==
                id(hdu.data._coldefs._arrays[0]))
        assert (id(hdu.data._coldefs.columns[0].array) ==
                id(hdu.columns.columns[0].array))
        assert (id(hdu.data._coldefs.columns[0].array) ==
                id(hdu.columns._arrays[0]))

        # Verify that the references in the original HDU are the same as the
        # references in the new HDU.
        assert (id(tbhdu1.data._coldefs.columns[0].array) ==
                id(hdu.data._coldefs._arrays[0]))

        # Verify that a change in the new HDU is reflected in both the new
        # and original HDU.

        hdu.data[0][1] = 213

        assert hdu.data[0][1] == 213
        assert hdu.data._coldefs._arrays[1][0] == 213
        assert hdu.data._coldefs.columns[1].array[0] == 213
        assert hdu.columns._arrays[1][0] == 213
        assert hdu.columns.columns[1].array[0] == 213
        assert tbhdu1.data[0][1] == 213
        assert tbhdu1.data._coldefs._arrays[1][0] == 213
        assert tbhdu1.data._coldefs.columns[1].array[0] == 213
        assert tbhdu1.columns._arrays[1][0] == 213
        assert tbhdu1.columns.columns[1].array[0] == 213

        hdu.data._coldefs._arrays[1][0] = 100

        assert hdu.data[0][1] == 100
        assert hdu.data._coldefs._arrays[1][0] == 100
        assert hdu.data._coldefs.columns[1].array[0] == 100
        assert hdu.columns._arrays[1][0] == 100
        assert hdu.columns.columns[1].array[0] == 100
        assert tbhdu1.data[0][1] == 100
        assert tbhdu1.data._coldefs._arrays[1][0] == 100
        assert tbhdu1.data._coldefs.columns[1].array[0] == 100
        assert tbhdu1.columns._arrays[1][0] == 100
        assert tbhdu1.columns.columns[1].array[0] == 100

        hdu.data._coldefs.columns[1].array[0] = 500
        assert hdu.data[0][1] == 500
        assert hdu.data._coldefs._arrays[1][0] == 500
        assert hdu.data._coldefs.columns[1].array[0] == 500
        assert hdu.columns._arrays[1][0] == 500
        assert hdu.columns.columns[1].array[0] == 500
        assert tbhdu1.data[0][1] == 500
        assert tbhdu1.data._coldefs._arrays[1][0] == 500
        assert tbhdu1.data._coldefs.columns[1].array[0] == 500
        assert tbhdu1.columns._arrays[1][0] == 500
        assert tbhdu1.columns.columns[1].array[0] == 500

        hdu.columns._arrays[1][0] = 600
        assert hdu.data[0][1] == 600
        assert hdu.data._coldefs._arrays[1][0] == 600
        assert hdu.data._coldefs.columns[1].array[0] == 600
        assert hdu.columns._arrays[1][0] == 600
        assert hdu.columns.columns[1].array[0] == 600
        assert tbhdu1.data[0][1] == 600
        assert tbhdu1.data._coldefs._arrays[1][0] == 600
        assert tbhdu1.data._coldefs.columns[1].array[0] == 600
        assert tbhdu1.columns._arrays[1][0] == 600
        assert tbhdu1.columns.columns[1].array[0] == 600

        hdu.columns.columns[1].array[0] = 800
        assert hdu.data[0][1] == 800
        assert hdu.data._coldefs._arrays[1][0] == 800
        assert hdu.data._coldefs.columns[1].array[0] == 800
        assert hdu.columns._arrays[1][0] == 800
        assert hdu.columns.columns[1].array[0] == 800
        assert tbhdu1.data[0][1] == 800
        assert tbhdu1.data._coldefs._arrays[1][0] == 800
        assert tbhdu1.data._coldefs.columns[1].array[0] == 800
        assert tbhdu1.columns._arrays[1][0] == 800
        assert tbhdu1.columns.columns[1].array[0] == 800

    def test_constructor_name_arg(self):
        """testConstructorNameArg

        Passing name='...' to the BinTableHDU and TableHDU constructors
        should set the .name attribute and 'EXTNAME' header keyword, and
        override any name in an existing 'EXTNAME' value.
        """

        for hducls in [fits.BinTableHDU, fits.TableHDU]:
            # First test some default assumptions
            hdu = hducls()
            assert hdu.name == ''
            assert 'EXTNAME' not in hdu.header
            hdu.name = 'FOO'
            assert hdu.name == 'FOO'
            assert hdu.header['EXTNAME'] == 'FOO'

            # Passing name to constructor
            hdu = hducls(name='FOO')
            assert hdu.name == 'FOO'
            assert hdu.header['EXTNAME'] == 'FOO'

            # And overriding a header with a different extname
            hdr = fits.Header()
            hdr['EXTNAME'] = 'EVENTS'
            hdu = hducls(header=hdr, name='FOO')
            assert hdu.name == 'FOO'
            assert hdu.header['EXTNAME'] == 'FOO'

    def test_constructor_ver_arg(self):
        for hducls in [fits.BinTableHDU, fits.TableHDU]:
            # First test some default assumptions
            hdu = hducls()
            assert hdu.ver == 1
            assert 'EXTVER' not in hdu.header
            hdu.ver = 2
            assert hdu.ver == 2
            assert hdu.header['EXTVER'] == 2

            # Passing name to constructor
            hdu = hducls(ver=3)
            assert hdu.ver == 3
            assert hdu.header['EXTVER'] == 3

            # And overriding a header with a different extver
            hdr = fits.Header()
            hdr['EXTVER'] = 4
            hdu = hducls(header=hdr, ver=5)
            assert hdu.ver == 5
            assert hdu.header['EXTVER'] == 5

    def test_unicode_colname(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/5204
        "Handle unicode FITS BinTable column names on Python 2"
        """
        col = fits.Column(name='spam', format='E', array=[42.])
        # This used to raise a TypeError, now it works
        fits.BinTableHDU.from_columns([col])

    def test_bin_table_with_logical_array(self):
        c1 = fits.Column(name='flag', format='2L',
                         array=[[True, False], [False, True]])
        coldefs = fits.ColDefs([c1])

        tbhdu1 = fits.BinTableHDU.from_columns(coldefs)

        assert (tbhdu1.data.field('flag')[0] ==
                np.array([True, False], dtype=bool)).all()
        assert (tbhdu1.data.field('flag')[1] ==
                np.array([False, True], dtype=bool)).all()

        tbhdu = fits.BinTableHDU.from_columns(tbhdu1.data)

        assert (tbhdu.data.field('flag')[0] ==
                np.array([True, False], dtype=bool)).all()
        assert (tbhdu.data.field('flag')[1] ==
                np.array([False, True], dtype=bool)).all()

    def test_fits_rec_column_access(self):
        tbdata = fits.getdata(self.data('table.fits'))
        assert (tbdata.V_mag == tbdata.field('V_mag')).all()
        assert (tbdata.V_mag == tbdata['V_mag']).all()

        # Table with scaling (c3) and tnull (c1)
        tbdata = fits.getdata(self.data('tb.fits'))
        for col in ('c1', 'c2', 'c3', 'c4'):
            data = getattr(tbdata, col)
            assert (data == tbdata.field(col)).all()
            assert (data == tbdata[col]).all()

        # ascii table
        tbdata = fits.getdata(self.data('ascii.fits'))
        for col in ('a', 'b'):
            data = getattr(tbdata, col)
            assert (data == tbdata.field(col)).all()
            assert (data == tbdata[col]).all()

        # with VLA column
        col1 = fits.Column(name='x', format='PI()',
                           array=np.array([[45, 56], [11, 12, 13]],
                                          dtype=np.object_))
        hdu = fits.BinTableHDU.from_columns([col1])
        assert type(hdu.data['x']) == type(hdu.data.x)  # noqa
        assert (hdu.data['x'][0] == hdu.data.x[0]).all()
        assert (hdu.data['x'][1] == hdu.data.x[1]).all()

    def test_table_with_zero_width_column(self):
        hdul = fits.open(self.data('zerowidth.fits'))
        tbhdu = hdul[2]  # This HDU contains a zero-width column 'ORBPARM'
        assert 'ORBPARM' in tbhdu.columns.names
        # The ORBPARM column should not be in the data, though the data should
        # be readable
        assert 'ORBPARM' in tbhdu.data.names
        assert 'ORBPARM' in tbhdu.data.dtype.names
        # Verify that some of the data columns are still correctly accessible
        # by name
        assert tbhdu.data[0]['ANNAME'] == 'VLA:_W16'
        assert comparefloats(
            tbhdu.data[0]['STABXYZ'],
            np.array([499.85566663, -1317.99231554, -735.18866164],
                     dtype=np.float64))
        assert tbhdu.data[0]['NOSTA'] == 1
        assert tbhdu.data[0]['MNTSTA'] == 0
        assert tbhdu.data[-1]['ANNAME'] == 'VPT:_OUT'
        assert comparefloats(
            tbhdu.data[-1]['STABXYZ'],
            np.array([0.0, 0.0, 0.0], dtype=np.float64))
        assert tbhdu.data[-1]['NOSTA'] == 29
        assert tbhdu.data[-1]['MNTSTA'] == 0
        hdul.writeto(self.temp('newtable.fits'))
        hdul.close()
        hdul = fits.open(self.temp('newtable.fits'))
        tbhdu = hdul[2]

        # Verify that the previous tests still hold after writing
        assert 'ORBPARM' in tbhdu.columns.names
        assert 'ORBPARM' in tbhdu.data.names
        assert 'ORBPARM' in tbhdu.data.dtype.names
        assert tbhdu.data[0]['ANNAME'] == 'VLA:_W16'
        assert comparefloats(
            tbhdu.data[0]['STABXYZ'],
            np.array([499.85566663, -1317.99231554, -735.18866164],
                     dtype=np.float64))
        assert tbhdu.data[0]['NOSTA'] == 1
        assert tbhdu.data[0]['MNTSTA'] == 0
        assert tbhdu.data[-1]['ANNAME'] == 'VPT:_OUT'
        assert comparefloats(
            tbhdu.data[-1]['STABXYZ'],
            np.array([0.0, 0.0, 0.0], dtype=np.float64))
        assert tbhdu.data[-1]['NOSTA'] == 29
        assert tbhdu.data[-1]['MNTSTA'] == 0
        hdul.close()

    def test_string_column_padding(self):
        a = ['img1', 'img2', 'img3a', 'p']
        s = 'img1\x00\x00\x00\x00\x00\x00' \
            'img2\x00\x00\x00\x00\x00\x00' \
            'img3a\x00\x00\x00\x00\x00' \
            'p\x00\x00\x00\x00\x00\x00\x00\x00\x00'

        acol = fits.Column(name='MEMNAME', format='A10',
                           array=chararray.array(a))
        ahdu = fits.BinTableHDU.from_columns([acol])
        assert ahdu.data.tobytes().decode('raw-unicode-escape') == s
        ahdu.writeto(self.temp('newtable.fits'))
        with fits.open(self.temp('newtable.fits')) as hdul:
            assert hdul[1].data.tobytes().decode('raw-unicode-escape') == s
            assert (hdul[1].data['MEMNAME'] == a).all()
        del hdul

        ahdu = fits.TableHDU.from_columns([acol])
        ahdu.writeto(self.temp('newtable.fits'), overwrite=True)

        with fits.open(self.temp('newtable.fits')) as hdul:
            assert (hdul[1].data.tobytes().decode('raw-unicode-escape') ==
                    s.replace('\x00', ' '))
            assert (hdul[1].data['MEMNAME'] == a).all()
            ahdu = fits.BinTableHDU.from_columns(hdul[1].data.copy())
        del hdul

        # Now serialize once more as a binary table; padding bytes should
        # revert to zeroes
        ahdu.writeto(self.temp('newtable.fits'), overwrite=True)
        with fits.open(self.temp('newtable.fits')) as hdul:
            assert hdul[1].data.tobytes().decode('raw-unicode-escape') == s
            assert (hdul[1].data['MEMNAME'] == a).all()

    def test_multi_dimensional_columns(self):
        """
        Tests the multidimensional column implementation with both numeric
        arrays and string arrays.
        """

        data = np.rec.array(
            [([0, 1, 2, 3, 4, 5], 'row1' * 2),
             ([6, 7, 8, 9, 0, 1], 'row2' * 2),
             ([2, 3, 4, 5, 6, 7], 'row3' * 2)], formats='6i4,a8')

        thdu = fits.BinTableHDU.from_columns(data)

        thdu.writeto(self.temp('newtable.fits'))

        with fits.open(self.temp('newtable.fits'), mode='update') as hdul:
            # Modify the TDIM fields to my own specification
            hdul[1].header['TDIM1'] = '(2,3)'
            hdul[1].header['TDIM2'] = '(4,2)'

        with fits.open(self.temp('newtable.fits')) as hdul:
            thdu = hdul[1]

            c1 = thdu.data.field(0)
            c2 = thdu.data.field(1)

            assert c1.shape == (3, 3, 2)
            assert c2.shape == (3, 2)
            assert (c1 == np.array([[[0, 1], [2, 3], [4, 5]],
                                    [[6, 7], [8, 9], [0, 1]],
                                    [[2, 3], [4, 5], [6, 7]]])).all()
            assert (c2 == np.array([['row1', 'row1'],
                                    ['row2', 'row2'],
                                    ['row3', 'row3']])).all()
        del c1
        del c2
        del thdu
        del hdul

        # Test setting the TDIMn header based on the column data
        data = np.zeros(3, dtype=[('x', 'f4'), ('s', 'S5', 4)])
        data['x'] = 1, 2, 3
        data['s'] = 'ok'
        fits.writeto(self.temp('newtable.fits'), data, overwrite=True)

        t = fits.getdata(self.temp('newtable.fits'))

        assert t.field(1).dtype.str[-1] == '5'
        assert t.field(1).shape == (3, 4)

        # Like the previous test, but with an extra dimension (a bit more
        # complicated)
        data = np.zeros(3, dtype=[('x', 'f4'), ('s', 'S5', (4, 3))])
        data['x'] = 1, 2, 3
        data['s'] = 'ok'

        del t

        fits.writeto(self.temp('newtable.fits'), data, overwrite=True)

        t = fits.getdata(self.temp('newtable.fits'))

        assert t.field(1).dtype.str[-1] == '5'
        assert t.field(1).shape == (3, 4, 3)

    def test_oned_array_single_element(self):
        # a table with rows that are 1d arrays of a single value
        data = np.array([(1, ), (2, )], dtype=([('x', 'i4', (1, ))]))
        thdu = fits.BinTableHDU.from_columns(data)

        thdu.writeto(self.temp('onedtable.fits'))

        with fits.open(self.temp('onedtable.fits')) as hdul:
            thdu = hdul[1]

            c = thdu.data.field(0)
            assert c.shape == (2, 1)
            assert thdu.header['TDIM1'] == '(1)'

    def test_bin_table_init_from_string_array_column(self):
        """
        Tests two ways of creating a new `BinTableHDU` from a column of
        string arrays.

        This tests for a couple different regressions, and ensures that
        both BinTableHDU(data=arr) and BinTableHDU.from_columns(arr) work
        equivalently.

        Some of this is redundant with the following test, but checks some
        subtly different cases.
        """

        data = [[b'abcd', b'efgh'],
                [b'ijkl', b'mnop'],
                [b'qrst', b'uvwx']]

        arr = np.array([(data,), (data,), (data,), (data,), (data,)],
                       dtype=[('S', '(3, 2)S4')])

        tbhdu1 = fits.BinTableHDU(data=arr)

        def test_dims_and_roundtrip(tbhdu):
            assert tbhdu.data['S'].shape == (5, 3, 2)
            assert tbhdu.data['S'].dtype.str.endswith('U4')

            tbhdu.writeto(self.temp('test.fits'), overwrite=True)

            with fits.open(self.temp('test.fits')) as hdul:
                tbhdu2 = hdul[1]
                assert tbhdu2.header['TDIM1'] == '(4,2,3)'
                assert tbhdu2.data['S'].shape == (5, 3, 2)
                assert tbhdu.data['S'].dtype.str.endswith('U4')
                assert np.all(tbhdu2.data['S'] == tbhdu.data['S'])

        test_dims_and_roundtrip(tbhdu1)

        tbhdu2 = fits.BinTableHDU.from_columns(arr)
        test_dims_and_roundtrip(tbhdu2)

    def test_columns_with_truncating_tdim(self):
        """
        According to the FITS standard (section 7.3.2):

            If the number of elements in the array implied by the TDIMn is less
            than the allocated size of the ar- ray in the FITS file, then the
            unused trailing elements should be interpreted as containing
            undefined fill values.

        *deep sigh* What this means is if a column has a repeat count larger
        than the number of elements indicated by its TDIM (ex: TDIM1 = '(2,2)',
        but TFORM1 = 6I), then instead of this being an outright error we are
        to take the first 4 elements as implied by the TDIM and ignore the
        additional two trailing elements.
        """

        # It's hard to even successfully create a table like this.  I think
        # it *should* be difficult, but once created it should at least be
        # possible to read.
        arr1 = [[b'ab', b'cd'], [b'ef', b'gh'], [b'ij', b'kl']]
        arr2 = [1, 2, 3, 4, 5]

        arr = np.array([(arr1, arr2), (arr1, arr2)],
                       dtype=[('a', '(3, 2)S2'), ('b', '5i8')])

        tbhdu = fits.BinTableHDU(data=arr)
        tbhdu.writeto(self.temp('test.fits'))

        with open(self.temp('test.fits'), 'rb') as f:
            raw_bytes = f.read()

        # Artificially truncate TDIM in the header; this seems to be the
        # easiest way to do this while getting around Astropy's insistence on the
        # data and header matching perfectly; again, we have no interest in
        # making it possible to write files in this format, only read them
        with open(self.temp('test.fits'), 'wb') as f:
            f.write(raw_bytes.replace(b'(2,2,3)', b'(2,2,2)'))

        with fits.open(self.temp('test.fits')) as hdul:
            tbhdu2 = hdul[1]
            assert tbhdu2.header['TDIM1'] == '(2,2,2)'
            assert tbhdu2.header['TFORM1'] == '12A'
            for row in tbhdu2.data:
                assert np.all(row['a'] == [['ab', 'cd'], ['ef', 'gh']])
                assert np.all(row['b'] == [1, 2, 3, 4, 5])

    def test_string_array_round_trip(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/201"""

        data = [['abc', 'def', 'ghi'],
                ['jkl', 'mno', 'pqr'],
                ['stu', 'vwx', 'yz ']]

        recarr = np.rec.array([(data,), (data,)], formats=['(3,3)S3'])

        t = fits.BinTableHDU(data=recarr)
        t.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            assert 'TDIM1' in h[1].header
            assert h[1].header['TDIM1'] == '(3,3,3)'
            assert len(h[1].data) == 2
            assert len(h[1].data[0]) == 1
            assert (h[1].data.field(0)[0] ==
                    np.char.decode(recarr.field(0)[0], 'ascii')).all()

        with fits.open(self.temp('test.fits')) as h:
            # Access the data; I think this is necessary to exhibit the bug
            # reported in https://aeon.stsci.edu/ssb/trac/pyfits/ticket/201
            h[1].data[:]
            h.writeto(self.temp('test2.fits'))

        with fits.open(self.temp('test2.fits')) as h:
            assert 'TDIM1' in h[1].header
            assert h[1].header['TDIM1'] == '(3,3,3)'
            assert len(h[1].data) == 2
            assert len(h[1].data[0]) == 1
            assert (h[1].data.field(0)[0] ==
                    np.char.decode(recarr.field(0)[0], 'ascii')).all()

    def test_new_table_with_nd_column(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/3
        """

        arra = np.array(['a', 'b'], dtype='|S1')
        arrb = np.array([['a', 'bc'], ['cd', 'e']], dtype='|S2')
        arrc = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])

        cols = [
            fits.Column(name='str', format='1A', array=arra),
            fits.Column(name='strarray', format='4A', dim='(2,2)',
                        array=arrb),
            fits.Column(name='intarray', format='4I', dim='(2, 2)',
                        array=arrc)
        ]

        hdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            # Need to force string arrays to byte arrays in order to compare
            # correctly on Python 3
            assert (h[1].data['str'].encode('ascii') == arra).all()
            assert (h[1].data['strarray'].encode('ascii') == arrb).all()
            assert (h[1].data['intarray'] == arrc).all()

    def test_mismatched_tform_and_tdim(self):
        """Normally the product of the dimensions listed in a TDIMn keyword
        must be less than or equal to the repeat count in the TFORMn keyword.

        This tests that this works if less than (treating the trailing bytes
        as unspecified fill values per the FITS standard) and fails if the
        dimensions specified by TDIMn are greater than the repeat count.
        """

        arra = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
        arrb = np.array([[[9, 10], [11, 12]], [[13, 14], [15, 16]]])

        cols = [fits.Column(name='a', format='20I', dim='(2,2)',
                            array=arra),
                fits.Column(name='b', format='4I', dim='(2,2)',
                            array=arrb)]

        # The first column has the mismatched repeat count
        hdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            assert h[1].header['TFORM1'] == '20I'
            assert h[1].header['TFORM2'] == '4I'
            assert h[1].header['TDIM1'] == h[1].header['TDIM2'] == '(2,2)'
            assert (h[1].data['a'] == arra).all()
            assert (h[1].data['b'] == arrb).all()
            assert h[1].data.itemsize == 48  # 16-bits times 24

        # If dims is more than the repeat count in the format specifier raise
        # an error
        pytest.raises(VerifyError, fits.Column, name='a', format='2I',
                      dim='(2,2)', array=arra)

    def test_tdim_of_size_one(self):
        """Regression test for https://github.com/astropy/astropy/pull/3580"""

        with fits.open(self.data('tdim.fits')) as hdulist:
            assert hdulist[1].data['V_mag'].shape == (3, 1, 1)

    def test_slicing(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/52"""

        with fits.open(self.data('table.fits')) as f:
            data = f[1].data
        targets = data.field('target')
        s = data[:]
        assert (s.field('target') == targets).all()
        for n in range(len(targets) + 2):
            s = data[:n]
            assert (s.field('target') == targets[:n]).all()
            s = data[n:]
            assert (s.field('target') == targets[n:]).all()
        s = data[::2]
        assert (s.field('target') == targets[::2]).all()
        s = data[::-1]
        assert (s.field('target') == targets[::-1]).all()

    def test_array_slicing(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/55"""

        with fits.open(self.data('table.fits')) as f:
            data = f[1].data
        s1 = data[data['target'] == 'NGC1001']
        s2 = data[np.where(data['target'] == 'NGC1001')]
        s3 = data[[0]]
        s4 = data[:1]
        for s in [s1, s2, s3, s4]:
            assert isinstance(s, fits.FITS_rec)

        assert comparerecords(s1, s2)
        assert comparerecords(s2, s3)
        assert comparerecords(s3, s4)

    def test_array_broadcasting(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/pull/48
        """

        with fits.open(self.data('table.fits')) as hdu:
            data = hdu[1].data
            data['V_mag'] = 0
            assert np.all(data['V_mag'] == 0)

            data['V_mag'] = 1
            assert np.all(data['V_mag'] == 1)

            for container in (list, tuple, np.array):
                data['V_mag'] = container([1, 2, 3])
                assert np.array_equal(data['V_mag'], np.array([1, 2, 3]))

    def test_array_slicing_readonly(self):
        """
        Like test_array_slicing but with the file opened in 'readonly' mode.
        Regression test for a crash when slicing readonly memmap'd tables.
        """

        with fits.open(self.data('table.fits'), mode='readonly') as f:
            data = f[1].data
        s1 = data[data['target'] == 'NGC1001']
        s2 = data[np.where(data['target'] == 'NGC1001')]
        s3 = data[[0]]
        s4 = data[:1]
        for s in [s1, s2, s3, s4]:
            assert isinstance(s, fits.FITS_rec)
        assert comparerecords(s1, s2)
        assert comparerecords(s2, s3)
        assert comparerecords(s3, s4)

    @pytest.mark.parametrize('tablename', ['table.fits', 'tb.fits'])
    def test_dump_load_round_trip(self, tablename):
        """
        A simple test of the dump/load methods; dump the data, column, and
        header files and try to reload the table from them.
        """

        with fits.open(self.data(tablename)) as hdul:
            tbhdu = hdul[1]
            datafile = self.temp('data.txt')
            cdfile = self.temp('coldefs.txt')
            hfile = self.temp('header.txt')

            tbhdu.dump(datafile, cdfile, hfile)

            new_tbhdu = fits.BinTableHDU.load(datafile, cdfile, hfile)

            assert comparerecords(tbhdu.data, new_tbhdu.data)

            _assert_attr_col(new_tbhdu, hdul[1])

    def test_dump_load_array_colums(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/22

        Ensures that a table containing a multi-value array column can be
        dumped and loaded successfully.
        """

        data = np.rec.array([('a', [1, 2, 3, 4], 0.1),
                             ('b', [5, 6, 7, 8], 0.2)],
                            formats='a1,4i4,f8')
        tbhdu = fits.BinTableHDU.from_columns(data)
        datafile = self.temp('data.txt')
        cdfile = self.temp('coldefs.txt')
        hfile = self.temp('header.txt')

        tbhdu.dump(datafile, cdfile, hfile)
        new_tbhdu = fits.BinTableHDU.load(datafile, cdfile, hfile)
        assert comparerecords(tbhdu.data, new_tbhdu.data)
        assert str(tbhdu.header) == str(new_tbhdu.header)

    def test_load_guess_format(self):
        """
        Tests loading a table dump with no supplied coldefs or header, so that
        the table format has to be guessed at.  There is of course no exact
        science to this; the table that's produced simply uses sensible guesses
        for that format.  Ideally this should never have to be used.
        """

        # Create a table containing a variety of data types.
        a0 = np.array([False, True, False], dtype=bool)
        c0 = fits.Column(name='c0', format='L', array=a0)

        # Format X currently not supported by the format
        # a1 = np.array([[0], [1], [0]], dtype=np.uint8)
        # c1 = fits.Column(name='c1', format='X', array=a1)

        a2 = np.array([1, 128, 255], dtype=np.uint8)
        c2 = fits.Column(name='c2', format='B', array=a2)
        a3 = np.array([-30000, 1, 256], dtype=np.int16)
        c3 = fits.Column(name='c3', format='I', array=a3)
        a4 = np.array([-123123123, 1234, 123123123], dtype=np.int32)
        c4 = fits.Column(name='c4', format='J', array=a4)
        a5 = np.array(['a', 'abc', 'ab'])
        c5 = fits.Column(name='c5', format='A3', array=a5)
        a6 = np.array([1.1, 2.2, 3.3], dtype=np.float64)
        c6 = fits.Column(name='c6', format='D', array=a6)
        a7 = np.array([1.1 + 2.2j, 3.3 + 4.4j, 5.5 + 6.6j],
                      dtype=np.complex128)
        c7 = fits.Column(name='c7', format='M', array=a7)
        a8 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.int32)
        c8 = fits.Column(name='c8', format='PJ()', array=a8)

        tbhdu = fits.BinTableHDU.from_columns([c0, c2, c3, c4, c5, c6, c7, c8])

        datafile = self.temp('data.txt')
        tbhdu.dump(datafile)

        new_tbhdu = fits.BinTableHDU.load(datafile)

        # In this particular case the record data at least should be equivalent
        assert comparerecords(tbhdu.data, new_tbhdu.data)

    def test_attribute_field_shadowing(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/86

        Numpy recarray objects have a poorly-considered feature of allowing
        field access by attribute lookup.  However, if a field name coincides
        with an existing attribute/method of the array, the existing name takes
        presence (making the attribute-based field lookup completely unreliable
        in general cases).

        This ensures that any FITS_rec attributes still work correctly even
        when there is a field with the same name as that attribute.
        """

        c1 = fits.Column(name='names', format='I', array=[1])
        c2 = fits.Column(name='formats', format='I', array=[2])
        c3 = fits.Column(name='other', format='I', array=[3])

        t = fits.BinTableHDU.from_columns([c1, c2, c3])
        assert t.data.names == ['names', 'formats', 'other']
        assert t.data.formats == ['I'] * 3
        assert (t.data['names'] == [1]).all()
        assert (t.data['formats'] == [2]).all()
        assert (t.data.other == [3]).all()

    def test_table_from_bool_fields(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/113

        Tests creating a table from a recarray containing numpy.bool columns.
        """

        array = np.rec.array([(True, False), (False, True)], formats='|b1,|b1')
        thdu = fits.BinTableHDU.from_columns(array)
        assert thdu.columns.formats == ['L', 'L']
        assert comparerecords(thdu.data, array)

        # Test round trip
        thdu.writeto(self.temp('table.fits'))
        data = fits.getdata(self.temp('table.fits'), ext=1)
        assert thdu.columns.formats == ['L', 'L']
        assert comparerecords(data, array)

    def test_table_from_bool_fields2(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/215

        Tests the case where a multi-field ndarray (not a recarray) containing
        a bool field is used to initialize a `BinTableHDU`.
        """

        arr = np.array([(False,), (True,), (False,)], dtype=[('a', '?')])
        hdu = fits.BinTableHDU(data=arr)
        assert (hdu.data['a'] == arr['a']).all()

    def test_bool_column_update(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/139"""

        c1 = fits.Column('F1', 'L', array=[True, False])
        c2 = fits.Column('F2', 'L', array=[False, True])
        thdu = fits.BinTableHDU.from_columns(fits.ColDefs([c1, c2]))
        thdu.writeto(self.temp('table.fits'))

        with fits.open(self.temp('table.fits'), mode='update') as hdul:
            hdul[1].data['F1'][1] = True
            hdul[1].data['F2'][0] = True

        with fits.open(self.temp('table.fits')) as hdul:
            assert (hdul[1].data['F1'] == [True, True]).all()
            assert (hdul[1].data['F2'] == [True, True]).all()

    def test_missing_tnull(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/197"""

        c = fits.Column('F1', 'A3', null='---',
                        array=np.array(['1.0', '2.0', '---', '3.0']),
                        ascii=True)
        table = fits.TableHDU.from_columns([c])
        table.writeto(self.temp('test.fits'))

        # Now let's delete the TNULL1 keyword, making this essentially
        # unreadable
        with fits.open(self.temp('test.fits'), mode='update') as h:
            h[1].header['TFORM1'] = 'E3'
            del h[1].header['TNULL1']

        with fits.open(self.temp('test.fits')) as h:
            pytest.raises(ValueError, lambda: h[1].data['F1'])

        try:
            with fits.open(self.temp('test.fits')) as h:
                h[1].data['F1']
        except ValueError as e:
            assert str(e).endswith(
                         "the header may be missing the necessary TNULL1 "
                         "keyword or the table contains invalid data")

    def test_blank_field_zero(self):
        """Regression test for https://github.com/astropy/astropy/issues/5134

        Blank values in numerical columns of ASCII tables should be replaced
        with zeros, so they can be loaded into numpy arrays.

        When a TNULL value is set and there are blank fields not equal to that
        value, they should be replaced with zeros.
        """

        # Test an integer column with blank string as null
        nullval1 = ' '

        c1 = fits.Column('F1', format='I8', null=nullval1,
                         array=np.array([0, 1, 2, 3, 4]),
                         ascii=True)
        table = fits.TableHDU.from_columns([c1])
        table.writeto(self.temp('ascii_null.fits'))

        # Replace the 1st col, 3rd row, with a null field.
        with open(self.temp('ascii_null.fits'), mode='r+') as h:
            nulled = h.read().replace('2       ', '        ')
            h.seek(0)
            h.write(nulled)

        with fits.open(self.temp('ascii_null.fits'), memmap=True) as f:
            assert f[1].data[2][0] == 0

        # Test a float column with a null value set and blank fields.
        nullval2 = 'NaN'
        c2 = fits.Column('F1', format='F12.8', null=nullval2,
                         array=np.array([1.0, 2.0, 3.0, 4.0]),
                         ascii=True)
        table = fits.TableHDU.from_columns([c2])
        table.writeto(self.temp('ascii_null2.fits'))

        # Replace the 1st col, 3rd row, with a null field.
        with open(self.temp('ascii_null2.fits'), mode='r+') as h:
            nulled = h.read().replace('3.00000000', '          ')
            h.seek(0)
            h.write(nulled)

        with fits.open(self.temp('ascii_null2.fits'), memmap=True) as f:
            # (Currently it should evaluate to 0.0, but if a TODO in fitsrec is
            # completed, then it should evaluate to NaN.)
            assert f[1].data[2][0] == 0.0 or np.isnan(f[1].data[2][0])

    def test_column_array_type_mismatch(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/218"""

        arr = [-99] * 20
        col = fits.Column('mag', format='E', array=arr)
        assert (arr == col.array).all()

    def test_table_none(self):
        """Regression test
        for https://github.com/spacetelescope/PyFITS/issues/27
        """

        with fits.open(self.data('tb.fits')) as h:
            h[1].data
            h[1].data = None
            assert isinstance(h[1].data, fits.FITS_rec)
            assert len(h[1].data) == 0
            h[1].writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            assert h[1].header['NAXIS'] == 2
            assert h[1].header['NAXIS1'] == 12
            assert h[1].header['NAXIS2'] == 0
            assert isinstance(h[1].data, fits.FITS_rec)
            assert len(h[1].data) == 0

    def test_unncessary_table_load(self):
        """Test unnecessary parsing and processing of FITS tables when writing
        directly from one FITS file to a new file without first reading the
        data for user manipulation.

        In other words, it should be possible to do a direct copy of the raw
        data without unnecessary processing of the data.
        """

        with fits.open(self.data('table.fits')) as h:
            h[1].writeto(self.temp('test.fits'))

        # Since this was a direct copy the h[1].data attribute should not have
        # even been accessed (since this means the data was read and parsed)
        assert 'data' not in h[1].__dict__

        with fits.open(self.data('table.fits')) as h1:
            with fits.open(self.temp('test.fits')) as h2:
                assert str(h1[1].header) == str(h2[1].header)
                assert comparerecords(h1[1].data, h2[1].data)

    def test_table_from_columns_of_other_table(self):
        """Tests a rare corner case where the columns of an existing table
        are used to create a new table with the new_table function.  In this
        specific case, however, the existing table's data has not been read
        yet, so new_table has to get at it through the Delayed proxy.

        Note: Although this previously tested new_table it now uses
        BinTableHDU.from_columns directly, around which new_table is a mere
        wrapper.
        """

        hdul = fits.open(self.data('table.fits'))

        # Make sure the column array is in fact delayed...
        assert isinstance(hdul[1].columns._arrays[0], Delayed)

        # Create a new table...
        t = fits.BinTableHDU.from_columns(hdul[1].columns)

        # The original columns should no longer be delayed...
        assert not isinstance(hdul[1].columns._arrays[0], Delayed)

        t.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul2:
            assert comparerecords(hdul[1].data, hdul2[1].data)

        hdul.close()

    def test_bintable_to_asciitable(self):
        """Tests initializing a TableHDU with the data from a BinTableHDU."""

        with fits.open(self.data('tb.fits')) as hdul:
            tbdata = hdul[1].data
            tbhdu = fits.TableHDU(data=tbdata)
            tbhdu.writeto(self.temp('test.fits'), overwrite=True)
            with fits.open(self.temp('test.fits')) as hdul2:
                tbdata2 = hdul2[1].data
                assert np.all(tbdata['c1'] == tbdata2['c1'])
                assert np.all(tbdata['c2'] == tbdata2['c2'])
                # c3 gets converted from float32 to float64 when writing
                # test.fits, so cast to float32 before testing that the correct
                # value is retrieved
                assert np.all(tbdata['c3'].astype(np.float32) ==
                              tbdata2['c3'].astype(np.float32))
                # c4 is a boolean column in the original table; we want ASCII
                # columns to convert these to columns of 'T'/'F' strings
                assert np.all(np.where(tbdata['c4'], 'T', 'F') ==
                              tbdata2['c4'])

    def test_pickle(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/1597

        Tests for pickling FITS_rec objects
        """

        # open existing FITS tables (images pickle by default, no test needed):
        with fits.open(self.data('tb.fits')) as btb:
            # Test column array is delayed and can pickle
            assert isinstance(btb[1].columns._arrays[0], Delayed)

            btb_pd = pickle.dumps(btb[1].data)
            btb_pl = pickle.loads(btb_pd)

            # It should not be delayed any more
            assert not isinstance(btb[1].columns._arrays[0], Delayed)

            assert comparerecords(btb_pl, btb[1].data)

        with fits.open(self.data('ascii.fits')) as asc:
            asc_pd = pickle.dumps(asc[1].data)
            asc_pl = pickle.loads(asc_pd)
            assert comparerecords(asc_pl, asc[1].data)

        with fits.open(self.data('random_groups.fits')) as rgr:
            rgr_pd = pickle.dumps(rgr[0].data)
            rgr_pl = pickle.loads(rgr_pd)
            assert comparerecords(rgr_pl, rgr[0].data)

        with fits.open(self.data('zerowidth.fits')) as zwc:
            # Doesn't pickle zero-width (_phanotm) column 'ORBPARM'
            zwc_pd = pickle.dumps(zwc[2].data)
            zwc_pl = pickle.loads(zwc_pd)
            with pytest.warns(UserWarning, match='Field 2 has a repeat count of 0'):
                assert comparerecords(zwc_pl, zwc[2].data)

    def test_zero_length_table(self):
        array = np.array([], dtype=[
            ('a', 'i8'),
            ('b', 'S64'),
            ('c', ('i4', (3, 2)))])
        hdu = fits.BinTableHDU(array)
        assert hdu.header['NAXIS1'] == 96
        assert hdu.header['NAXIS2'] == 0
        assert hdu.header['TDIM3'] == '(2,3)'

        field = hdu.data.field(1)
        assert field.shape == (0,)

    def test_dim_column_byte_order_mismatch(self):
        """
        When creating a table column with non-trivial TDIMn, and
        big-endian array data read from an existing FITS file, the data
        should not be unnecessarily byteswapped.

        Regression test for https://github.com/astropy/astropy/issues/3561
        """

        data = fits.getdata(self.data('random_groups.fits'))['DATA']
        col = fits.Column(name='TEST', array=data, dim='(3,1,128,1,1)',
                          format='1152E')
        thdu = fits.BinTableHDU.from_columns([col])
        thdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            assert np.all(hdul[1].data['TEST'] == data)

    def test_fits_rec_from_existing(self):
        """
        Tests creating a `FITS_rec` object with `FITS_rec.from_columns`
        from an existing `FITS_rec` object read from a FITS file.

        This ensures that the per-column arrays are updated properly.

        Regression test for https://github.com/spacetelescope/PyFITS/issues/99
        """

        # The use case that revealed this problem was trying to create a new
        # table from an existing table, but with additional rows so that we can
        # append data from a second table (with the same column structure)

        data1 = fits.getdata(self.data('tb.fits'))
        data2 = fits.getdata(self.data('tb.fits'))
        nrows = len(data1) + len(data2)

        merged = fits.FITS_rec.from_columns(data1, nrows=nrows)
        merged[len(data1):] = data2
        mask = merged['c1'] > 1
        masked = merged[mask]

        # The test table only has two rows, only the second of which is > 1 for
        # the 'c1' column
        assert comparerecords(data1[1:], masked[:1])
        assert comparerecords(data1[1:], masked[1:])

        # Double check that the original data1 table hasn't been affected by
        # its use in creating the "merged" table
        assert comparerecords(data1, fits.getdata(self.data('tb.fits')))

    def test_update_string_column_inplace(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/4452

        Ensure that changes to values in a string column are saved when
        a file is opened in ``mode='update'``.
        """

        data = np.array([('abc',)], dtype=[('a', 'S3')])
        fits.writeto(self.temp('test.fits'), data)

        with fits.open(self.temp('test.fits'), mode='update') as hdul:
            hdul[1].data['a'][0] = 'XYZ'
            assert hdul[1].data['a'][0] == 'XYZ'

        with fits.open(self.temp('test.fits')) as hdul:
            assert hdul[1].data['a'][0] == 'XYZ'

        # Test update but with a non-trivial TDIMn
        data = np.array([([['abc', 'def', 'geh'],
                           ['ijk', 'lmn', 'opq']],)],
                        dtype=[('a', ('S3', (2, 3)))])

        fits.writeto(self.temp('test2.fits'), data)

        expected = [['abc', 'def', 'geh'],
                    ['ijk', 'XYZ', 'opq']]

        with fits.open(self.temp('test2.fits'), mode='update') as hdul:
            assert hdul[1].header['TDIM1'] == '(3,3,2)'
            # Note: Previously I wrote data['a'][0][1, 1] to address
            # the single row.  However, this is broken for chararray because
            # data['a'][0] does *not* return a view of the original array--this
            # is a bug in chararray though and not a bug in any FITS-specific
            # code so we'll roll with it for now...
            # (by the way the bug in question is fixed in newer Numpy versions)
            hdul[1].data['a'][0, 1, 1] = 'XYZ'
            assert np.all(hdul[1].data['a'][0] == expected)

        with fits.open(self.temp('test2.fits')) as hdul:
            assert hdul[1].header['TDIM1'] == '(3,3,2)'
            assert np.all(hdul[1].data['a'][0] == expected)

    @pytest.mark.skipif('not HAVE_OBJGRAPH')
    def test_reference_leak(self):
        """Regression test for https://github.com/astropy/astropy/pull/520"""

        def readfile(filename):
            with fits.open(filename) as hdul:
                data = hdul[1].data.copy()

            for colname in data.dtype.names:
                data[colname]

        with _refcounting('FITS_rec'):
            readfile(self.data('memtest.fits'))

    @pytest.mark.skipif('not HAVE_OBJGRAPH')
    @pytest.mark.slow
    def test_reference_leak2(self, tmpdir):
        """
        Regression test for https://github.com/astropy/astropy/pull/4539

        This actually re-runs a small set of tests that I found, during
        careful testing, exhibited the reference leaks fixed by #4539, but
        now with reference counting around each test to ensure that the
        leaks are fixed.
        """

        from .test_core import TestCore
        from .test_connect import TestMultipleHDU

        t1 = TestCore()
        t1.setup()
        try:
            with _refcounting('FITS_rec'):
                t1.test_add_del_columns2()
        finally:
            t1.teardown()
        del t1

        t2 = self.__class__()
        for test_name in ['test_recarray_to_bintablehdu',
                          'test_numpy_ndarray_to_bintablehdu',
                          'test_new_table_from_recarray',
                          'test_new_fitsrec']:
            t2.setup()
            try:
                with _refcounting('FITS_rec'):
                    getattr(t2, test_name)()
            finally:
                t2.teardown()
        del t2

        t3 = TestMultipleHDU()
        t3.setup_class()
        try:
            with _refcounting('FITS_rec'):
                t3.test_read(tmpdir)
        finally:
            t3.teardown_class()
        del t3

    def test_dump_overwrite(self):
        with fits.open(self.data('table.fits')) as hdul:
            tbhdu = hdul[1]
            datafile = self.temp('data.txt')
            cdfile = self.temp('coldefs.txt')
            hfile = self.temp('header.txt')
            tbhdu.dump(datafile, cdfile, hfile)
            msg = (r"File .* already exists\.  File .* already exists\.  File "
                   r".* already exists\.  If you mean to replace the "
                   r"file\(s\) then use the argument 'overwrite=True'\.")
            with pytest.raises(OSError, match=msg):
                tbhdu.dump(datafile, cdfile, hfile)
            tbhdu.dump(datafile, cdfile, hfile, overwrite=True)

    def test_pseudo_unsigned_ints(self):
        """
        Tests updating a table column containing pseudo-unsigned ints.
        """

        data = np.array([1, 2, 3], dtype=np.uint32)
        col = fits.Column(name='A', format='1J', bzero=2**31, array=data)
        thdu = fits.BinTableHDU.from_columns([col])
        thdu.writeto(self.temp('test.fits'))

        # Test that the file wrote out correctly
        with fits.open(self.temp('test.fits'), uint=True) as hdul:
            hdu = hdul[1]
            assert 'TZERO1' in hdu.header
            assert hdu.header['TZERO1'] == 2**31
            assert hdu.data['A'].dtype == np.dtype('uint32')
            assert np.all(hdu.data['A'] == data)

            # Test updating the unsigned int data
            hdu.data['A'][0] = 99
            hdu.writeto(self.temp('test2.fits'))

        with fits.open(self.temp('test2.fits'), uint=True) as hdul:
            hdu = hdul[1]
            assert 'TZERO1' in hdu.header
            assert hdu.header['TZERO1'] == 2**31
            assert hdu.data['A'].dtype == np.dtype('uint32')
            assert np.all(hdu.data['A'] == [99, 2, 3])

    def test_column_with_scaling(self):
        """Check that a scaled column if correctly saved once it is modified.
        Regression test for https://github.com/astropy/astropy/issues/6887
        """
        c1 = fits.Column(name='c1', array=np.array([1], dtype='>i2'),
                         format='1I', bscale=1, bzero=32768)
        S = fits.HDUList([fits.PrimaryHDU(),
                          fits.BinTableHDU.from_columns([c1])])

        # Change value in memory
        S[1].data['c1'][0] = 2
        S.writeto(self.temp("a.fits"))
        assert S[1].data['c1'] == 2

        # Read and change value in memory
        with fits.open(self.temp("a.fits")) as X:
            X[1].data['c1'][0] = 10
            assert X[1].data['c1'][0] == 10

            # Write back to file
            X.writeto(self.temp("b.fits"))

        # Now check the file
        with fits.open(self.temp("b.fits")) as hdul:
            assert hdul[1].data['c1'][0] == 10

    def test_ascii_inttypes(self):
        """
        Test correct integer dtypes according to ASCII table field widths.
        Regression for https://github.com/astropy/astropy/issues/9899
        """
        i08 = np.array([2**3, 2**23, -2**22, 10, 2**23], dtype='i4')
        i10 = np.array([2**8, 2**31-1, -2**29, 30, 2**31-1], dtype='i8')
        i20 = np.array([2**16, 2**63-1, -2**63, 40, 2**63-1], dtype='i8')
        i02 = np.array([2**8, 2**13, -2**9, 50, 2**13], dtype='i2')
        t0 = Table([i08, i08*2, i10, i20, i02])

        t1 = Table.read(self.data('ascii_i4-i20.fits'))
        assert t1.dtype == t0.dtype
        assert comparerecords(t1, t0)


@contextlib.contextmanager
def _refcounting(type_):
    """
    Perform the body of a with statement with reference counting for the
    given type (given by class name)--raises an assertion error if there
    are more unfreed objects of the given type than when we entered the
    with statement.
    """

    gc.collect()
    refcount = len(objgraph.by_type(type_))
    yield refcount
    gc.collect()
    assert len(objgraph.by_type(type_)) <= refcount, \
        "More {0!r} objects still in memory than before."


class TestVLATables(FitsTestCase):
    """Tests specific to tables containing variable-length arrays."""

    def test_variable_length_columns(self):
        def test(format_code):
            col = fits.Column(name='QUAL_SPE', format=format_code,
                              array=[[0] * 1571] * 225)
            tb_hdu = fits.BinTableHDU.from_columns([col])
            pri_hdu = fits.PrimaryHDU()
            hdu_list = fits.HDUList([pri_hdu, tb_hdu])
            hdu_list.writeto(self.temp('toto.fits'), overwrite=True)

            with fits.open(self.temp('toto.fits')) as toto:
                q = toto[1].data.field('QUAL_SPE')
                assert (q[0][4:8] ==
                        np.array([0, 0, 0, 0], dtype=np.uint8)).all()
                assert toto[1].columns[0].format.endswith('J(1571)')

        for code in ('PJ()', 'QJ()'):
            test(code)

    def test_extend_variable_length_array(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/54"""

        def test(format_code):
            arr = [[1] * 10] * 10
            col1 = fits.Column(name='TESTVLF', format=format_code, array=arr)
            col2 = fits.Column(name='TESTSCA', format='J', array=[1] * 10)
            tb_hdu = fits.BinTableHDU.from_columns([col1, col2], nrows=15)
            # This asserts that the normal 'scalar' column's length was extended
            assert len(tb_hdu.data['TESTSCA']) == 15
            # And this asserts that the VLF column was extended in the same manner
            assert len(tb_hdu.data['TESTVLF']) == 15
            # We can't compare the whole array since the _VLF is an array of
            # objects, but comparing just the edge case rows should suffice
            assert (tb_hdu.data['TESTVLF'][0] == arr[0]).all()
            assert (tb_hdu.data['TESTVLF'][9] == arr[9]).all()
            assert (tb_hdu.data['TESTVLF'][10] == ([0] * 10)).all()
            assert (tb_hdu.data['TESTVLF'][-1] == ([0] * 10)).all()

        for code in ('PJ()', 'QJ()'):
            test(code)

    def test_variable_length_table_format_pd_from_object_array(self):
        def test(format_code):
            a = np.array([np.array([7.2e-20, 7.3e-20]), np.array([0.0]),
                          np.array([0.0])], 'O')
            acol = fits.Column(name='testa', format=format_code, array=a)
            tbhdu = fits.BinTableHDU.from_columns([acol])
            tbhdu.writeto(self.temp('newtable.fits'), overwrite=True)
            with fits.open(self.temp('newtable.fits')) as tbhdu1:
                assert tbhdu1[1].columns[0].format.endswith('D(2)')
                for j in range(3):
                    for i in range(len(a[j])):
                        assert tbhdu1[1].data.field(0)[j][i] == a[j][i]

        for code in ('PD()', 'QD()'):
            test(code)

    def test_variable_length_table_format_pd_from_list(self):
        def test(format_code):
            a = [np.array([7.2e-20, 7.3e-20]), np.array([0.0]),
                 np.array([0.0])]
            acol = fits.Column(name='testa', format=format_code, array=a)
            tbhdu = fits.BinTableHDU.from_columns([acol])
            tbhdu.writeto(self.temp('newtable.fits'), overwrite=True)

            with fits.open(self.temp('newtable.fits')) as tbhdu1:
                assert tbhdu1[1].columns[0].format.endswith('D(2)')
                for j in range(3):
                    for i in range(len(a[j])):
                        assert tbhdu1[1].data.field(0)[j][i] == a[j][i]

        for code in ('PD()', 'QD()'):
            test(code)

    def test_variable_length_table_format_pa_from_object_array(self):
        def test(format_code):
            a = np.array([np.array(['a', 'b', 'c']), np.array(['d', 'e']),
                          np.array(['f'])], 'O')
            acol = fits.Column(name='testa', format=format_code, array=a)
            tbhdu = fits.BinTableHDU.from_columns([acol])
            tbhdu.writeto(self.temp('newtable.fits'), overwrite=True)

            with fits.open(self.temp('newtable.fits')) as hdul:
                assert hdul[1].columns[0].format.endswith('A(3)')
                for j in range(3):
                    for i in range(len(a[j])):
                        assert hdul[1].data.field(0)[j][i] == a[j][i]

        for code in ('PA()', 'QA()'):
            test(code)

    def test_variable_length_table_format_pa_from_list(self):
        def test(format_code):
            a = ['a', 'ab', 'abc']
            acol = fits.Column(name='testa', format=format_code, array=a)
            tbhdu = fits.BinTableHDU.from_columns([acol])
            tbhdu.writeto(self.temp('newtable.fits'), overwrite=True)

            with fits.open(self.temp('newtable.fits')) as hdul:
                assert hdul[1].columns[0].format.endswith('A(3)')
                for j in range(3):
                    for i in range(len(a[j])):
                        assert hdul[1].data.field(0)[j][i] == a[j][i]

        for code in ('PA()', 'QA()'):
            test(code)

    def test_getdata_vla(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/200"""

        def test(format_code):
            col = fits.Column(name='QUAL_SPE', format=format_code,
                              array=[np.arange(1572)] * 225)
            tb_hdu = fits.BinTableHDU.from_columns([col])
            pri_hdu = fits.PrimaryHDU()
            hdu_list = fits.HDUList([pri_hdu, tb_hdu])
            hdu_list.writeto(self.temp('toto.fits'), overwrite=True)

            data = fits.getdata(self.temp('toto.fits'))

            # Need to compare to the original data row by row since the FITS_rec
            # returns an array of _VLA objects
            for row_a, row_b in zip(data['QUAL_SPE'], col.array):
                assert (row_a == row_b).all()

        for code in ('PJ()', 'QJ()'):
            test(code)

    @pytest.mark.skipif(not NUMPY_LT_1_22 and NUMPY_LT_1_22_1 and sys.platform == 'win32',
                        reason='https://github.com/numpy/numpy/issues/20699')
    def test_copy_vla(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/47
        """

        # Make a file containing a couple of VLA tables
        arr1 = [np.arange(n + 1) for n in range(255)]
        arr2 = [np.arange(255, 256 + n) for n in range(255)]

        # A dummy non-VLA column needed to reproduce issue #47
        c = fits.Column('test', format='J', array=np.arange(255))
        c1 = fits.Column('A', format='PJ', array=arr1)
        c2 = fits.Column('B', format='PJ', array=arr2)
        t1 = fits.BinTableHDU.from_columns([c, c1])
        t2 = fits.BinTableHDU.from_columns([c, c2])

        hdul = fits.HDUList([fits.PrimaryHDU(), t1, t2])
        hdul.writeto(self.temp('test.fits'), overwrite=True)

        # Just test that the test file wrote out correctly
        with fits.open(self.temp('test.fits')) as h:
            assert h[1].header['TFORM2'] == 'PJ(255)'
            assert h[2].header['TFORM2'] == 'PJ(255)'
            assert comparerecords(h[1].data, t1.data)
            assert comparerecords(h[2].data, t2.data)

        # Try copying the second VLA and writing to a new file
        with fits.open(self.temp('test.fits')) as h:
            new_hdu = fits.BinTableHDU(data=h[2].data, header=h[2].header)
            new_hdu.writeto(self.temp('test3.fits'))

        with fits.open(self.temp('test3.fits')) as h2:
            assert comparerecords(h2[1].data, t2.data)

        new_hdul = fits.HDUList([fits.PrimaryHDU()])
        new_hdul.writeto(self.temp('test2.fits'))

        # Open several copies of the test file and append copies of the second
        # VLA table
        with fits.open(self.temp('test2.fits'), mode='append') as new_hdul:
            for _ in range(2):
                with fits.open(self.temp('test.fits')) as h:
                    new_hdul.append(h[2])
                    new_hdul.flush()

        # Test that all the VLA copies wrote correctly
        with fits.open(self.temp('test2.fits')) as new_hdul:
            for idx in range(1, 3):
                assert comparerecords(new_hdul[idx].data, t2.data)

    def test_vla_with_gap(self):
        hdul = fits.open(self.data('theap-gap.fits'))
        data = hdul[1].data
        assert data.shape == (500,)
        assert data['i'][497] == 497
        assert np.array_equal(data['arr'][497], [0, 1, 2, 3, 4])
        hdul.close()

    def test_tolist(self):
        col = fits.Column(
            name='var', format='PI()',
            array=np.array([[1, 2, 3], [11, 12]], dtype=np.object_))
        hdu = fits.BinTableHDU.from_columns([col])
        assert hdu.data.tolist() == [[[1, 2, 3]], [[11, 12]]]
        assert hdu.data['var'].tolist() == [[1, 2, 3], [11, 12]]

    def test_tolist_from_file(self):
        filename = self.data('variable_length_table.fits')

        with fits.open(filename) as hdul:
            hdu = hdul[1]
            assert hdu.data.tolist() == [[[45, 56], [11, 3]], [[11, 12, 13], [12, 4]]]
            assert hdu.data['var'].tolist() == [[45, 56], [11, 12, 13]]


# These are tests that solely test the Column and ColDefs interfaces and
# related functionality without directly involving full tables; currently there
# are few of these but I expect there to be more as I improve the test coverage
class TestColumnFunctions(FitsTestCase):
    def test_column_format_interpretation(self):
        """
        Test to ensure that when Numpy-style record formats are passed in to
        the Column constructor for the format argument, they are recognized so
        long as it's unambiguous (where "unambiguous" here is questionable
        since Numpy is case insensitive when parsing the format codes.  But
        their "proper" case is lower-case, so we can accept that.  Basically,
        actually, any key in the NUMPY2FITS dict should be accepted.
        """

        for recformat, fitsformat in NUMPY2FITS.items():
            c = fits.Column('TEST', np.dtype(recformat))
            c.format == fitsformat
            c = fits.Column('TEST', recformat)
            c.format == fitsformat
            c = fits.Column('TEST', fitsformat)
            c.format == fitsformat

        # Test a few cases that are ambiguous in that they *are* valid binary
        # table formats though not ones that are likely to be used, but are
        # also valid common ASCII table formats
        c = fits.Column('TEST', 'I4')
        assert c.format == 'I4'
        assert c.format.format == 'I'
        assert c.format.width == 4

        c = fits.Column('TEST', 'F15.8')
        assert c.format == 'F15.8'
        assert c.format.format == 'F'
        assert c.format.width == 15
        assert c.format.precision == 8

        c = fits.Column('TEST', 'E15.8')
        assert c.format.format == 'E'
        assert c.format.width == 15
        assert c.format.precision == 8

        c = fits.Column('TEST', 'D15.8')
        assert c.format.format == 'D'
        assert c.format.width == 15
        assert c.format.precision == 8

        # zero-precision should be allowed as well, for float types
        # https://github.com/astropy/astropy/issues/3422
        c = fits.Column('TEST', 'F10.0')
        assert c.format.format == 'F'
        assert c.format.width == 10
        assert c.format.precision == 0

        c = fits.Column('TEST', 'E10.0')
        assert c.format.format == 'E'
        assert c.format.width == 10
        assert c.format.precision == 0

        c = fits.Column('TEST', 'D10.0')
        assert c.format.format == 'D'
        assert c.format.width == 10
        assert c.format.precision == 0

        # These are a couple cases where the format code is a valid binary
        # table format, and is not strictly a valid ASCII table format but
        # could be *interpreted* as one by appending a default width.  This
        # will only happen either when creating an ASCII table or when
        # explicitly specifying ascii=True when the column is created
        c = fits.Column('TEST', 'I')
        assert c.format == 'I'
        assert c.format.recformat == 'i2'
        c = fits.Column('TEST', 'I', ascii=True)
        assert c.format == 'I10'
        assert c.format.recformat == 'i4'

        # With specified widths, integer precision should be set appropriately
        c = fits.Column('TEST', 'I4', ascii=True)
        assert c.format == 'I4'
        assert c.format.recformat == 'i2'
        c = fits.Column('TEST', 'I9', ascii=True)
        assert c.format == 'I9'
        assert c.format.recformat == 'i4'
        c = fits.Column('TEST', 'I12', ascii=True)
        assert c.format == 'I12'
        assert c.format.recformat == 'i8'

        c = fits.Column('TEST', 'E')
        assert c.format == 'E'
        assert c.format.recformat == 'f4'
        c = fits.Column('TEST', 'E', ascii=True)
        assert c.format == 'E15.7'

        # F is not a valid binary table format so it should be unambiguously
        # treated as an ASCII column
        c = fits.Column('TEST', 'F')
        assert c.format == 'F16.7'

        c = fits.Column('TEST', 'D')
        assert c.format == 'D'
        assert c.format.recformat == 'f8'
        c = fits.Column('TEST', 'D', ascii=True)
        assert c.format == 'D25.17'

    def test_zero_precision_float_column(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/3422
        """

        c = fits.Column('TEST', 'F5.0', array=[1.1, 2.2, 3.3])
        # The decimal places will be clipped
        t = fits.TableHDU.from_columns([c])
        t.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            assert hdul[1].header['TFORM1'] == 'F5.0'
            assert hdul[1].data['TEST'].dtype == np.dtype('float64')
            assert np.all(hdul[1].data['TEST'] == [1.0, 2.0, 3.0])

            # Check how the raw data looks
            raw = np.rec.recarray.field(hdul[1].data, 'TEST')
            assert raw.tobytes() == b'   1.   2.   3.'

    def test_column_array_type_mismatch(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/218"""

        arr = [-99] * 20
        col = fits.Column('mag', format='E', array=arr)
        assert (arr == col.array).all()

    def test_new_coldefs_with_invalid_seqence(self):
        """Test that a TypeError is raised when a ColDefs is instantiated with
        a sequence of non-Column objects.
        """

        pytest.raises(TypeError, fits.ColDefs, [1, 2, 3])

    def test_coldefs_init_from_array(self):
        """Test that ColDefs._init_from_array works with single element data-
        types as well as multi-element data-types
        """
        nd_array = np.ndarray((1,), dtype=[('A', '<u4', (2,)), ('B', '>u2')])
        col_defs = fits.column.ColDefs(nd_array)
        assert 2**31 == col_defs['A'].bzero
        assert 2**15 == col_defs['B'].bzero

    def test_pickle(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/1597

        Tests for pickling FITS_rec objects
        """

        # open existing FITS tables (images pickle by default, no test needed):
        with fits.open(self.data('tb.fits')) as btb:
            # Test column array is delayed and can pickle
            assert isinstance(btb[1].columns._arrays[0], Delayed)

            btb_pd = pickle.dumps(btb[1].data)
            btb_pl = pickle.loads(btb_pd)

            # It should not be delayed any more
            assert not isinstance(btb[1].columns._arrays[0], Delayed)

            assert comparerecords(btb_pl, btb[1].data)

        with fits.open(self.data('ascii.fits')) as asc:
            asc_pd = pickle.dumps(asc[1].data)
            asc_pl = pickle.loads(asc_pd)
            assert comparerecords(asc_pl, asc[1].data)

        with fits.open(self.data('random_groups.fits')) as rgr:
            rgr_pd = pickle.dumps(rgr[0].data)
            rgr_pl = pickle.loads(rgr_pd)
            assert comparerecords(rgr_pl, rgr[0].data)

        with fits.open(self.data('zerowidth.fits')) as zwc:
            # Doesn't pickle zero-width (_phanotm) column 'ORBPARM'
            zwc_pd = pickle.dumps(zwc[2].data)
            zwc_pl = pickle.loads(zwc_pd)
            with pytest.warns(UserWarning, match=r'Field 2 has a repeat count '
                              r'of 0 in its format code'):
                assert comparerecords(zwc_pl, zwc[2].data)

    def test_column_lookup_by_name(self):
        """Tests that a `ColDefs` can be indexed by column name."""

        a = fits.Column(name='a', format='D')
        b = fits.Column(name='b', format='D')

        cols = fits.ColDefs([a, b])

        assert cols['a'] == cols[0]
        assert cols['b'] == cols[1]

    def test_column_attribute_change_after_removal(self):
        """
        This is a test of the column attribute change notification system.

        After a column has been removed from a table (but other references
        are kept to that same column) changes to that column's attributes
        should not trigger a notification on the table it was removed from.
        """

        # One way we can check this is to ensure there are no further changes
        # to the header
        table = fits.BinTableHDU.from_columns([
            fits.Column('a', format='D'),
            fits.Column('b', format='D')])

        b = table.columns['b']

        table.columns.del_col('b')
        assert table.data.dtype.names == ('a',)

        b.name = 'HELLO'

        assert b.name == 'HELLO'
        assert 'TTYPE2' not in table.header
        assert table.header['TTYPE1'] == 'a'
        assert table.columns.names == ['a']

        with pytest.raises(KeyError):
            table.columns['b']

        # Make sure updates to the remaining column still work
        table.columns.change_name('a', 'GOODBYE')
        with pytest.raises(KeyError):
            table.columns['a']

        assert table.columns['GOODBYE'].name == 'GOODBYE'
        assert table.data.dtype.names == ('GOODBYE',)
        assert table.columns.names == ['GOODBYE']
        assert table.data.columns.names == ['GOODBYE']

        table.columns['GOODBYE'].name = 'foo'
        with pytest.raises(KeyError):
            table.columns['GOODBYE']

        assert table.columns['foo'].name == 'foo'
        assert table.data.dtype.names == ('foo',)
        assert table.columns.names == ['foo']
        assert table.data.columns.names == ['foo']

    def test_x_column_deepcopy(self):
        """
        Regression test for https://github.com/astropy/astropy/pull/4514

        Tests that columns with the X (bit array) format can be deep-copied.
        """

        c = fits.Column('xcol', format='5X', array=[1, 0, 0, 1, 0])
        c2 = copy.deepcopy(c)
        assert c2.name == c.name
        assert c2.format == c.format
        assert np.all(c2.array == c.array)

    def test_p_column_deepcopy(self):
        """
        Regression test for https://github.com/astropy/astropy/pull/4514

        Tests that columns with the P/Q formats (variable length arrays) can be
        deep-copied.
        """

        c = fits.Column('pcol', format='PJ', array=[[1, 2], [3, 4, 5]])
        c2 = copy.deepcopy(c)
        assert c2.name == c.name
        assert c2.format == c.format
        assert np.all(c2.array[0] == c.array[0])
        assert np.all(c2.array[1] == c.array[1])

        c3 = fits.Column('qcol', format='QJ', array=[[1, 2], [3, 4, 5]])
        c4 = copy.deepcopy(c3)
        assert c4.name == c3.name
        assert c4.format == c3.format
        assert np.all(c4.array[0] == c3.array[0])
        assert np.all(c4.array[1] == c3.array[1])

    def test_column_verify_keywords(self):
        """
        Test that the keyword arguments used to initialize a Column, specifically
        those that typically read from a FITS header (so excluding array),
        are verified to have a valid value.
        """

        with pytest.raises(AssertionError) as err:
            _ = fits.Column(1, format='I', array=[1, 2, 3, 4, 5])
        assert 'Column name must be a string able to fit' in str(err.value)

        with pytest.raises(VerifyError) as err:
            _ = fits.Column('col', format=0, null='Nan', disp=1, coord_type=1,
                            coord_unit=2, coord_inc='1', time_ref_pos=1,
                            coord_ref_point='1', coord_ref_value='1')
        err_msgs = ['keyword arguments to Column were invalid',
                    'TFORM', 'TNULL', 'TDISP', 'TCTYP', 'TCUNI', 'TCRPX',
                    'TCRVL', 'TCDLT', 'TRPOS']
        for msg in err_msgs:
            assert msg in str(err.value)

    def test_column_verify_start(self):
        """
        Regression test for https://github.com/astropy/astropy/pull/6359

        Test the validation of the column start position option (ASCII table only),
        corresponding to ``TBCOL`` keyword.
        Test whether the VerifyError message generated is the one with highest priority,
        i.e. the order of error messages to be displayed is maintained.
        """

        with pytest.raises(VerifyError) as err:
            _ = fits.Column('a', format='B', start='a', array=[1, 2, 3])
        assert "start option (TBCOLn) is not allowed for binary table columns" in str(err.value)

        with pytest.raises(VerifyError) as err:
            _ = fits.Column('a', format='I', start='a', array=[1, 2, 3])
        assert "start option (TBCOLn) must be a positive integer (got 'a')." in str(err.value)

        with pytest.raises(VerifyError) as err:
            _ = fits.Column('a', format='I', start='-56', array=[1, 2, 3])
        assert "start option (TBCOLn) must be a positive integer (got -56)." in str(err.value)

    @pytest.mark.parametrize('keys',
                             [{'TFORM': 'Z', 'TDISP': 'E'},
                              {'TFORM': '2', 'TDISP': '2E'},
                              {'TFORM': 3, 'TDISP': 6.3},
                              {'TFORM': float, 'TDISP': np.float64},
                              {'TFORM': '', 'TDISP': 'E.5'}])
    def test_column_verify_formats(self, keys):
        """
        Additional tests for verification of 'TFORM' and 'TDISP' keyword
        arguments used to initialize a Column.
        """
        with pytest.raises(VerifyError) as err:
            _ = fits.Column('col', format=keys['TFORM'], disp=keys['TDISP'])

        for key in keys.keys():
            assert key in str(err.value)
            assert str(keys[key]) in str(err.value)


def test_regression_5383():

    # Regression test for an undefined variable

    x = np.array([1, 2, 3])
    col = fits.Column(name='a', array=x, format='E')
    hdu = fits.BinTableHDU.from_columns([col])
    del hdu._header['TTYPE1']
    hdu.columns[0].name = 'b'


def test_table_to_hdu():
    from astropy.table import Table
    table = Table([[1, 2, 3], ['a', 'b', 'c'], [2.3, 4.5, 6.7]],
                  names=['a', 'b', 'c'], dtype=['i', 'U1', 'f'])
    table['a'].unit = 'm/s'
    table['b'].unit = 'not-a-unit'
    table.meta['foo'] = 'bar'

    with pytest.warns(UnitsWarning, match="'not-a-unit' did not parse as"
                      " fits unit") as w:
        hdu = fits.BinTableHDU(table, header=fits.Header({'TEST': 1}))
    assert len(w) == 1

    for name in 'abc':
        assert np.array_equal(table[name], hdu.data[name])

    # Check that TUNITn cards appear in the correct order
    # (https://github.com/astropy/astropy/pull/5720)
    assert hdu.header.index('TUNIT1') < hdu.header.index('TTYPE2')

    assert hdu.header['FOO'] == 'bar'
    assert hdu.header['TEST'] == 1


def test_regression_scalar_indexing():
    # Indexing a FITS_rec with a tuple that returns a scalar record
    # should work
    x = np.array([(1.0, 2), (3.0, 4)],
                 dtype=[('x', float), ('y', int)]).view(fits.FITS_rec)
    x1a = x[1]
    # this should succeed.
    x1b = x[(1,)]
    # FITS_record does not define __eq__; so test elements.
    assert all(a == b for a, b in zip(x1a, x1b))


def test_new_column_attributes_preserved(tmpdir):

    # Regression test for https://github.com/astropy/astropy/issues/7145
    # This makes sure that for now we don't clear away keywords that have
    # newly been recognized (in Astropy 3.0) as special column attributes but
    # instead just warn that we might do so in future. The new keywords are:
    # TCTYP, TCUNI, TCRPX, TCRVL, TCDLT, TRPOS

    col = []
    col.append(fits.Column(name="TIME", format="1E", unit="s"))
    col.append(fits.Column(name="RAWX", format="1I", unit="pixel"))
    col.append(fits.Column(name="RAWY", format="1I"))
    cd = fits.ColDefs(col)

    hdr = fits.Header()

    # Keywords that will get ignored in favor of these in the data
    hdr['TUNIT1'] = 'pixel'
    hdr['TUNIT2'] = 'm'
    hdr['TUNIT3'] = 'm'

    # Keywords that were added in Astropy 3.0 that should eventually be
    # ignored and set on the data instead
    hdr['TCTYP2'] = 'RA---TAN'
    hdr['TCTYP3'] = 'ANGLE'
    hdr['TCRVL2'] = -999.0
    hdr['TCRVL3'] = -999.0
    hdr['TCRPX2'] = 1.0
    hdr['TCRPX3'] = 1.0
    hdr['TALEN2'] = 16384
    hdr['TALEN3'] = 1024
    hdr['TCUNI2'] = 'angstrom'
    hdr['TCUNI3'] = 'deg'

    # Other non-relevant keywords
    hdr['RA'] = 1.5
    hdr['DEC'] = 3.0

    with pytest.warns(AstropyDeprecationWarning) as warning_list:
        hdu = fits.BinTableHDU.from_columns(cd, hdr)
    assert str(warning_list[0].message).startswith(
        "The following keywords are now recognized as special")

    # First, check that special keywords such as TUNIT are ignored in the header
    # We may want to change that behavior in future, but this is the way it's
    # been for a while now.

    assert hdu.columns[0].unit == 's'
    assert hdu.columns[1].unit == 'pixel'
    assert hdu.columns[2].unit is None

    assert hdu.header['TUNIT1'] == 's'
    assert hdu.header['TUNIT2'] == 'pixel'
    assert 'TUNIT3' not in hdu.header  # TUNIT3 was removed

    # Now, check that the new special keywords are actually still there
    # but weren't used to set the attributes on the data

    assert hdu.columns[0].coord_type is None
    assert hdu.columns[1].coord_type is None
    assert hdu.columns[2].coord_type is None

    assert 'TCTYP1' not in hdu.header
    assert hdu.header['TCTYP2'] == 'RA---TAN'
    assert hdu.header['TCTYP3'] == 'ANGLE'

    # Make sure that other keywords are still there

    assert hdu.header['RA'] == 1.5
    assert hdu.header['DEC'] == 3.0

    # Now we can write this HDU to a file and re-load. Re-loading *should*
    # cause the special column attribtues to be picked up (it's just that when a
    # header is manually specified, these values are ignored)

    filename = tmpdir.join('test.fits').strpath

    hdu.writeto(filename)

    # Make sure we don't emit a warning in this case
    with warnings.catch_warnings(record=True) as warning_list:
        with fits.open(filename) as hdul:
            hdu2 = hdul[1]
    assert len(warning_list) == 0

    # Check that column attributes are now correctly set

    assert hdu2.columns[0].unit == 's'
    assert hdu2.columns[1].unit == 'pixel'
    assert hdu2.columns[2].unit is None

    assert hdu2.header['TUNIT1'] == 's'
    assert hdu2.header['TUNIT2'] == 'pixel'
    assert 'TUNIT3' not in hdu2.header  # TUNIT3 was removed

    # Now, check that the new special keywords are actually still there
    # but weren't used to set the attributes on the data

    assert hdu2.columns[0].coord_type is None
    assert hdu2.columns[1].coord_type == 'RA---TAN'
    assert hdu2.columns[2].coord_type == 'ANGLE'

    assert 'TCTYP1' not in hdu2.header
    assert hdu2.header['TCTYP2'] == 'RA---TAN'
    assert hdu2.header['TCTYP3'] == 'ANGLE'

    # Make sure that other keywords are still there

    assert hdu2.header['RA'] == 1.5
    assert hdu2.header['DEC'] == 3.0


def test_empty_table(tmpdir):
    ofile = str(tmpdir.join('emptytable.fits'))
    hdu = fits.BinTableHDU(header=None, data=None, name='TEST')
    hdu.writeto(ofile)

    with fits.open(ofile) as hdul:
        assert hdul['TEST'].data.size == 0

    ofile = str(tmpdir.join('emptytable.fits.gz'))
    hdu = fits.BinTableHDU(header=None, data=None, name='TEST')
    hdu.writeto(ofile, overwrite=True)

    with fits.open(ofile) as hdul:
        assert hdul['TEST'].data.size == 0


def test_a3dtable(tmpdir):
    testfile = str(tmpdir.join('test.fits'))
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='FOO', format='J', array=np.arange(10))
    ])
    hdu.header['XTENSION'] = 'A3DTABLE'
    hdu.writeto(testfile, output_verify='ignore')

    with fits.open(testfile) as hdul:
        assert hdul[1].header['XTENSION'] == 'A3DTABLE'

        with pytest.warns(AstropyUserWarning) as w:
            hdul.verify('fix')

        assert str(w[0].message) == 'Verification reported errors:'
        assert str(w[2].message).endswith(
            'Converted the XTENSION keyword to BINTABLE.')

        assert hdul[1].header['XTENSION'] == 'BINTABLE'


def test_invalid_file(tmp_path):
    hdu = fits.BinTableHDU()
    # little trick to write an invalid card ...
    hdu.header['FOO'] = None
    hdu.header.cards['FOO']._value = np.nan

    testfile = tmp_path / 'test.fits'
    hdu.writeto(testfile, output_verify='ignore')
    with fits.open(testfile) as hdul:
        assert hdul[1].data is not None


def test_unit_parse_strict(tmp_path):
    path = tmp_path / 'invalid_unit.fits'

    # this is a unit parseable by the generic format but invalid for FITS
    invalid_unit = '1 / (MeV sr s)'
    unit = Unit(invalid_unit)

    t = Table({'a': [1, 2, 3]})
    t.write(path)
    with fits.open(path, mode='update') as hdul:
        hdul[1].header['TUNIT1'] = invalid_unit

    # default is "warn"
    with pytest.warns(UnitsWarning):
        t = Table.read(path)

    assert isinstance(t['a'].unit, UnrecognizedUnit)

    t = Table.read(path, unit_parse_strict='silent')
    assert isinstance(t['a'].unit, UnrecognizedUnit)

    with pytest.raises(ValueError):
        Table.read(path, unit_parse_strict='raise')

    with pytest.warns(UnitsWarning):
        Table.read(path, unit_parse_strict='warn')
