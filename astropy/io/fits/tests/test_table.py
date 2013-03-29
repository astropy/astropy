# Licensed under a 3-clause BSD style license - see PYFITS.rst

import numpy as np
from numpy import char as chararray

from ....io import fits
from ....tests.helper import pytest

from ..util import decode_ascii
from . import FitsTestCase
from .util import ignore_warnings


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
        print "number of fields don't match"
        return False
    for i in range(nfieldsa):
        fielda = a.field(i)
        fieldb = b.field(i)
        if fielda.dtype.char == 'S':
            fielda = decode_ascii(fielda)
        if fieldb.dtype.char == 'S':
            fieldb = decode_ascii(fieldb)
        if (type(fielda) != type(fieldb) and not
            (issubclass(type(fielda), type(fieldb)) or
             issubclass(type(fieldb), type(fielda)))):
            print "type(fielda): ", type(fielda), " fielda: ", fielda
            print "type(fieldb): ", type(fieldb), " fieldb: ", fieldb
            print 'field %d type differs' % i
            return False
        if isinstance(fielda[0], np.floating):
            if not comparefloats(fielda, fieldb):
                print "fielda: ", fielda
                print "fieldb: ", fieldb
                print 'field %d differs' % i
                return False
        elif (isinstance(fielda, fits.column._VLF) or
              isinstance(fieldb, fits.column._VLF)):
            for row in range(len(fielda)):
                if np.any(fielda[row] != fieldb[row]):
                    print 'fielda[%d]: %s' % (row, fielda[row])
                    print 'fieldb[%d]: %s' % (row, fieldb[row])
                    print 'field %d differs in row %d' (i, row)
        else:
            if np.any(fielda != fieldb):
                print "fielda: ", fielda
                print "fieldb: ", fieldb
                print 'field %d differs' % i
                return False
    return True


class TestTableFunctions(FitsTestCase):
    def test_constructor_copies_header(self):
       """
       Regression test for #153.  Ensure that a header from one HDU is copied
       when used to initialize new HDU.

       This is like the test of the same name in test_image, but tests this for
       tables as well.
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

        # create a new binary table HDU object by using the new_table function

        tbhdu = fits.new_table(x)

        # another way to create a table is by using existing table's
        # information:

        x2 = fits.ColDefs(tt[1])
        t2 = fits.new_table(x2, nrows=2)
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
        assert tbhdu.data[1][7].all() == a8[1].all()

        # and a column like this:
        assert str(tbhdu.data.field('abc')) == "['abc' 'def' 'xx']"

        # An alternative way to create a column-definitions object is from an
        # existing table.
        xx = fits.ColDefs(tt[1])

        # now we write out the newly created table HDU to a FITS file:
        fout = fits.HDUList(fits.PrimaryHDU())
        fout.append(tbhdu)
        fout.writeto(self.temp('tableout1.fits'), clobber=True)

        with fits.open(self.temp('tableout1.fits')) as f2:
            temp = f2[1].data.field(7)
            assert (temp[0] == [True, True, False, True, False, True, True,
                                True, False, False, True]).all()

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
                'dim': ['', '', '', '']}

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
        hdu = fits.new_table([c2, c1, c3], tbtype='TableHDU')

        assert (dict(hdu.data.dtype.fields) ==
                {'abc': (np.dtype('|S3'), 18),
                 'def': (np.dtype('|S15'), 2),
                 't1': (np.dtype('|S10'), 21)})
        hdu.writeto(self.temp('toto.fits'), clobber=True)
        hdul = fits.open(self.temp('toto.fits'))
        assert comparerecords(hdu.data, hdul[1].data)
        hdul.close()
        a.close()

    def test_variable_length_columns(self):
        col = fits.Column(name='QUAL_SPE', format='PJ()',
                            array=[[0] * 1571] * 225)
        tb_hdu = fits.new_table([col])
        pri_hdu = fits.PrimaryHDU()
        hdu_list = fits.HDUList([pri_hdu, tb_hdu])
        hdu_list.writeto(self.temp('toto.fits'), clobber=True)
        toto = fits.open(self.temp('toto.fits'))
        q = toto[1].data.field('QUAL_SPE')
        assert (q[0][4:8].all() ==
                np.array([0, 0, 0, 0], dtype=np.uint8)).all()
        toto.close()

    def test_extend_variable_length_array(self):
        """Regression test for issue #54."""

        arr = [[1] * 10] * 10
        col1 = fits.Column(name='TESTVLF', format='PJ()', array=arr)
        col2 = fits.Column(name='TESTSCA', format='J', array=[1] * 10)
        tb_hdu = fits.new_table([col1, col2], nrows=15)
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

    def test_endianness(self):
        x = np.ndarray((1,), dtype=object)
        channelsIn = np.array([3], dtype='uint8')
        x[0] = channelsIn
        col = fits.Column(name="Channels", format="PB()", array=x)
        cols = fits.ColDefs([col])
        tbhdu = fits.new_table(cols)
        tbhdu.name = "RFI"
        tbhdu.writeto(self.temp('testendian.fits'), clobber=True)
        hduL = fits.open(self.temp('testendian.fits'))
        rfiHDU = hduL['RFI']
        data = rfiHDU.data
        channelsOut = data.field('Channels')[0]
        assert (channelsIn == channelsOut).all()
        hduL.close()

    def test_column_endianness(self):
        """
        Regression test for #77 [PyFITS doesn't preserve byte order of
        non-native order column arrays]
        """

        a = [1., 2., 3., 4.]
        a1 = np.array(a, dtype='<f8')
        a2 = np.array(a, dtype='>f8')

        col1 = fits.Column(name='a', format='D', array=a1)
        col2 = fits.Column(name='b', format='D', array=a2)
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.new_table(cols)

        assert (tbhdu.data['a'] == a1).all()
        assert (tbhdu.data['b'] == a2).all()

        # Double check that the array is converted to the correct byte-order
        # for FITS (big-endian).
        tbhdu.writeto(self.temp('testendian.fits'), clobber=True)
        with fits.open(self.temp('testendian.fits')) as hdul:
            assert (hdul[1].data['a'] == a2).all()
            assert (hdul[1].data['b'] == a2).all()

    def test_recarray_to_bintablehdu(self):
        bright = np.rec.array([(1, 'Serius', -1.45, 'A1V'),
                               (2, 'Canopys', -0.73, 'F0Ib'),
                               (3, 'Rigil Kent', -0.1, 'G2V')],
                              formats='int16,a20,float32,a10',
                              names='order,name,mag,Sp')
        hdu = fits.BinTableHDU(bright)
        assert comparerecords(hdu.data, bright)
        hdu.writeto(self.temp('toto.fits'), clobber=True)
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
        hdu.writeto(self.temp('toto.fits'), clobber=True)
        hdul = fits.open(self.temp('toto.fits'))
        assert comparerecords(hdu.data, hdul[1].data)
        hdul.close()

    def test_new_table_from_recarray(self):
        bright = np.rec.array([(1, 'Serius', -1.45, 'A1V'),
                               (2, 'Canopys', -0.73, 'F0Ib'),
                               (3, 'Rigil Kent', -0.1, 'G2V')],
                              formats='int16,a20,float32,a10',
                              names='order,name,mag,Sp')
        hdu = fits.new_table(bright, nrows=2, tbtype='TableHDU')

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

        assert (hdu.data.field(0) == np.array([800, 2], dtype=np.int16)).all()
        assert hdu.data[0][1] == 'Serius'
        assert hdu.data[1][1] == 'Canopys'
        assert (hdu.data.field(2) ==
                np.array([-1.45, -0.73], dtype=np.float32)).all()
        assert hdu.data[0][3] == 'A1V'
        assert hdu.data[1][3] == 'F0Ib'
        with ignore_warnings():
            hdu.writeto(self.temp('toto.fits'), clobber=True)
        with fits.open(self.temp('toto.fits')) as hdul:
            assert (hdul[1].data.field(0) ==
                    np.array([800, 2], dtype=np.int16)).all()
            assert hdul[1].data[0][1] == 'Serius'
            assert hdul[1].data[1][1] == 'Canopys'
            assert (hdul[1].data.field(2) ==
                    np.array([-1.45, -0.73], dtype=np.float32)).all()
            assert hdul[1].data[0][3] == 'A1V'
            assert hdul[1].data[1][3] == 'F0Ib'
        del hdul

        hdu = fits.new_table(bright, nrows=2)
        tmp = np.rec.array([(1, 'Serius', -1.45, 'A1V'),
                            (2, 'Canopys', -0.73, 'F0Ib')],
                           formats='int16,a20,float32,a10',
                           names='order,name,mag,Sp')
        assert comparerecords(hdu.data, tmp)
        with ignore_warnings():
            hdu.writeto(self.temp('toto.fits'), clobber=True)
        with fits.open(self.temp('toto.fits')) as hdul:
            assert comparerecords(hdu.data, hdul[1].data)

    def test_new_fitsrec(self):
        """
        Tests creating a new FITS_rec object from a multi-field ndarray.
        """

        h = fits.open(self.data('tb.fits'))
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
        tbhdu = fits.new_table(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        counts = np.array([412, 434, 408, 417])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[0, 1, 0, 0])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.new_table(coldefs)
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

        assert t1[1].columns._arrays[1] is t1[1].columns.columns[1].array

        # Create a new table that consists of the data from the first table
        # but has enough space in the ndarray to hold the data from both tables
        hdu = fits.new_table(t1[1].columns, nrows=nrows)

        # For each column in the tables append the data from table 2 after the
        # data from table 1.
        for i in range(len(t1[1].columns)):
            hdu.data.field(i)[nrows1:] = t2[1].data.field(i)

        hdu.writeto(self.temp('newtable.fits'))

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 19, '8R x 5C', '[10A, J, 10A, 5E, L]',
                 '')]

        assert fits.info(self.temp('newtable.fits'), output=False) == info

        z = np.array([0.,  0.,  0.,  0.,  0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '0.0', z, True),
             ('NGC2', 334, '0.0', z, False),
             ('NGC3', 308, '0.0', z, True),
             ('NCG4', 317, '0.0', z, True),
             ('NGC5', 412, '0.0', z, False),
             ('NGC6', 434, '0.0', z, True),
             ('NGC7', 408, '0.0', z, False),
             ('NCG8', 417, '0.0', z, False)],
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
        tbhdu = fits.new_table(coldefs)

        assert tbhdu.columns.names == ['target', 'counts', 'notes', 'spectrum']
        coldefs1 = coldefs + c5

        tbhdu1 = fits.new_table(coldefs1)
        assert tbhdu1.columns.names == ['target', 'counts', 'notes',
                                        'spectrum', 'flag']

        z = np.array([0.,  0.,  0.,  0.,  0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '0.0', z, True),
             ('NGC2', 334, '0.0', z, False),
             ('NGC3', 308, '0.0', z, True),
             ('NCG4', 317, '0.0', z, True)],
             formats='a10,u4,a10,5f4,l')
        assert comparerecords(tbhdu1.data, array)

    def test_merge_tables(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.new_table(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        counts = np.array([412, 434, 408, 417])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = fits.Column(name='target1', format='10A', array=names)
        c2 = fits.Column(name='counts1', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes1', format='A10')
        c4 = fits.Column(name='spectrum1', format='5E')
        c5 = fits.Column(name='flag1', format='L', array=[0, 1, 0, 0])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = fits.new_table(coldefs)
        tbhdu.writeto(self.temp('table2.fits'))

        # Merge the columns of table 2 after the columns of table 1
        # The column names are assumed to be different

        # Open the two files we want to append
        t1 = fits.open(self.temp('table1.fits'))
        t2 = fits.open(self.temp('table2.fits'))

        hdu = fits.new_table(t1[1].columns + t2[1].columns)

        z = np.array([0.,  0.,  0.,  0.,  0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '0.0', z, True, 'NGC5', 412, '0.0', z, False),
             ('NGC2', 334, '0.0', z, False, 'NGC6', 434, '0.0', z, True),
             ('NGC3', 308, '0.0', z, True, 'NGC7', 408, '0.0', z, False),
             ('NCG4', 317, '0.0', z, True, 'NCG8', 417, '0.0', z, False)],
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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 30, '4R x 10C',
                 '[10A, J, 10A, 5E, L, 10A, J, 10A, 5E, L]', '')]

        assert fits.info(self.temp('newtable.fits'), output=False) == info

        hdul = fits.open(self.temp('newtable.fits'))
        hdu = hdul[1]

        assert hdu.columns.names == ['target', 'counts', 'notes', 'spectrum',
                                     'flag', 'target1', 'counts1', 'notes1',
                                     'spectrum1', 'flag1']

        z = np.array([0.,  0.,  0.,  0.,  0.], dtype=np.float32)
        array = np.rec.array(
            [('NGC1', 312, '0.0', z, True, 'NGC5', 412, '0.0', z, False),
             ('NGC2', 334, '0.0', z, False, 'NGC6', 434, '0.0', z, True),
             ('NGC3', 308, '0.0', z, True, 'NGC7', 408, '0.0', z, False),
             ('NCG4', 317, '0.0', z, True, 'NCG8', 417, '0.0', z, False)],
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

    def test_mask_array(self):
        t = fits.open(self.data('table.fits'))
        tbdata = t[1].data
        mask = tbdata.field('V_mag') > 12
        newtbdata = tbdata[mask]
        hdu = fits.BinTableHDU(newtbdata)
        hdu.writeto(self.temp('newtable.fits'))

        hdul = fits.open(self.temp('newtable.fits'))

        assert str(hdu.data) == "[('NGC1002', 12.3) ('NGC1003', 15.2)]"

        assert str(hdul[1].data) == "[('NGC1002', 12.3) ('NGC1003', 15.2)]"

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
        tbhdu = fits.new_table(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        t1 = fits.open(self.temp('table1.fits'))
        row = t1[1].data[2]
        assert row['counts'] == 308
        a, b, c = row[1:4]
        assert a == counts[2]
        assert b == '0.0'
        assert (c == np.array([0., 0.,  0.,  0., 0.], dtype=np.float32)).all()
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

        # Test stepping for #59
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
        tbhdu = fits.new_table(coldefs)
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

        tbhdu1 = fits.new_table(coldefs)

        c1 = fits.Column(name='target', format='10A')
        c2 = fits.Column(name='counts', format='J', unit='DN')
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L')
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.new_table(coldefs, nrows=5)

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
        assert tbhdu.columns.columns[2].array[0] == '0.0'
        assert (tbhdu.columns.columns[3].array[0] ==
                np.array([0., 0., 0., 0., 0.], dtype=np.float32)).all()
        assert tbhdu.columns.columns[4].array[0] == True

        assert tbhdu.data[3][1] == 33
        assert tbhdu.data._coldefs._arrays[1][3] == 33
        assert tbhdu.data._coldefs.columns[1].array[3] == 33
        assert tbhdu.columns._arrays[1][3] == 33
        assert tbhdu.columns.columns[1].array[3] == 33
        assert tbhdu.columns.columns[0].array[3] == 'JIM1'
        assert tbhdu.columns.columns[2].array[3] == 'A Note'
        assert (tbhdu.columns.columns[3].array[3] ==
                np.array([1., 2., 3., 4., 5.], dtype=np.float32)).all()
        assert tbhdu.columns.columns[4].array[3] == True

    def test_assign_multiple_rows_to_table(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu1 = fits.new_table(coldefs)

        counts = np.array([112, 134, 108, 117])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[0, 1, 0, 0])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.new_table(coldefs)
        tbhdu.data[0][3] = np.array([1., 2., 3., 4., 5.], dtype=np.float32)

        tbhdu2 = fits.new_table(tbhdu1.data, nrows=9)

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
        assert tbhdu2.columns.columns[2].array[0] == '0.0'
        assert (tbhdu2.columns.columns[3].array[0] ==
                np.array([0., 0., 0., 0., 0.], dtype=np.float32)).all()
        assert tbhdu2.columns.columns[4].array[0] == True

        assert tbhdu2.data[4][1] == 112
        assert tbhdu2.data._coldefs._arrays[1][4] == 112
        assert tbhdu2.data._coldefs.columns[1].array[4] == 112
        assert tbhdu2.columns._arrays[1][4] == 112
        assert tbhdu2.columns.columns[1].array[4] == 112
        assert tbhdu2.columns.columns[0].array[4] == 'NGC5'
        assert tbhdu2.columns.columns[2].array[4] == '0.0'
        assert (tbhdu2.columns.columns[3].array[4] ==
                np.array([1., 2., 3., 4., 5.], dtype=np.float32)).all()
        assert tbhdu2.columns.columns[4].array[4] == False

        assert tbhdu2.columns.columns[1].array[8] == 0
        assert tbhdu2.columns.columns[0].array[8] == '0.0'
        assert tbhdu2.columns.columns[2].array[8] == '0.0'
        assert (tbhdu2.columns.columns[3].array[8] ==
                np.array([0., 0., 0., 0., 0.], dtype=np.float32)).all()
        assert tbhdu2.columns.columns[4].array[8] == False

    def test_verify_data_references(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = fits.Column(name='target', format='10A', array=names)
        c2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = fits.Column(name='notes', format='A10')
        c4 = fits.Column(name='spectrum', format='5E')
        c5 = fits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = fits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = fits.new_table(coldefs)

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

        tbhdu = fits.new_table(coldefs)

        tbhdu1 = fits.new_table(tbhdu.data.view(np.ndarray))

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

        tbhdu = fits.new_table(coldefs)

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

        tbhdu1 = fits.new_table(fr)

        i = 0
        for row in tbhdu1.data:
            for j in xrange(len(row)):
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

        tbhdu1 = fits.new_table(coldefs)

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

    def test_bin_table_with_logical_array(self):
        c1 = fits.Column(name='flag', format='2L',
                           array=[[True, False], [False, True]])
        coldefs = fits.ColDefs([c1])

        tbhdu1 = fits.new_table(coldefs)

        assert (tbhdu1.data.field('flag')[0] ==
                np.array([True, False], dtype=np.bool)).all()
        assert (tbhdu1.data.field('flag')[1] ==
                np.array([False, True], dtype=np.bool)).all()

        tbhdu = fits.new_table(tbhdu1.data)

        assert (tbhdu.data.field('flag')[0] ==
                np.array([True, False], dtype=np.bool)).all()
        assert (tbhdu.data.field('flag')[1] ==
                np.array([False, True], dtype=np.bool)).all()

    def test_variable_length_table_format_pd_from_object_array(self):
        a = np.array([np.array([7.2e-20, 7.3e-20]), np.array([0.0]),
                      np.array([0.0])], 'O')
        acol = fits.Column(name='testa', format='PD()', array=a)
        tbhdu = fits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        tbhdu1 = fits.open(self.temp('newtable.fits'))

        for j in range(0, 3):
            for i in xrange(len(a[j])):
                assert tbhdu1[1].data.field(0)[j][i] == a[j][i]

        tbhdu1.close()

    def test_variable_length_table_format_pd_from_list(self):
        a = [np.array([7.2e-20, 7.3e-20]), np.array([0.0]), np.array([0.0])]
        acol = fits.Column(name='testa', format='PD()', array=a)
        tbhdu = fits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        tbhdu1 = fits.open(self.temp('newtable.fits'))

        for j in range(0, 3):
            for i in xrange(len(a[j])):
                assert tbhdu1[1].data.field(0)[j][i] == a[j][i]

        tbhdu1.close()

    def test_variable_length_table_format_pa_from_object_array(self):
        a = np.array([np.array(['a', 'b', 'c']), np.array(['d', 'e']),
                      np.array(['f'])], 'O')
        acol = fits.Column(name='testa', format='PA()', array=a)
        tbhdu = fits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        hdul = fits.open(self.temp('newtable.fits'))

        for j in range(0, 3):
            for i in xrange(len(a[j])):
                assert hdul[1].data.field(0)[j][i] == a[j][i]

        hdul.close()

    def test_variable_length_table_format_pa_from_list(self):
        a = ['a', 'ab', 'abc']
        acol = fits.Column(name='testa', format='PA()', array=a)
        tbhdu = fits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        hdul = fits.open(self.temp('newtable.fits'))

        for j in range(0, 3):
            for i in xrange(len(a[j])):
                assert hdul[1].data.field(0)[j][i] == a[j][i]

        hdul.close()

    def test_fits_rec_column_access(self):
        t = fits.open(self.data('table.fits'))
        tbdata = t[1].data
        assert (tbdata.V_mag == tbdata.field('V_mag')).all()
        assert (tbdata.V_mag == tbdata['V_mag']).all()

        t.close()

    def test_table_with_zero_width_column(self):
        hdul = fits.open(self.data('zerowidth.fits'))
        tbhdu = hdul[2]  # This HDU contains a zero-width column 'ORBPARM'
        assert 'ORBPARM' in tbhdu.columns.names
        # The ORBPARM column should not be in the data, though the data should
        # be readable
        assert 'ORBPARM' not in tbhdu.data.names
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
        assert 'ORBPARM' not in tbhdu.data.names
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
        ahdu = fits.new_table([acol])
        assert ahdu.data.tostring().decode('raw-unicode-escape') == s
        ahdu.writeto(self.temp('newtable.fits'))
        with fits.open(self.temp('newtable.fits')) as hdul:
            assert hdul[1].data.tostring().decode('raw-unicode-escape') == s
            assert (hdul[1].data['MEMNAME'] == a).all()
        del hdul

        ahdu = fits.new_table([acol], tbtype='TableHDU')
        with ignore_warnings():
            ahdu.writeto(self.temp('newtable.fits'), clobber=True)
        with fits.open(self.temp('newtable.fits')) as hdul:
            assert (hdul[1].data.tostring().decode('raw-unicode-escape') ==
                    s.replace('\x00', ' '))
            assert (hdul[1].data['MEMNAME'] == a).all()
            ahdu = fits.new_table(hdul[1].data.copy())
        del hdul

        # Now serialize once more as a binary table; padding bytes should
        # revert to zeroes
        ahdu.writeto(self.temp('newtable.fits'), clobber=True)
        with fits.open(self.temp('newtable.fits')) as hdul:
            assert hdul[1].data.tostring().decode('raw-unicode-escape') == s
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

        thdu = fits.new_table(data)
        # Modify the TDIM fields to my own specification
        thdu.header['TDIM1'] = '(2,3)'
        thdu.header['TDIM2'] = '(4,2)'

        thdu.writeto(self.temp('newtable.fits'))

        with fits.open(self.temp('newtable.fits')) as hdul:
            thdu = hdul[1]

            c1 = thdu.data.field(0)
            c2 = thdu.data.field(1)

            hdul.close()

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
        with ignore_warnings():
            fits.writeto(self.temp('newtable.fits'), data, clobber=True)

        t = fits.getdata(self.temp('newtable.fits'))

        assert t.field(1).dtype.str[-1] == '5'
        assert t.field(1).shape == (3, 4)

        # Like the previous test, but with an extra dimension (a bit more
        # complicated)
        data = np.zeros(3, dtype=[('x', 'f4'), ('s', 'S5', (4, 3))])
        data['x'] = 1, 2, 3
        data['s'] = 'ok'

        del t

        with ignore_warnings():
            fits.writeto(self.temp('newtable.fits'), data, clobber=True)

        t = fits.getdata(self.temp('newtable.fits'))

        assert t.field(1).dtype.str[-1] == '5'
        assert t.field(1).shape == (3, 4, 3)

    def test_string_array_round_trip(self):
        """Regression test for #201."""

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
                    recarr.field(0)[0].decode('ascii')).all()

        with fits.open(self.temp('test.fits')) as h:
            # Access the data; I think this is necessary to exhibit the bug
            # reported in #201
            h[1].data[:]
            h.writeto(self.temp('test2.fits'))

        with fits.open(self.temp('test2.fits')) as h:
            assert 'TDIM1' in h[1].header
            assert h[1].header['TDIM1'] == '(3,3,3)'
            assert len(h[1].data) == 2
            assert len(h[1].data[0]) == 1
            assert (h[1].data.field(0)[0] ==
                    recarr.field(0)[0].decode('ascii')).all()

    def test_slicing(self):
        """Regression test for #52."""

        f = fits.open(self.data('table.fits'))
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
        """Regression test for #55."""

        f = fits.open(self.data('table.fits'))
        data = f[1].data
        s1 = data[data['target'] == 'NGC1001']
        s2 = data[np.where(data['target'] == 'NGC1001')]
        s3 = data[[0]]
        s4 = data[:1]
        for s in [s1, s2, s3, s4]:
            assert isinstance(s, fits.FITS_rec)
        assert (s1 == s2).all()
        assert (s2 == s3).all()
        assert (s3 == s4).all()

    def test_array_slicing_readonly(self):
        """
        Like test_array_slicing but with the file opened in 'readonly' mode.
        Regression test for a crash when slicing readonly memmap'd tables.
        """

        f = fits.open(self.data('table.fits'), mode='readonly')
        data = f[1].data
        s1 = data[data['target'] == 'NGC1001']
        s2 = data[np.where(data['target'] == 'NGC1001')]
        s3 = data[[0]]
        s4 = data[:1]
        for s in [s1, s2, s3, s4]:
            assert isinstance(s, fits.FITS_rec)
        assert (s1 == s2).all()
        assert (s2 == s3).all()
        assert (s3 == s4).all()

    def test_dump_load_round_trip(self):
        """
        A simple test of the dump/load methods; dump the data, column, and
        header files and try to reload the table from them.
        """

        hdul = fits.open(self.data('table.fits'))
        tbhdu = hdul[1]
        datafile = self.temp('data.txt')
        cdfile = self.temp('coldefs.txt')
        hfile = self.temp('header.txt')

        tbhdu.dump(datafile, cdfile, hfile)

        new_tbhdu = fits.BinTableHDU.load(datafile, cdfile, hfile)

        assert comparerecords(tbhdu.data, new_tbhdu.data)

        # Double check that the headers are equivalent
        assert str(tbhdu.header) == str(new_tbhdu.header)

    def test_load_guess_format(self):
        """
        Tests loading a table dump with no supplied coldefs or header, so that
        the table format has to be guessed at.  There is of course no exact
        science to this; the table that's produced simply uses sensible guesses
        for that format.  Ideally this should never have to be used.
        """

        # Create a table containing a variety of data types.
        a0 = np.array([False, True, False], dtype=np.bool)
        c0 = fits.Column(name='c0', format='L', array=a0)

        # Format X currently not supported by the format
        #a1 = np.array([[0], [1], [0]], dtype=np.uint8)
        #c1 = fits.Column(name='c1', format='X', array=a1)

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

        tbhdu = fits.new_table([c0, c2, c3, c4, c5, c6, c7, c8])

        datafile = self.temp('data.txt')
        tbhdu.dump(datafile)

        new_tbhdu = fits.BinTableHDU.load(datafile)

        # In this particular case the record data at least should be equivalent
        assert comparerecords(tbhdu.data, new_tbhdu.data)

    def test_attribute_field_shadowing(self):
        """
        Regression test for #86.

        Numpy recarray objects have a poorly-considered feature of allowing
        field access by attribute lookup.  However, if a field name conincides
        with an existing attribute/method of the array, the existing name takes
        precence (making the attribute-based field lookup completely unreliable
        in general cases).

        This ensures that any FITS_rec attributes still work correctly even
        when there is a field with the same name as that attribute.
        """

        c1 = fits.Column(name='names', format='I', array=[1])
        c2 = fits.Column(name='formats', format='I', array=[2])
        c3 = fits.Column(name='other', format='I', array=[3])

        t = fits.new_table([c1, c2, c3])
        assert t.data.names == ['names', 'formats', 'other']
        assert t.data.formats == ['I'] * 3
        assert (t.data['names'] == [1]).all()
        assert (t.data['formats'] == [2]).all()
        assert (t.data.other == [3]).all()

    def test_table_from_bool_fields(self):
        """
        Regression test for #113.

        Tests creating a table from a recarray containing numpy.bool columns.
        """

        array = np.rec.array([(True, False), (False, True)], formats='|b1,|b1')
        thdu = fits.new_table(array)
        assert thdu.columns.formats == ['L', 'L']
        assert comparerecords(thdu.data, array)

        # Test round trip
        thdu.writeto(self.temp('table.fits'))
        data = fits.getdata(self.temp('table.fits'), ext=1)
        assert thdu.columns.formats == ['L', 'L']
        assert comparerecords(data, array)

    def test_bool_column_update(self):
        """Regression test for #139."""

        c1 = fits.Column('F1', 'L', array=[True, False])
        c2 = fits.Column('F2', 'L', array=[False, True])
        thdu = fits.new_table(fits.ColDefs([c1, c2]))
        thdu.writeto(self.temp('table.fits'))

        with fits.open(self.temp('table.fits'), mode='update') as hdul:
            hdul[1].data['F1'][1] = True
            hdul[1].data['F2'][0] = True

        with fits.open(self.temp('table.fits')) as hdul:
            assert (hdul[1].data['F1'] == [True, True]).all()
            assert (hdul[1].data['F2'] == [True, True]).all()

    def test_missing_tnull(self):
        """Regression test for #197."""

        c = fits.Column('F1', 'A3', null='---',
                          array=np.array(['1.0', '2.0', '---', '3.0']))
        table = fits.new_table([c], tbtype='TableHDU')
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
        except ValueError, e:
            assert str(e).endswith(
                         "the header may be missing the necessary TNULL1 "
                         "keyword or the table contains invalid data")
