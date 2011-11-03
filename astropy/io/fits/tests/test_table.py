from __future__ import division # confidence high
from __future__ import with_statement

import numpy as np
from numpy import char as chararray

import pyfits
from pyfits.util import decode_ascii
from pyfits.tests import PyfitsTestCase

from nose.tools import (assert_equal, assert_not_equal, assert_raises,
                        assert_true)


def comparefloats(a, b):
    """
    Compare two float scalars or arrays and see if they are consistent

    Consistency is determined ensuring the difference is less than the
    expected amount. Return True if consistent, False if any differences.
    """

    aa = a
    bb = b
    # compute expected precision
    if aa.dtype.name=="float32" or bb.dtype.name=='float32':
        precision = 0.000001
    else:
        precision = 0.0000000000000001
    precision = 0.00001 # until precision problem is fixed in pyfits
    diff = np.absolute(aa - bb)
    mask0 = aa == 0
    masknz = aa != 0.
    if np.any(mask0):
        if diff[mask0].max() != 0.:
            return False
    if np.any(masknz):
        if (diff[masknz]/np.absolute(aa[masknz])).max() > precision:
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
            print "type(fielda): ",type(fielda)," fielda: ",fielda
            print "type(fieldb): ",type(fieldb)," fieldb: ",fieldb
            print 'field %d type differs' % i
            return False
        if isinstance(fielda[0], np.floating):
            if not comparefloats(fielda, fieldb):
                print "fielda: ",fielda
                print "fieldb: ",fieldb
                print 'field %d differs' % i
                return False
        elif (isinstance(fielda, pyfits.column._VLF) or
              isinstance(fieldb, pyfits.column._VLF)):
            for row in range(len(fielda)):
                if np.any(fielda[row] != fieldb[row]):
                    print 'fielda[%d]: %s' % (row, fielda[row])
                    print 'fieldb[%d]: %s' % (row, fieldb[row])
                    print 'field %d differs in row %d' (i, row)
        else:
            if np.any(fielda != fieldb):
                print "fielda: ",fielda
                print "fieldb: ",fieldb
                print 'field %d differs' % i
                return False
    return True


class TestTableFunctions(PyfitsTestCase):
    def test_open(self):
        # open some existing FITS files:
        tt = pyfits.open(self.data('tb.fits'))
        fd = pyfits.open(self.data('test0.fits'))

        # create some local arrays
        a1 = chararray.array(['abc', 'def', 'xx'])
        r1 = np.array([11.,12.,13.], dtype=np.float32)

        # create a table from scratch, using a mixture of columns from existing
        # tables and locally created arrays:

        # first, create individual column definitions

        c1 = pyfits.Column(name='abc', format='3A', array=a1)
        c2 = pyfits.Column(name='def', format='E', array=r1)
        a3 = np.array([3,4,5], dtype='i2')
        c3 = pyfits.Column(name='xyz', format='I', array=a3)
        a4 = np.array([1,2,3], dtype='i2')
        c4 = pyfits.Column(name='t1', format='I', array=a4)
        a5 = np.array([3+3j,4+4j,5+5j], dtype='c8')
        c5 = pyfits.Column(name='t2', format='C', array=a5)

        # Note that X format must be two-D array
        a6 = np.array([[0], [1], [0]], dtype=np.uint8)
        c6 = pyfits.Column(name='t3', format='X', array=a6)
        a7 = np.array([101, 102, 103],dtype='i4')
        c7 = pyfits.Column(name='t4', format='J', array=a7)
        a8 = np.array([[1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1],
                       [0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0],
                       [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1]], dtype=np.uint8)
        c8=pyfits.Column(name='t5', format='11X', array=a8)

        # second, create a column-definitions object for all columns in a table

        x = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8])

        # create a new binary table HDU object by using the new_table function

        tbhdu = pyfits.new_table(x)

        # another way to create a table is by using existing table's information:

        x2 = pyfits.ColDefs(tt[1])
        t2 = pyfits.new_table(x2, nrows=2)
        ra = np.rec.array([
            (1, 'abc', 3.7000002861022949, 0),
            (2, 'xy ', 6.6999998092651367, 1)], names='c1, c2, c3, c4')

        assert_equal(comparerecords(t2.data, ra), True)

        # the table HDU's data is a subclass of a record array, so we can access
        # one row like this:

        assert_equal(tbhdu.data[1][0], a1[1])
        assert_equal(tbhdu.data[1][1], r1[1])
        assert_equal(tbhdu.data[1][2], a3[1])
        assert_equal(tbhdu.data[1][3], a4[1])
        assert_equal(tbhdu.data[1][4], a5[1])
        assert_true((tbhdu.data[1][5] == a6[1].view('bool')).all())
        assert_equal(tbhdu.data[1][6], a7[1])
        assert_equal(tbhdu.data[1][7].all(), a8[1].all())

        # and a column like this:
        assert_equal(str(tbhdu.data.field('abc')), "['abc' 'def' 'xx']")

        # An alternative way to create a column-definitions object is from an
        # existing table.
        xx = pyfits.ColDefs(tt[1])

        # now we write out the newly created table HDU to a FITS file:
        fout = pyfits.HDUList(pyfits.PrimaryHDU())
        fout.append(tbhdu)
        fout.writeto(self.temp('tableout1.fits'), clobber=True)

        f2 = pyfits.open(self.temp('tableout1.fits'))
        temp = f2[1].data.field(7)
        assert_true((temp[0] == [True, True, False, True, False, True, True,
                                 True, False, False, True]).all())
        f2.close()

        # An alternative way to create an output table FITS file:
        fout2 = pyfits.open(self.temp('tableout2.fits'), 'append')
        fout2.append(fd[0])
        fout2.append(tbhdu)
        fout2.close()
        tt.close()
        fd.close()

    def test_binary_table(self):
        # binary table:
        t = pyfits.open(self.data('tb.fits'))
        assert_equal(t[1].header['tform1'], '1J')

        info = {'name': ['c1', 'c2', 'c3', 'c4'],
                'format': ['1J', '3A', '1E', '1L'],
                'unit': ['', '', '', ''],
                'null': [-2147483647, '', '', ''],
                'bscale': ['', '', 3, ''],
                'bzero': ['', '', 0.4, ''],
                'disp': ['I11', 'A3', 'G15.7', 'L6'],
                'start': ['', '', '', ''],
                'dim': ['', '', '', '']}

        assert_equal(t[1].columns.info(output=False), info)

        ra = np.rec.array([
            (1, 'abc', 3.7000002861022949, 0),
            (2, 'xy ', 6.6999998092651367, 1)], names='c1, c2, c3, c4')

        assert_equal(comparerecords(t[1].data, ra[:2]), True)

        # Change scaled field and scale back to the original array
        t[1].data.field('c4')[0] = 1
        t[1].data._scale_back()
        assert_equal(str(np.rec.recarray.field(t[1].data,'c4')), "[84 84]")

        # look at data column-wise
        assert_equal(t[1].data.field(0).all(), np.array([1, 2]).all())

        # When there are scaled columns, the raw data are in data._parent

        t.close()

    def test_ascii_table(self):
        # ASCII table
        a = pyfits.open(self.data('ascii.fits'))
        ra1 = np.rec.array([
            (10.123000144958496, 37),
            (5.1999998092651367, 23),
            (15.609999656677246, 17),
            (0.0, 0),
            (345.0, 345)], names='c1, c2')
        assert_equal(comparerecords(a[1].data, ra1), True)

        # Test slicing
        a2 = a[1].data[2:][2:]
        ra2 = np.rec.array([(345.0,345)],names='c1, c2')

        assert_equal(comparerecords(a2, ra2), True)

        assert_equal(a2.field(1).all(),np.array([345]).all())

        ra3 = np.rec.array([
            (10.123000144958496, 37),
            (15.609999656677246, 17),
            (345.0, 345)
            ], names='c1, c2')

        assert_equal(comparerecords(a[1].data[::2], ra3), True)

        # Test Start Column

        a1 = chararray.array(['abcd','def'])
        r1 = np.array([11.,12.])
        c1 = pyfits.Column(name='abc',format='A3',start=19, array=a1)
        c2 = pyfits.Column(name='def',format='E',start=3, array=r1)
        c3 = pyfits.Column(name='t1',format='I',array=[91, 92, 93])
        hdu = pyfits.new_table([c2, c1, c3],tbtype='TableHDU')


        assert_equal(dict(hdu.data.dtype.fields),
                         {'abc': (np.dtype('|S3'), 18),
                          'def': (np.dtype('|S15'), 2),
                          't1': (np.dtype('|S10'), 21)})
        hdu.writeto(self.temp('toto.fits'), clobber=True)
        hdul = pyfits.open(self.temp('toto.fits'))
        assert_equal(comparerecords(hdu.data,hdul[1].data), True)
        hdul.close()
        a.close()

    def test_variable_length_columns(self):
        col = pyfits.Column(name='QUAL_SPE', format='PJ()',
                            array=[[0]*1571]*225)
        tb_hdu = pyfits.new_table([col])
        pri_hdu = pyfits.PrimaryHDU()
        hdu_list = pyfits.HDUList([pri_hdu,tb_hdu])
        hdu_list.writeto(self.temp('toto.fits'), clobber=True)
        toto = pyfits.open(self.temp('toto.fits'))
        q = toto[1].data.field('QUAL_SPE')
        assert_equal(q[0][4:8].all(),
                         np.array([0, 0, 0, 0],dtype=np.uint8).all())
        toto.close()

    def test_extend_variable_length_array(self):
        """Regression test for issue #54."""

        arr = [[1] * 10] * 10
        col1 = pyfits.Column(name='TESTVLF', format='PJ()', array=arr)
        col2 = pyfits.Column(name='TESTSCA', format='J', array=[1] * 10)
        tb_hdu = pyfits.new_table([col1, col2], nrows=15)
        # This asserts that the normal 'scalar' column's length was extended
        assert_equal(len(tb_hdu.data['TESTSCA']), 15)
        # And this asserts that the VLF column was extended in the same manner
        assert_equal(len(tb_hdu.data['TESTVLF']), 15)
        # We can't compare the whole array since the _VLF is an array of
        # objects, but comparing just the edge case rows should suffice
        assert_true((tb_hdu.data['TESTVLF'][0] == arr[0]).all())
        assert_true((tb_hdu.data['TESTVLF'][9] == arr[9]).all())
        assert_true((tb_hdu.data['TESTVLF'][10] == ([0] * 10)).all())
        assert_true((tb_hdu.data['TESTVLF'][-1] == ([0] * 10)).all())

    def test_endianness(self):
        x = np.ndarray((1,), dtype=object)
        channelsIn = np.array([3], dtype='uint8')
        x[0] = channelsIn
        col = pyfits.Column(name="Channels", format="PB()", array=x)
        cols = pyfits.ColDefs([col])
        tbhdu = pyfits.new_table(cols)
        tbhdu.name = "RFI"
        tbhdu.writeto(self.temp('testendian.fits'), clobber=True)
        hduL = pyfits.open(self.temp('testendian.fits'))
        rfiHDU = hduL['RFI']
        data = rfiHDU.data
        channelsOut = data.field('Channels')[0]
        assert_equal(channelsIn.all(),channelsOut.all())
        hduL.close()

    def test_column_endianness(self):
        """
        Regression test for #77 [PyFITS doesn't preserve byte order of
        non-native order column arrays]
        """

        a = [1., 2., 3., 4.]
        a1 = np.array(a, dtype='<f8')
        a2 = np.array(a, dtype='>f8')

        col1 = pyfits.Column(name='a', format='D', array=a1)
        col2 = pyfits.Column(name='b', format='D', array=a2)
        cols = pyfits.ColDefs([col1, col2])
        tbhdu = pyfits.new_table(cols)

        assert_true((tbhdu.data['a'] == a1).all())
        assert_true((tbhdu.data['b'] == a2).all())

        # Double check that the array is converted to the correct byte-order
        # for FITS (big-endian).
        tbhdu.writeto(self.temp('testendian.fits'), clobber=True)
        hdul = pyfits.open(self.temp('testendian.fits'))
        assert_true((hdul[1].data['a'] == a2).all())
        assert_true((hdul[1].data['b'] == a2).all())

    def test_recarray_to_bintablehdu(self):
        bright=np.rec.array([(1,'Serius',-1.45,'A1V'),\
                             (2,'Canopys',-0.73,'F0Ib'),\
                             (3,'Rigil Kent',-0.1,'G2V')],\
                            formats='int16,a20,float32,a10',\
                            names='order,name,mag,Sp')
        hdu=pyfits.BinTableHDU(bright)
        assert_equal(comparerecords(hdu.data, bright), True)
        hdu.writeto(self.temp('toto.fits'), clobber=True)
        hdul = pyfits.open(self.temp('toto.fits'))
        assert_equal(comparerecords(hdu.data,hdul[1].data),True)
        assert_equal(comparerecords(bright,hdul[1].data),True)
        hdul.close()

    def test_numpy_ndarray_to_bintablehdu(self):
        desc = np.dtype({'names': ['order','name','mag','Sp'],
                         'formats': ['int','S20','float32','S10']})
        a = np.array([(1,'Serius',-1.45,'A1V'),
                      (2,'Canopys',-0.73,'F0Ib'),
                      (3,'Rigil Kent',-0.1,'G2V')], dtype=desc)
        hdu=pyfits.BinTableHDU(a)
        assert_equal(comparerecords(hdu.data, a.view(pyfits.FITS_rec)),
                         True)
        hdu.writeto(self.temp('toto.fits'), clobber=True)
        hdul = pyfits.open(self.temp('toto.fits'))
        assert_equal(comparerecords(hdu.data,hdul[1].data),True)
        hdul.close()

    def test_new_table_from_recarray(self):
        bright = np.rec.array([(1,'Serius',-1.45,'A1V'),
                            (2,'Canopys',-0.73,'F0Ib'),
                            (3,'Rigil Kent',-0.1,'G2V')],
                           formats='int16,a20,float32,a10',
                           names='order,name,mag,Sp')
        hdu=pyfits.new_table(bright,nrows=2,tbtype='TableHDU')

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert_equal(id(hdu.data._coldefs.columns[0].array),
                         id(hdu.data._coldefs._arrays[0]))
        assert_equal(id(hdu.data._coldefs.columns[0].array),
                         id(hdu.columns.data[0].array))
        assert_equal(id(hdu.data._coldefs.columns[0].array),
                         id(hdu.columns._arrays[0]))

        # Ensure I can change the value of one data element and it effects
        # all of the others.
        hdu.data[0][0] = 213

        assert_equal(hdu.data[0][0], 213)
        assert_equal(hdu.data._coldefs._arrays[0][0], 213)
        assert_equal(hdu.data._coldefs.columns[0].array[0], 213)
        assert_equal(hdu.columns._arrays[0][0], 213)
        assert_equal(hdu.columns.data[0].array[0], 213)

        hdu.data._coldefs._arrays[0][0] = 100

        assert_equal(hdu.data[0][0], 100)
        assert_equal(hdu.data._coldefs._arrays[0][0], 100)
        assert_equal(hdu.data._coldefs.columns[0].array[0], 100)
        assert_equal(hdu.columns._arrays[0][0], 100)
        assert_equal(hdu.columns.data[0].array[0], 100)

        hdu.data._coldefs.columns[0].array[0] = 500
        assert_equal(hdu.data[0][0], 500)
        assert_equal(hdu.data._coldefs._arrays[0][0], 500)
        assert_equal(hdu.data._coldefs.columns[0].array[0], 500)
        assert_equal(hdu.columns._arrays[0][0], 500)
        assert_equal(hdu.columns.data[0].array[0], 500)

        hdu.columns._arrays[0][0] = 600
        assert_equal(hdu.data[0][0], 600)
        assert_equal(hdu.data._coldefs._arrays[0][0], 600)
        assert_equal(hdu.data._coldefs.columns[0].array[0], 600)
        assert_equal(hdu.columns._arrays[0][0], 600)
        assert_equal(hdu.columns.data[0].array[0], 600)

        hdu.columns.data[0].array[0] = 800
        assert_equal(hdu.data[0][0], 800)
        assert_equal(hdu.data._coldefs._arrays[0][0], 800)
        assert_equal(hdu.data._coldefs.columns[0].array[0], 800)
        assert_equal(hdu.columns._arrays[0][0], 800)
        assert_equal(hdu.columns.data[0].array[0], 800)

        assert_equal(hdu.data.field(0).all(),
                         np.array([1, 2],dtype=np.int16).all())
        assert_equal(hdu.data[0][1], 'Serius')
        assert_equal(hdu.data[1][1], 'Canopys')
        assert_equal(hdu.data.field(2).all(),
                         np.array([-1.45, -0.73], dtype=np.float32).all())
        assert_equal(hdu.data[0][3], 'A1V')
        assert_equal(hdu.data[1][3], 'F0Ib')
        hdu.writeto(self.temp('toto.fits'), clobber=True)
        hdul = pyfits.open(self.temp('toto.fits'))
        assert_equal(hdul[1].data.field(0).all(),
                         np.array([1, 2], dtype=np.int16).all())
        assert_equal(hdul[1].data[0][1], 'Serius')
        assert_equal(hdul[1].data[1][1], 'Canopys')
        assert_equal(hdul[1].data.field(2).all(),
                         np.array([-1.45, -0.73], dtype=np.float32).all())
        assert_equal(hdul[1].data[0][3], 'A1V')
        assert_equal(hdul[1].data[1][3], 'F0Ib')

        hdul.close()

        hdu=pyfits.new_table(bright,nrows=2)
        tmp=np.rec.array([(1,'Serius',-1.45,'A1V'),
                          (2,'Canopys',-0.73,'F0Ib')],
                         formats='int16,a20,float32,a10',
                         names='order,name,mag,Sp')
        assert_equal(comparerecords(hdu.data,tmp), True)
        hdu.writeto(self.temp('toto.fits'), clobber=True)
        hdul = pyfits.open(self.temp('toto.fits'))
        assert_equal(comparerecords(hdu.data,hdul[1].data),True)
        hdul.close()

    def test_new_fitsrec(self):
        """
        Tests creating a new FITS_rec object from a multi-field ndarray.
        """

        h = pyfits.open(self.data('tb.fits'))
        data = h[1].data
        new_data = np.array([(3, 'qwe', 4.5, False)], dtype=data.dtype)
        appended = np.append(data, new_data).view(pyfits.FITS_rec)
        assert_true(repr(appended).startswith('FITS_rec('))
        # This test used to check the entire string representation of FITS_rec,
        # but that has problems between different numpy versions.  Instead just
        # check that the FITS_rec was created, and we'll let subsequent tests
        # worry about checking values and such

    def test_appending_a_column(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu=pyfits.new_table(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        counts = np.array([412, 434, 408, 417])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[0, 1, 0, 0])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.writeto(self.temp('table2.fits'))

        # Append the rows of table 2 after the rows of table 1
        # The column definitions are assumed to be the same

        # Open the two files we want to append
        t1 = pyfits.open(self.temp('table1.fits'))
        t2 = pyfits.open(self.temp('table2.fits'))

        # Get the number of rows in the table from the first file
        nrows1 = t1[1].data.shape[0]

        # Get the total number of rows in the resulting appended table
        nrows = t1[1].data.shape[0] + t2[1].data.shape[0]

        assert_equal(t1[1].columns._arrays[1] is
                         t1[1].columns.columns[1].array, True)

        # Create a new table that consists of the data from the first table
        # but has enough space in the ndarray to hold the data from both tables
        hdu = pyfits.new_table(t1[1].columns, nrows=nrows)

        # For each column in the tables append the data from table 2 after the
        # data from table 1.
        for i in range(len(t1[1].columns)):
            hdu.data.field(i)[nrows1:] = t2[1].data.field(i)

        hdu.writeto(self.temp('newtable.fits'))

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 19, '8R x 5C', '[10A, J, 10A, 5E, L]',
                 '')]

        assert_equal(pyfits.info(self.temp('newtable.fits'), output=False), info)

        array = np.rec.array(
            [('NGC1', 312, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NGC2', 334, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NGC3', 308, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NCG4', 317, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NGC5', 412, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NGC6', 434, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NGC7', 408, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NCG8', 417, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False)],
             formats='a10,u4,a10,5f4,l')

        assert_true(comparerecords(hdu.data, array))

        # Verify that all of the references to the data point to the same
        # numarray
        hdu.data[0][1] = 300
        assert_equal(hdu.data._coldefs._arrays[1][0], 300)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 300)
        assert_equal(hdu.columns._arrays[1][0], 300)
        assert_equal(hdu.columns.data[1].array[0], 300)
        assert_equal(hdu.data[0][1], 300)

        hdu.data._coldefs._arrays[1][0] = 200
        assert_equal(hdu.data._coldefs._arrays[1][0], 200)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 200)
        assert_equal(hdu.columns._arrays[1][0], 200)
        assert_equal(hdu.columns.data[1].array[0], 200)
        assert_equal(hdu.data[0][1], 200)

        hdu.data._coldefs.columns[1].array[0] = 100
        assert_equal(hdu.data._coldefs._arrays[1][0], 100)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 100)
        assert_equal(hdu.columns._arrays[1][0], 100)
        assert_equal(hdu.columns.data[1].array[0], 100)
        assert_equal(hdu.data[0][1], 100)

        hdu.columns._arrays[1][0] = 90
        assert_equal(hdu.data._coldefs._arrays[1][0], 90)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 90)
        assert_equal(hdu.columns._arrays[1][0], 90)
        assert_equal(hdu.columns.data[1].array[0], 90)
        assert_equal(hdu.data[0][1], 90)

        hdu.columns.data[1].array[0] = 80
        assert_equal(hdu.data._coldefs._arrays[1][0], 80)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 80)
        assert_equal(hdu.columns._arrays[1][0], 80)
        assert_equal(hdu.columns.data[1].array[0], 80)
        assert_equal(hdu.data[0][1], 80)

        # Same verification from the file
        hdul = pyfits.open(self.temp('newtable.fits'))
        hdu = hdul[1]
        hdu.data[0][1] = 300
        assert_equal(hdu.data._coldefs._arrays[1][0], 300)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 300)
        assert_equal(hdu.columns._arrays[1][0], 300)
        assert_equal(hdu.columns.data[1].array[0], 300)
        assert_equal(hdu.data[0][1], 300)

        hdu.data._coldefs._arrays[1][0] = 200
        assert_equal(hdu.data._coldefs._arrays[1][0], 200)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 200)
        assert_equal(hdu.columns._arrays[1][0], 200)
        assert_equal(hdu.columns.data[1].array[0], 200)
        assert_equal(hdu.data[0][1], 200)

        hdu.data._coldefs.columns[1].array[0] = 100
        assert_equal(hdu.data._coldefs._arrays[1][0], 100)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 100)
        assert_equal(hdu.columns._arrays[1][0], 100)
        assert_equal(hdu.columns.data[1].array[0], 100)
        assert_equal(hdu.data[0][1], 100)

        hdu.columns._arrays[1][0] = 90
        assert_equal(hdu.data._coldefs._arrays[1][0], 90)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 90)
        assert_equal(hdu.columns._arrays[1][0], 90)
        assert_equal(hdu.columns.data[1].array[0], 90)
        assert_equal(hdu.data[0][1], 90)

        hdu.columns.data[1].array[0] = 80
        assert_equal(hdu.data._coldefs._arrays[1][0], 80)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 80)
        assert_equal(hdu.columns._arrays[1][0], 80)
        assert_equal(hdu.columns.data[1].array[0], 80)
        assert_equal(hdu.data[0][1], 80)

        t1.close()
        t2.close()
        hdul.close()

    def test_adding_a_column(self):
        # Tests adding a column to a table.
        counts = np.array([312,334,308,317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4])
        tbhdu = pyfits.new_table(coldefs)

        assert_equal(tbhdu.columns.names,
                         ['target', 'counts', 'notes', 'spectrum'])
        coldefs1 = coldefs + c5

        tbhdu1=pyfits.new_table(coldefs1)
        assert_equal(tbhdu1.columns.names,
                         ['target', 'counts', 'notes', 'spectrum', 'flag'])

        array = np.rec.array(
            [('NGC1', 312, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NGC2', 334, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NGC3', 308, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NCG4', 317, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True)],
             formats='a10,u4,a10,5f4,l')
        assert_true(comparerecords(tbhdu1.data, array))

    def test_merge_tables(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        counts = np.array([412, 434, 408, 417])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = pyfits.Column(name='target1', format='10A', array=names)
        c2 = pyfits.Column(name='counts1', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes1', format='A10')
        c4 = pyfits.Column(name='spectrum1',format='5E')
        c5 = pyfits.Column(name='flag1',format='L',array=[0,1,0,0])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.writeto(self.temp('table2.fits'))

        # Merge the columns of table 2 after the columns of table 1
        # The column names are assumed to be different

        # Open the two files we want to append
        t1 = pyfits.open(self.temp('table1.fits'))
        t2 = pyfits.open(self.temp('table2.fits'))

        hdu = pyfits.new_table(t1[1].columns+t2[1].columns)

        array = np.rec.array(
            [('NGC1', 312, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True, 'NGC5', 412, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NGC2', 334, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False, 'NGC6', 434, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NGC3', 308, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True, 'NGC7', 408, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NCG4', 317, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True, 'NCG8', 417, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False)],
             formats='a10,u4,a10,5f4,l,a10,u4,a10,5f4,l')
        assert_true(comparerecords(hdu.data, array))

        hdu.writeto(self.temp('newtable.fits'))

        # Verify that all of the references to the data point to the same
        # numarray
        hdu.data[0][1] = 300
        assert_equal(hdu.data._coldefs._arrays[1][0], 300)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 300)
        assert_equal(hdu.columns._arrays[1][0], 300)
        assert_equal(hdu.columns.data[1].array[0], 300)
        assert_equal(hdu.data[0][1], 300)

        hdu.data._coldefs._arrays[1][0] = 200
        assert_equal(hdu.data._coldefs._arrays[1][0], 200)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 200)
        assert_equal(hdu.columns._arrays[1][0], 200)
        assert_equal(hdu.columns.data[1].array[0], 200)
        assert_equal(hdu.data[0][1], 200)

        hdu.data._coldefs.columns[1].array[0] = 100
        assert_equal(hdu.data._coldefs._arrays[1][0], 100)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 100)
        assert_equal(hdu.columns._arrays[1][0], 100)
        assert_equal(hdu.columns.data[1].array[0], 100)
        assert_equal(hdu.data[0][1], 100)

        hdu.columns._arrays[1][0] = 90
        assert_equal(hdu.data._coldefs._arrays[1][0], 90)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 90)
        assert_equal(hdu.columns._arrays[1][0], 90)
        assert_equal(hdu.columns.data[1].array[0], 90)
        assert_equal(hdu.data[0][1], 90)

        hdu.columns.data[1].array[0] = 80
        assert_equal(hdu.data._coldefs._arrays[1][0], 80)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 80)
        assert_equal(hdu.columns._arrays[1][0], 80)
        assert_equal(hdu.columns.data[1].array[0], 80)
        assert_equal(hdu.data[0][1], 80)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 30, '4R x 10C',
                 '[10A, J, 10A, 5E, L, 10A, J, 10A, 5E, L]', '')]

        assert_equal(pyfits.info(self.temp('newtable.fits'), output=False), info)

        hdul = pyfits.open(self.temp('newtable.fits'))
        hdu = hdul[1]

        assert_equal(hdu.columns.names,
                         ['target', 'counts', 'notes', 'spectrum', 'flag',
                          'target1', 'counts1', 'notes1', 'spectrum1', 'flag1'])

        array = np.rec.array(
            [('NGC1', 312, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True, 'NGC5', 412, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NGC2', 334, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False, 'NGC6', 434, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True),
             ('NGC3', 308, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True, 'NGC7', 408, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False),
             ('NCG4', 317, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), True, 'NCG8', 417, '0.0', np.array([ 0.,  0.,  0.,  0.,  0.], dtype=np.float32), False)],
             formats='a10,u4,a10,5f4,l,a10,u4,a10,5f4,l')
        assert_true(comparerecords(hdu.data, array))

        # Same verification from the file
        hdu.data[0][1] = 300
        assert_equal(hdu.data._coldefs._arrays[1][0], 300)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 300)
        assert_equal(hdu.columns._arrays[1][0], 300)
        assert_equal(hdu.columns.data[1].array[0], 300)
        assert_equal(hdu.data[0][1], 300)

        hdu.data._coldefs._arrays[1][0] = 200
        assert_equal(hdu.data._coldefs._arrays[1][0], 200)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 200)
        assert_equal(hdu.columns._arrays[1][0], 200)
        assert_equal(hdu.columns.data[1].array[0], 200)
        assert_equal(hdu.data[0][1], 200)

        hdu.data._coldefs.columns[1].array[0] = 100
        assert_equal(hdu.data._coldefs._arrays[1][0], 100)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 100)
        assert_equal(hdu.columns._arrays[1][0], 100)
        assert_equal(hdu.columns.data[1].array[0], 100)
        assert_equal(hdu.data[0][1], 100)

        hdu.columns._arrays[1][0] = 90
        assert_equal(hdu.data._coldefs._arrays[1][0], 90)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 90)
        assert_equal(hdu.columns._arrays[1][0], 90)
        assert_equal(hdu.columns.data[1].array[0], 90)
        assert_equal(hdu.data[0][1], 90)

        hdu.columns.data[1].array[0] = 80
        assert_equal(hdu.data._coldefs._arrays[1][0], 80)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 80)
        assert_equal(hdu.columns._arrays[1][0], 80)
        assert_equal(hdu.columns.data[1].array[0], 80)
        assert_equal(hdu.data[0][1], 80)

        t1.close()
        t2.close()
        hdul.close()

    def test_mask_array(self):
        t = pyfits.open(self.data('table.fits'))
        tbdata = t[1].data
        mask = tbdata.field('V_mag') > 12
        newtbdata = tbdata[mask]
        hdu = pyfits.BinTableHDU(newtbdata)
        hdu.writeto(self.temp('newtable.fits'))

        hdul = pyfits.open(self.temp('newtable.fits'))

        assert_equal(str(hdu.data),
                         "[('NGC1002', 12.3) ('NGC1003', 15.2)]")

        assert_equal(str(hdul[1].data),
                         "[('NGC1002', 12.3) ('NGC1003', 15.2)]")

        t.close()
        hdul.close()

    def test_slice_a_row(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        t1=pyfits.open(self.temp('table1.fits'))
        row = t1[1].data[2]
        assert_equal(row['counts'], 308)
        a,b,c = row[1:4]
        assert_equal(a, counts[2])
        assert_equal(b, '0.0')
        assert_equal(c.all(), np.array([ 0.,  0.,  0.,  0.,  0.],
                                           dtype=np.float32).all())
        row['counts'] = 310
        assert_equal(row['counts'], 310)

        row[1] = 315
        assert_equal(row['counts'], 315)

        assert_equal(row[1:4]['counts'], 315)

        assert_raises(KeyError, lambda r: r[1:4]['flag'], row)

        row[1:4]['counts'] = 300
        assert_equal(row[1:4]['counts'], 300)
        assert_equal(row['counts'], 300)

        row[1:4][0] = 400
        assert_equal(row[1:4]['counts'], 400)
        row[1:4]['counts'] = 300
        assert_equal(row[1:4]['counts'], 300)

        # Test stepping for #59
        row[1:4][::-1][-1] = 500
        assert_equal(row[1:4]['counts'], 500)
        row[1:4:2][0] = 300
        assert_equal(row[1:4]['counts'], 300)

        assert_raises(KeyError, lambda r: r[1:4]['flag'], row)

        assert_equal(row[1:4].field(0), 300)
        assert_equal(row[1:4].field('counts'), 300)

        assert_raises(KeyError, row[1:4].field, 'flag')

        row[1:4].setfield('counts', 500)
        assert_equal(row[1:4].field(0), 500)

        assert_raises(KeyError, row[1:4].setfield, 'flag', False)

        assert_equal(t1[1].data._coldefs._arrays[1][2], 500)
        assert_equal(t1[1].data._coldefs.columns[1].array[2], 500)
        assert_equal(t1[1].columns._arrays[1][2], 500)
        assert_equal(t1[1].columns.data[1].array[2], 500)
        assert_equal(t1[1].data[2][1], 500)

        t1.close()

    def test_fits_record_len(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
        tbhdu = pyfits.new_table(coldefs)
        tbhdu.writeto(self.temp('table1.fits'))

        t1 = pyfits.open(self.temp('table1.fits'))

        assert_equal(len(t1[1].data[0]), 5)
        assert_equal(len(t1[1].data[0][0:4]), 4)
        assert_equal(len(t1[1].data[0][0:5]), 5)
        assert_equal(len(t1[1].data[0][0:6]), 5)
        assert_equal(len(t1[1].data[0][0:7]), 5)
        assert_equal(len(t1[1].data[0][1:4]), 3)
        assert_equal(len(t1[1].data[0][1:5]), 4)
        assert_equal(len(t1[1].data[0][1:6]), 4)
        assert_equal(len(t1[1].data[0][1:7]), 4)

        t1.close()

    def test_add_data_by_rows(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu1=pyfits.new_table(coldefs)

        c1 = pyfits.Column(name='target', format='10A')
        c2 = pyfits.Column(name='counts', format='J', unit='DN')
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L')
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = pyfits.new_table(coldefs, nrows=5)

        # Test assigning data to a tables row using a FITS_record
        tbhdu.data[0] = tbhdu1.data[0]
        tbhdu.data[4] = tbhdu1.data[3]

        # Test assigning data to a tables row using a tuple
        tbhdu.data[2] = ('NGC1', 312, 'A Note',
                         np.array([1.1, 2.2, 3.3, 4.4, 5.5], dtype=np.float32),
                         True)

        # Test assigning data to a tables row using a list
        tbhdu.data[3] = ['JIM1', '33', 'A Note',
                         np.array([1., 2., 3., 4., 5.],dtype=np.float32),True]

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert_equal(id(tbhdu.data._coldefs.columns[0].array),
                         id(tbhdu.data._coldefs._arrays[0]))
        assert_equal(id(tbhdu.data._coldefs.columns[0].array),
                         id(tbhdu.columns.data[0].array))
        assert_equal(id(tbhdu.data._coldefs.columns[0].array),
                         id(tbhdu.columns._arrays[0]))

        assert_equal(tbhdu.data[0][1], 312)
        assert_equal(tbhdu.data._coldefs._arrays[1][0], 312)
        assert_equal(tbhdu.data._coldefs.columns[1].array[0], 312)
        assert_equal(tbhdu.columns._arrays[1][0], 312)
        assert_equal(tbhdu.columns.data[1].array[0], 312)
        assert_equal(tbhdu.columns.data[0].array[0], 'NGC1')
        assert_equal(tbhdu.columns.data[2].array[0], '0.0')
        assert_equal(tbhdu.columns.data[3].array[0].all(),
                         np.array([0., 0., 0., 0., 0.],dtype=np.float32).all())
        assert_equal(tbhdu.columns.data[4].array[0], True)

        assert_equal(tbhdu.data[3][1], 33)
        assert_equal(tbhdu.data._coldefs._arrays[1][3], 33)
        assert_equal(tbhdu.data._coldefs.columns[1].array[3], 33)
        assert_equal(tbhdu.columns._arrays[1][3], 33)
        assert_equal(tbhdu.columns.data[1].array[3], 33)
        assert_equal(tbhdu.columns.data[0].array[3], 'JIM1')
        assert_equal(tbhdu.columns.data[2].array[3], 'A Note')
        assert_equal(tbhdu.columns.data[3].array[3].all(),
                         np.array([1., 2., 3., 4., 5.],dtype=np.float32).all())
        assert_equal(tbhdu.columns.data[4].array[3], True)

    def test_assign_multiple_rows_to_table(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu1 = pyfits.new_table(coldefs)

        counts = np.array([112, 134, 108, 117])
        names = np.array(['NGC5', 'NGC6', 'NGC7', 'NCG8'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[0, 1, 0, 0])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = pyfits.new_table(coldefs)
        tbhdu.data[0][3] = np.array([1., 2., 3., 4., 5.], dtype=np.float32)

        tbhdu2 = pyfits.new_table(tbhdu1.data, nrows=9)

        # Assign the 4 rows from the second table to rows 5 thru 8 of the
        # new table.  Note that the last row of the new table will still be
        # initialized to the default values.
        tbhdu2.data[4:] = tbhdu.data

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert_equal(id(tbhdu2.data._coldefs.columns[0].array),
                         id(tbhdu2.data._coldefs._arrays[0]))
        assert_equal(id(tbhdu2.data._coldefs.columns[0].array),
                         id(tbhdu2.columns.data[0].array))
        assert_equal(id(tbhdu2.data._coldefs.columns[0].array),
                         id(tbhdu2.columns._arrays[0]))

        assert_equal(tbhdu2.data[0][1], 312)
        assert_equal(tbhdu2.data._coldefs._arrays[1][0], 312)
        assert_equal(tbhdu2.data._coldefs.columns[1].array[0], 312)
        assert_equal(tbhdu2.columns._arrays[1][0], 312)
        assert_equal(tbhdu2.columns.data[1].array[0], 312)
        assert_equal(tbhdu2.columns.data[0].array[0], 'NGC1')
        assert_equal(tbhdu2.columns.data[2].array[0], '0.0')
        assert_equal(tbhdu2.columns.data[3].array[0].all(),
                         np.array([0., 0., 0., 0., 0.],dtype=np.float32).all())
        assert_equal(tbhdu2.columns.data[4].array[0], True)

        assert_equal(tbhdu2.data[4][1], 112)
        assert_equal(tbhdu2.data._coldefs._arrays[1][4], 112)
        assert_equal(tbhdu2.data._coldefs.columns[1].array[4], 112)
        assert_equal(tbhdu2.columns._arrays[1][4], 112)
        assert_equal(tbhdu2.columns.data[1].array[4], 112)
        assert_equal(tbhdu2.columns.data[0].array[4], 'NGC5')
        assert_equal(tbhdu2.columns.data[2].array[4], '0.0')
        assert_equal(tbhdu2.columns.data[3].array[4].all(),
                         np.array([1., 2., 3., 4., 5.],dtype=np.float32).all())
        assert_equal(tbhdu2.columns.data[4].array[4], False)

        assert_equal(tbhdu2.columns.data[1].array[8], 0)
        assert_equal(tbhdu2.columns.data[0].array[8], '0.0')
        assert_equal(tbhdu2.columns.data[2].array[8], '0.0')
        assert_equal(tbhdu2.columns.data[3].array[8].all(),
                         np.array([0., 0., 0., 0., 0.],dtype=np.float32).all())
        assert_equal(tbhdu2.columns.data[4].array[8], False)

    def test_verify_data_references(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = pyfits.new_table(coldefs)

        # Verify that original ColDefs object has independent Column
        # objects.
        assert_not_equal(id(coldefs.columns[0]), id(c1))

        # Verify that original ColDefs object has independent ndarray
        # objects.
        assert_not_equal(id(coldefs.columns[0].array), id(names))

        # Verify that original ColDefs object references the same data
        # object as the original Column object.
        assert_equal(id(coldefs.columns[0].array), id(c1.array))
        assert_equal(id(coldefs.columns[0].array), id(coldefs._arrays[0]))

        # Verify new HDU has an independent ColDefs object.
        assert_not_equal(id(coldefs), id(tbhdu.columns))

        # Verify new HDU has independent Column objects.
        assert_not_equal(id(coldefs.columns[0]), id(tbhdu.columns.data[0]))

        # Verify new HDU has independent ndarray objects.
        assert_not_equal(id(coldefs.columns[0].array),
                         id(tbhdu.columns.data[0].array))

        # Verify that both ColDefs objects in the HDU reference the same
        # Coldefs object.
        assert_equal(id(tbhdu.columns), id(tbhdu.data._coldefs))

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert_equal(id(tbhdu.data._coldefs.columns[0].array),
                         id(tbhdu.data._coldefs._arrays[0]))
        assert_equal(id(tbhdu.data._coldefs.columns[0].array),
                         id(tbhdu.columns.data[0].array))
        assert_equal(id(tbhdu.data._coldefs.columns[0].array),
                         id(tbhdu.columns._arrays[0]))

        tbhdu.writeto(self.temp('table1.fits'))

        t1 = pyfits.open(self.temp('table1.fits'))

        t1[1].data[0][1] = 213

        assert_equal(t1[1].data[0][1], 213)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 213)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 213)
        assert_equal(t1[1].columns._arrays[1][0], 213)
        assert_equal(t1[1].columns.data[1].array[0], 213)

        t1[1].data._coldefs._arrays[1][0] = 100

        assert_equal(t1[1].data[0][1], 100)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 100)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 100)
        assert_equal(t1[1].columns._arrays[1][0], 100)
        assert_equal(t1[1].columns.data[1].array[0], 100)

        t1[1].data._coldefs.columns[1].array[0] = 500
        assert_equal(t1[1].data[0][1], 500)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 500)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 500)
        assert_equal(t1[1].columns._arrays[1][0], 500)
        assert_equal(t1[1].columns.data[1].array[0], 500)

        t1[1].columns._arrays[1][0] = 600
        assert_equal(t1[1].data[0][1], 600)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 600)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 600)
        assert_equal(t1[1].columns._arrays[1][0], 600)
        assert_equal(t1[1].columns.data[1].array[0], 600)

        t1[1].columns.data[1].array[0] = 800
        assert_equal(t1[1].data[0][1], 800)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 800)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 800)
        assert_equal(t1[1].columns._arrays[1][0], 800)
        assert_equal(t1[1].columns.data[1].array[0], 800)

        t1.close()

    def test_new_table_with_ndarray(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu = pyfits.new_table(coldefs)

        tbhdu1 = pyfits.new_table(tbhdu.data.view(np.ndarray))

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert_equal(id(tbhdu1.data._coldefs.columns[0].array),
                         id(tbhdu1.data._coldefs._arrays[0]))
        assert_equal(id(tbhdu1.data._coldefs.columns[0].array),
                         id(tbhdu1.columns.data[0].array))
        assert_equal(id(tbhdu1.data._coldefs.columns[0].array),
                         id(tbhdu1.columns._arrays[0]))

        # Ensure I can change the value of one data element and it effects
        # all of the others.
        tbhdu1.data[0][1] = 213

        assert_equal(tbhdu1.data[0][1], 213)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 213)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 213)
        assert_equal(tbhdu1.columns._arrays[1][0], 213)
        assert_equal(tbhdu1.columns.data[1].array[0], 213)

        tbhdu1.data._coldefs._arrays[1][0] = 100

        assert_equal(tbhdu1.data[0][1], 100)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 100)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 100)
        assert_equal(tbhdu1.columns._arrays[1][0], 100)
        assert_equal(tbhdu1.columns.data[1].array[0], 100)

        tbhdu1.data._coldefs.columns[1].array[0] = 500
        assert_equal(tbhdu1.data[0][1], 500)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 500)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 500)
        assert_equal(tbhdu1.columns._arrays[1][0], 500)
        assert_equal(tbhdu1.columns.data[1].array[0], 500)

        tbhdu1.columns._arrays[1][0] = 600
        assert_equal(tbhdu1.data[0][1], 600)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 600)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 600)
        assert_equal(tbhdu1.columns._arrays[1][0], 600)
        assert_equal(tbhdu1.columns.data[1].array[0], 600)

        tbhdu1.columns.data[1].array[0] = 800
        assert_equal(tbhdu1.data[0][1], 800)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 800)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 800)
        assert_equal(tbhdu1.columns._arrays[1][0], 800)
        assert_equal(tbhdu1.columns.data[1].array[0], 800)

        tbhdu1.writeto(self.temp('table1.fits'))

        t1=pyfits.open(self.temp('table1.fits'))

        t1[1].data[0][1] = 213

        assert_equal(t1[1].data[0][1], 213)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 213)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 213)
        assert_equal(t1[1].columns._arrays[1][0], 213)
        assert_equal(t1[1].columns.data[1].array[0], 213)

        t1[1].data._coldefs._arrays[1][0] = 100

        assert_equal(t1[1].data[0][1], 100)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 100)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 100)
        assert_equal(t1[1].columns._arrays[1][0], 100)
        assert_equal(t1[1].columns.data[1].array[0], 100)

        t1[1].data._coldefs.columns[1].array[0] = 500
        assert_equal(t1[1].data[0][1], 500)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 500)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 500)
        assert_equal(t1[1].columns._arrays[1][0], 500)
        assert_equal(t1[1].columns.data[1].array[0], 500)

        t1[1].columns._arrays[1][0] = 600
        assert_equal(t1[1].data[0][1], 600)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 600)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 600)
        assert_equal(t1[1].columns._arrays[1][0], 600)
        assert_equal(t1[1].columns.data[1].array[0], 600)

        t1[1].columns.data[1].array[0] = 800
        assert_equal(t1[1].data[0][1], 800)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 800)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 800)
        assert_equal(t1[1].columns._arrays[1][0], 800)
        assert_equal(t1[1].columns.data[1].array[0], 800)

        t1.close()

    def test_new_table_with_fits_rec(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu=pyfits.new_table(coldefs)

        tbhdu.data[0][1] = 213

        assert_equal(tbhdu.data[0][1], 213)
        assert_equal(tbhdu.data._coldefs._arrays[1][0], 213)
        assert_equal(tbhdu.data._coldefs.columns[1].array[0], 213)
        assert_equal(tbhdu.columns._arrays[1][0], 213)
        assert_equal(tbhdu.columns.data[1].array[0], 213)

        tbhdu.data._coldefs._arrays[1][0] = 100

        assert_equal(tbhdu.data[0][1], 100)
        assert_equal(tbhdu.data._coldefs._arrays[1][0], 100)
        assert_equal(tbhdu.data._coldefs.columns[1].array[0], 100)
        assert_equal(tbhdu.columns._arrays[1][0], 100)
        assert_equal(tbhdu.columns.data[1].array[0], 100)

        tbhdu.data._coldefs.columns[1].array[0] = 500
        assert_equal(tbhdu.data[0][1], 500)
        assert_equal(tbhdu.data._coldefs._arrays[1][0], 500)
        assert_equal(tbhdu.data._coldefs.columns[1].array[0], 500)
        assert_equal(tbhdu.columns._arrays[1][0], 500)
        assert_equal(tbhdu.columns.data[1].array[0], 500)

        tbhdu.columns._arrays[1][0] = 600
        assert_equal(tbhdu.data[0][1], 600)
        assert_equal(tbhdu.data._coldefs._arrays[1][0], 600)
        assert_equal(tbhdu.data._coldefs.columns[1].array[0], 600)
        assert_equal(tbhdu.columns._arrays[1][0], 600)
        assert_equal(tbhdu.columns.data[1].array[0], 600)

        tbhdu.columns.data[1].array[0] = 800
        assert_equal(tbhdu.data[0][1], 800)
        assert_equal(tbhdu.data._coldefs._arrays[1][0], 800)
        assert_equal(tbhdu.data._coldefs.columns[1].array[0], 800)
        assert_equal(tbhdu.columns._arrays[1][0], 800)
        assert_equal(tbhdu.columns.data[1].array[0], 800)

        tbhdu.columns.data[1].array[0] = 312

        tbhdu.writeto(self.temp('table1.fits'))

        t1=pyfits.open(self.temp('table1.fits'))

        t1[1].data[0][1] = 1
        fr = t1[1].data
        assert_equal(t1[1].data[0][1], 1)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 1)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 1)
        assert_equal(t1[1].columns._arrays[1][0], 1)
        assert_equal(t1[1].columns.data[1].array[0], 1)
        assert_equal(fr[0][1], 1)
        assert_equal(fr._coldefs._arrays[1][0], 1)
        assert_equal(fr._coldefs.columns[1].array[0], 1)

        fr._coldefs.columns[1].array[0] = 312

        tbhdu1 = pyfits.new_table(fr)
        #tbhdu1 = pyfits.new_table(t1[1].data)

        i = 0
        for row in tbhdu1.data:
            for j in range(0,len(row)):
                if isinstance(row[j], np.ndarray):
                    assert_equal(row[j].all(), tbhdu.data[i][j].all())
                else:
                    assert_equal(row[j], tbhdu.data[i][j])
            i = i + 1

        tbhdu1.data[0][1] = 213

        assert_equal(t1[1].data[0][1], 312)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 312)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 312)
        assert_equal(t1[1].columns._arrays[1][0], 312)
        assert_equal(t1[1].columns.data[1].array[0], 312)
        assert_equal(fr[0][1], 312)
        assert_equal(fr._coldefs._arrays[1][0], 312)
        assert_equal(fr._coldefs.columns[1].array[0], 312)
        assert_equal(tbhdu1.data[0][1], 213)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 213)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 213)
        assert_equal(tbhdu1.columns._arrays[1][0], 213)
        assert_equal(tbhdu1.columns.data[1].array[0], 213)

        t1[1].data[0][1] = 10

        assert_equal(t1[1].data[0][1], 10)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 10)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 10)
        assert_equal(t1[1].columns._arrays[1][0], 10)
        assert_equal(t1[1].columns.data[1].array[0], 10)
        assert_equal(fr[0][1], 10)
        assert_equal(fr._coldefs._arrays[1][0], 10)
        assert_equal(fr._coldefs.columns[1].array[0], 10)
        assert_equal(tbhdu1.data[0][1], 213)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 213)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 213)
        assert_equal(tbhdu1.columns._arrays[1][0], 213)
        assert_equal(tbhdu1.columns.data[1].array[0], 213)

        tbhdu1.data._coldefs._arrays[1][0] = 666

        assert_equal(t1[1].data[0][1], 10)
        assert_equal(t1[1].data._coldefs._arrays[1][0], 10)
        assert_equal(t1[1].data._coldefs.columns[1].array[0], 10)
        assert_equal(t1[1].columns._arrays[1][0], 10)
        assert_equal(t1[1].columns.data[1].array[0], 10)
        assert_equal(fr[0][1], 10)
        assert_equal(fr._coldefs._arrays[1][0], 10)
        assert_equal(fr._coldefs.columns[1].array[0], 10)
        assert_equal(tbhdu1.data[0][1], 666)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 666)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 666)
        assert_equal(tbhdu1.columns._arrays[1][0], 666)
        assert_equal(tbhdu1.columns.data[1].array[0], 666)

        t1.close()

    def test_bin_table_hdu_constructor(self):
        counts = np.array([312, 334, 308, 317])
        names = np.array(['NGC1', 'NGC2', 'NGC3', 'NCG4'])
        c1 = pyfits.Column(name='target', format='10A', array=names)
        c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)
        c3 = pyfits.Column(name='notes', format='A10')
        c4 = pyfits.Column(name='spectrum', format='5E')
        c5 = pyfits.Column(name='flag', format='L', array=[1, 0, 1, 1])
        coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])

        tbhdu1=pyfits.new_table(coldefs)

        hdu = pyfits.BinTableHDU(tbhdu1.data)

        # Verify that all ndarray objects within the HDU reference the
        # same ndarray.
        assert_equal(id(hdu.data._coldefs.columns[0].array),
                         id(hdu.data._coldefs._arrays[0]))
        assert_equal(id(hdu.data._coldefs.columns[0].array),
                         id(hdu.columns.data[0].array))
        assert_equal(id(hdu.data._coldefs.columns[0].array),
                         id(hdu.columns._arrays[0]))

        # Verify that the references in the original HDU are the same as the
        # references in the new HDU.
        assert_equal(id(tbhdu1.data._coldefs.columns[0].array),
                         id(hdu.data._coldefs._arrays[0]))


        # Verify that a change in the new HDU is reflected in both the new
        # and original HDU.

        hdu.data[0][1] = 213

        assert_equal(hdu.data[0][1], 213)
        assert_equal(hdu.data._coldefs._arrays[1][0], 213)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 213)
        assert_equal(hdu.columns._arrays[1][0], 213)
        assert_equal(hdu.columns.data[1].array[0], 213)
        assert_equal(tbhdu1.data[0][1], 213)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 213)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 213)
        assert_equal(tbhdu1.columns._arrays[1][0], 213)
        assert_equal(tbhdu1.columns.data[1].array[0], 213)

        hdu.data._coldefs._arrays[1][0] = 100

        assert_equal(hdu.data[0][1], 100)
        assert_equal(hdu.data._coldefs._arrays[1][0], 100)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 100)
        assert_equal(hdu.columns._arrays[1][0], 100)
        assert_equal(hdu.columns.data[1].array[0], 100)
        assert_equal(tbhdu1.data[0][1], 100)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 100)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 100)
        assert_equal(tbhdu1.columns._arrays[1][0], 100)
        assert_equal(tbhdu1.columns.data[1].array[0], 100)

        hdu.data._coldefs.columns[1].array[0] = 500
        assert_equal(hdu.data[0][1], 500)
        assert_equal(hdu.data._coldefs._arrays[1][0], 500)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 500)
        assert_equal(hdu.columns._arrays[1][0], 500)
        assert_equal(hdu.columns.data[1].array[0], 500)
        assert_equal(tbhdu1.data[0][1], 500)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 500)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 500)
        assert_equal(tbhdu1.columns._arrays[1][0], 500)
        assert_equal(tbhdu1.columns.data[1].array[0], 500)

        hdu.columns._arrays[1][0] = 600
        assert_equal(hdu.data[0][1], 600)
        assert_equal(hdu.data._coldefs._arrays[1][0], 600)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 600)
        assert_equal(hdu.columns._arrays[1][0], 600)
        assert_equal(hdu.columns.data[1].array[0], 600)
        assert_equal(tbhdu1.data[0][1], 600)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 600)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 600)
        assert_equal(tbhdu1.columns._arrays[1][0], 600)
        assert_equal(tbhdu1.columns.data[1].array[0], 600)

        hdu.columns.data[1].array[0] = 800
        assert_equal(hdu.data[0][1], 800)
        assert_equal(hdu.data._coldefs._arrays[1][0], 800)
        assert_equal(hdu.data._coldefs.columns[1].array[0], 800)
        assert_equal(hdu.columns._arrays[1][0], 800)
        assert_equal(hdu.columns.data[1].array[0], 800)
        assert_equal(tbhdu1.data[0][1], 800)
        assert_equal(tbhdu1.data._coldefs._arrays[1][0], 800)
        assert_equal(tbhdu1.data._coldefs.columns[1].array[0], 800)
        assert_equal(tbhdu1.columns._arrays[1][0], 800)
        assert_equal(tbhdu1.columns.data[1].array[0], 800)

    def test_constructor_name_arg(self):
        """testConstructorNameArg

        Passing name='...' to the BinTableHDU and TableHDU constructors
        should set the .name attribute and 'EXTNAME' header keyword, and
        override any name in an existing 'EXTNAME' value.
        """

        for hducls in [pyfits.BinTableHDU, pyfits.TableHDU]:
            # First test some default assumptions
            hdu = hducls()
            assert_equal(hdu.name, '')
            assert_true('EXTNAME' not in hdu.header)
            hdu.name = 'FOO'
            assert_equal(hdu.name, 'FOO')
            assert_equal(hdu.header['EXTNAME'], 'FOO')

            # Passing name to constructor
            hdu = hducls(name='FOO')
            assert_equal(hdu.name, 'FOO')
            assert_equal(hdu.header['EXTNAME'], 'FOO')

            # And overriding a header with a different extname
            hdr = pyfits.Header()
            hdr.update('EXTNAME', 'EVENTS')
            hdu = hducls(header=hdr, name='FOO')
            assert_equal(hdu.name, 'FOO')
            assert_equal(hdu.header['EXTNAME'], 'FOO')


    def test_bin_table_with_logical_array(self):
        c1 = pyfits.Column(name='flag', format='2L',
                           array=[[True, False], [False, True]])
        coldefs = pyfits.ColDefs([c1])

        tbhdu1 = pyfits.new_table(coldefs)

        assert_equal(tbhdu1.data.field('flag')[0].all(),
                         np.array([True, False],
                                  dtype = np.bool).all())
        assert_equal(tbhdu1.data.field('flag')[1].all(),
                         np.array([False, True],
                                  dtype = np.bool).all())

        tbhdu = pyfits.new_table(tbhdu1.data)

        assert_equal(tbhdu.data.field('flag')[0].all(),
                         np.array([True, False],
                                  dtype = np.bool).all())
        assert_equal(tbhdu.data.field('flag')[1].all(),
                         np.array([False, True],
                                  dtype = np.bool).all())

    def test_variable_length_table_format_pd_from_object_array(self):
        a = np.array([np.array([7.2e-20, 7.3e-20]), np.array([0.0]),
                      np.array([0.0])], 'O')
        acol = pyfits.Column(name='testa', format='PD()', array=a)
        tbhdu = pyfits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        tbhdu1 = pyfits.open(self.temp('newtable.fits'))

        for j in range(0,3):
            for i in range(0,len(a[j])):
                assert_equal(tbhdu1[1].data.field(0)[j][i], a[j][i])

        tbhdu1.close()

    def test_variable_length_table_format_pd_from_list(self):
        a = [np.array([7.2e-20,7.3e-20]),np.array([0.0]),np.array([0.0])]
        acol = pyfits.Column(name='testa',format='PD()',array=a)
        tbhdu = pyfits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        tbhdu1 = pyfits.open(self.temp('newtable.fits'))

        for j in range(0,3):
            for i in range(0,len(a[j])):
                assert_equal(tbhdu1[1].data.field(0)[j][i], a[j][i])

        tbhdu1.close()

    def test_variable_length_table_format_pa_from_object_array(self):
        a = np.array([np.array(['a', 'b', 'c']), np.array(['d', 'e']),
                      np.array(['f'])], 'O')
        acol = pyfits.Column(name='testa', format='PA()', array=a)
        tbhdu = pyfits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        hdul = pyfits.open(self.temp('newtable.fits'))

        for j in range(0,3):
            for i in range(0,len(a[j])):
                assert_equal(hdul[1].data.field(0)[j][i], a[j][i])

        hdul.close()

    def test_variable_length_table_format_pa_from_list(self):
        a = ['a', 'ab', 'abc']
        acol = pyfits.Column(name='testa', format='PA()', array=a)
        tbhdu = pyfits.new_table([acol])
        tbhdu.writeto(self.temp('newtable.fits'))
        hdul = pyfits.open(self.temp('newtable.fits'))

        for j in range(0,3):
            for i in range(0,len(a[j])):
                assert_equal(hdul[1].data.field(0)[j][i], a[j][i])

        hdul.close()

    def test_fits_rec_column_access(self):
        t=pyfits.open(self.data('table.fits'))
        tbdata = t[1].data
        assert_equal(tbdata.V_mag.all(), tbdata.field('V_mag').all())
        assert_equal(tbdata.V_mag.all(), tbdata['V_mag'].all())

        t.close()

    def test_table_with_zero_width_column(self):
        hdul = pyfits.open(self.data('zerowidth.fits'))
        tbhdu = hdul[2] # This HDU contains a zero-width column 'ORBPARM'
        assert_true('ORBPARM' in tbhdu.columns.names)
        # The ORBPARM column should not be in the data, though the data should
        # be readable
        assert_true('ORBPARM' not in tbhdu.data.names)
        # Verify that some of the data columns are still correctly accessible
        # by name
        assert_true(comparefloats(
            tbhdu.data[0]['STABXYZ'],
            np.array([499.85566663, -1317.99231554, -735.18866164],
                     dtype=np.float64)))
        assert_equal(tbhdu.data[0]['NOSTA'], 1)
        assert_equal(tbhdu.data[0]['MNTSTA'], 0)
        hdul.writeto(self.temp('newtable.fits'))
        hdul.close()
        hdul = pyfits.open(self.temp('newtable.fits'))
        tbhdu = hdul[2]
        # Verify that the previous tests still hold after writing
        assert_true('ORBPARM' in tbhdu.columns.names)
        assert_true('ORBPARM' not in tbhdu.data.names)
        assert_true(comparefloats(
            tbhdu.data[0]['STABXYZ'],
            np.array([499.85566663, -1317.99231554, -735.18866164],
                     dtype=np.float64)))
        assert_equal(tbhdu.data[0]['NOSTA'], 1)
        assert_equal(tbhdu.data[0]['MNTSTA'], 0)
        hdul.close()

    def test_string_column_padding(self):
        a = ['img1', 'img2', 'img3a', 'p']
        s = 'img1\x00\x00\x00\x00\x00\x00' \
            'img2\x00\x00\x00\x00\x00\x00' \
            'img3a\x00\x00\x00\x00\x00' \
            'p\x00\x00\x00\x00\x00\x00\x00\x00\x00'

        acol = pyfits.Column(name='MEMNAME', format='A10',
                             array=chararray.array(a))
        ahdu = pyfits.new_table([acol])
        assert_equal(ahdu.data.tostring().decode('raw-unicode-escape'), s)

        ahdu = pyfits.new_table([acol], tbtype='TableHDU')
        assert_equal(ahdu.data.tostring().decode('raw-unicode-escape'),
                         s.replace('\x00', ' '))

    def test_multi_dimensional_columns(self):
        """
        Tests the multidimensional column implementation with both numeric
        arrays and string arrays.
        """

        data = np.rec.array(
            [([0, 1, 2, 3, 4, 5], 'row1' * 2),
             ([6, 7, 8, 9, 0, 1], 'row2' * 2),
             ([2, 3, 4, 5, 6, 7], 'row3' * 2)], formats='6i4,a8')

        thdu = pyfits.new_table(data)
        # Modify the TDIM fields to my own specification
        thdu.header.update('TDIM1', '(2,3)')
        thdu.header.update('TDIM2', '(4,2)')

        thdu.writeto(self.temp('newtable.fits'))

        hdul = pyfits.open(self.temp('newtable.fits'))
        thdu = hdul[1]

        c1 = thdu.data.field(0)
        c2 = thdu.data.field(1)

        hdul.close()

        assert_equal(c1.shape, (3, 3, 2))
        assert_equal(c2.shape, (3, 2))
        assert_true((c1 == np.array([[[0, 1], [2, 3], [4, 5]],
                                     [[6, 7], [8, 9], [0, 1]],
                                     [[2, 3], [4, 5], [6, 7]]])).all())
        assert_true((c2 == np.array([['row1', 'row1'],
                                     ['row2', 'row2'],
                                     ['row3', 'row3']])).all())

        # Test setting the TDIMn header based on the column data
        data = np.zeros(3, dtype=[('x', 'f4'), ('s', 'S5', 4)])
        data['x'] = 1, 2, 3
        data['s'] = 'ok'
        pyfits.writeto(self.temp('newtable.fits'), data, clobber=True)

        t = pyfits.getdata(self.temp('newtable.fits'))

        assert_equal(t.field(1).dtype.str[-1], '5')
        assert_equal(t.field(1).shape, (3, 4))

        # Like the previous test, but with an extra dimension (a bit more
        # complicated)
        data = np.zeros(3, dtype=[('x', 'f4'), ('s', 'S5', (4, 3))])
        data['x'] = 1, 2, 3
        data['s'] = 'ok'
        pyfits.writeto(self.temp('newtable.fits'), data, clobber=True)

        t = pyfits.getdata(self.temp('newtable.fits'))

        assert_equal(t.field(1).dtype.str[-1], '5')
        assert_equal(t.field(1).shape, (3, 4, 3))

    def test_slicing(self):
        """Regression test for #52."""

        f = pyfits.open(self.data('table.fits'))
        data = f[1].data
        targets = data.field('target')
        s = data[:]
        assert_true((s.field('target') == targets).all())
        for n in range(len(targets) + 2):
            s = data[:n]
            assert_true((s.field('target') == targets[:n]).all())
            s = data[n:]
            assert_true((s.field('target') == targets[n:]).all())
        s = data[::2]
        assert_true((s.field('target') == targets[::2]).all())
        s = data[::-1]
        assert_true((s.field('target') == targets[::-1]).all())

    def test_array_slicing(self):
        """Regression test for #55."""

        f = pyfits.open(self.data('table.fits'))
        data = f[1].data
        s1 = data[data['target'] == 'NGC1001']
        s2 = data[np.where(data['target'] == 'NGC1001')]
        s3 = data[[0]]
        s4 = data[:1]
        for s in [s1, s2, s3, s4]:
            assert_true(isinstance(s, pyfits.FITS_rec))
        assert_true((s1 == s2).all())
        assert_true((s2 == s3).all())
        assert_true((s3 == s4).all())

    def test_dump_load_round_trip(self):
        """
        A simple test of the dump/load methods; dump the data, column, and
        header files and try to reload the table from them.
        """

        hdul = pyfits.open(self.data('table.fits'))
        tbhdu = hdul[1]
        datafile = self.temp('data.txt')
        cdfile = self.temp('coldefs.txt')
        hfile = self.temp('header.txt')

        tbhdu.dump(datafile, cdfile, hfile)

        new_tbhdu = pyfits.BinTableHDU.load(datafile, cdfile, hfile)

        assert_true(comparerecords(tbhdu.data, new_tbhdu.data))

        # Double check that the headers are equivalent
        assert_equal(str(tbhdu.header.ascard), str(new_tbhdu.header.ascard))

    def test_load_guess_format(self):
        """
        Tests loading a table dump with no supplied coldefs or header, so that
        the table format has to be guessed at.  There is of course no exact
        science to this; the table that's produced simply uses sensible guesses
        for that format.  Ideally this should never have to be used.
        """

        # Create a table containing a variety of data types.
        a0 = np.array([False, True, False], dtype=np.bool)
        c0 = pyfits.Column(name='c0', format='L', array=a0)

        # Format X currently not supported by the format
        #a1 = np.array([[0], [1], [0]], dtype=np.uint8)
        #c1 = pyfits.Column(name='c1', format='X', array=a1)

        a2 = np.array([1, 128, 255], dtype=np.uint8)
        c2 = pyfits.Column(name='c2', format='B', array=a2)
        a3 = np.array([-30000, 1, 256], dtype=np.int16)
        c3 = pyfits.Column(name='c3', format='I', array=a3)
        a4 = np.array([-123123123, 1234, 123123123], dtype=np.int32)
        c4 = pyfits.Column(name='c4', format='J', array=a4)
        a5 = np.array(['a', 'abc', 'ab'])
        c5 = pyfits.Column(name='c5', format='A3', array=a5)
        a6 = np.array([1.1, 2.2, 3.3], dtype=np.float64)
        c6 = pyfits.Column(name='c6', format='D', array=a6)
        a7 = np.array([1.1+2.2j, 3.3+4.4j, 5.5+6.6j], dtype=np.complex128)
        c7 = pyfits.Column(name='c7', format='M', array=a7)
        a8 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.int32)
        c8 = pyfits.Column(name='c8', format='PJ()', array=a8)

        tbhdu = pyfits.new_table([c0, c2, c3, c4, c5, c6, c7, c8])

        datafile = self.temp('data.txt')
        tbhdu.dump(datafile)

        new_tbhdu = pyfits.BinTableHDU.load(datafile)

        # In this particular case the record data at least should be equivalent
        assert_true(comparerecords(tbhdu.data, new_tbhdu.data))

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

        c1 = pyfits.Column(name='names', format='I', array=[1])
        c2 = pyfits.Column(name='formats', format='I', array=[2])
        c3 = pyfits.Column(name='other', format='I', array=[3])

        t = pyfits.new_table([c1, c2, c3])
        assert_equal(t.data.names, ['names', 'formats', 'other'])
        assert_equal(t.data.formats, ['I'] * 3)
        assert_true((t.data['names'] == [1]).all())
        assert_true((t.data['formats'] == [2]).all())
        assert_true((t.data.other == [3]).all())
