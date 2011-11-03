from __future__ import division # confidence high

import sys

import numpy as np

import pyfits
from pyfits.tests import PyfitsTestCase


def compare_arrays(arr1in, arr2in, verbose=False):
    """
    Compare the values field-by-field in two sets of numpy arrays or
    recarrays.
    """

    arr1 = arr1in.view(np.ndarray)
    arr2 = arr2in.view(np.ndarray)

    nfail = 0
    for n2 in arr2.dtype.names:
        n1 = n2
        if n1 not in arr1.dtype.names:
            n1 = n1.lower()
            if n1 not in arr1.dtype.names:
                n1 = n1.upper()
                if n1 not in arr1.dtype.names:
                    raise ValueError('field name %s not found in array 1' % n2)

        if verbose:
            sys.stdout.write("    testing field: '%s'\n" % n2)
            sys.stdout.write('        shape...........')
        if arr2[n2].shape != arr1[n1].shape:
            nfail += 1
            if verbose:
                sys.stdout.write('shapes differ\n')
        else:
            if verbose:
                sys.stdout.write('OK\n')
                sys.stdout.write('        elements........')
            w, = np.where(arr1[n1].ravel() != arr2[n2].ravel())
            if w.size > 0:
                nfail += 1
                if verbose:
                    sys.stdout.write('\n        '+\
                            '%s elements in field %s differ\n' % (w.size,n2))
            else:
                if verbose:
                    sys.stdout.write('OK\n')

    if nfail == 0:
        if verbose:
            sys.stdout.write('All tests passed\n')
        return True
    else:
        if verbose:
            sys.stdout.write('%d differences found\n' % nfail)
        return False


def get_test_data(verbose=False):
    st = np.zeros(3, [('f1', 'i4'), ('f2', 'S6'), ('f3', '>2f8')])

    np.random.seed(35)
    st['f1'] = [1, 3, 5]
    st['f2'] = ['hello', 'world', 'byebye']
    st['f3'] = np.random.random(st['f3'].shape)
    print st.dtype.descr
    print st

    return st


class TestStructured(PyfitsTestCase):
    def test_structured(self):
        fname = self.data('stddata.fits')

        print 'Reading from ', fname
        data1, h1 = pyfits.getdata(fname, ext=1, header=True)
        data2, h2 = pyfits.getdata(fname, ext=2, header=True)

        st = get_test_data()

        outfile = self.temp('test.fits')
        print 'Writing to file data1:', outfile
        pyfits.writeto(outfile, data1, clobber=True)
        print 'Appending to file: data2', outfile
        pyfits.append(outfile, data2)

        print 'Appending to file: st', outfile
        pyfits.append(outfile, st)
        print st.dtype.descr
        print st
        assert st.dtype.isnative
        assert np.all(st['f1'] == [1,3,5])

        print 'Reading data back'
        data1check, h1check = pyfits.getdata(outfile, ext=1, header=True)
        data2check, h2check = pyfits.getdata(outfile, ext=2, header=True)
        stcheck, sthcheck = pyfits.getdata(outfile, ext=3, header=True)

        if not compare_arrays(data1, data1check, verbose=True):
            raise ValueError('Fail')
        if not compare_arrays(data2, data2check, verbose=True):
            raise ValueError('Fail')
        print st, stcheck
        if not compare_arrays(st, stcheck, verbose=True):
            raise ValueError('Fail')

        # try reading with view
        print 'Reading with ndarray view'
        dataviewcheck, hviewcheck = pyfits.getdata(outfile, ext=2, header=True,
                                                   view=np.ndarray)
        if not compare_arrays(data2, dataviewcheck, verbose=True):
            raise ValueError('Fail')
