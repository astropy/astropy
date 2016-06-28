# Licensed under a 3-clause BSD style license - see PYFITS.rst

import sys

import numpy as np

from ....io import fits
from . import FitsTestCase


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
                    raise ValueError('field name {} not found in array 1'.format(n2))

        if verbose:
            sys.stdout.write("    testing field: '{}'\n".format(n2))
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
                    sys.stdout.write(
                        '\n        {} elements in field {} differ\n'.format(
                            w.size, n2))
            else:
                if verbose:
                    sys.stdout.write('OK\n')

    if nfail == 0:
        if verbose:
            sys.stdout.write('All tests passed\n')
        return True
    else:
        if verbose:
            sys.stdout.write('{} differences found\n'.format(nfail))
        return False


def get_test_data(verbose=False):
    st = np.zeros(3, [('f1', 'i4'), ('f2', 'S6'), ('f3', '>2f8')])

    np.random.seed(35)
    st['f1'] = [1, 3, 5]
    st['f2'] = ['hello', 'world', 'byebye']
    st['f3'] = np.random.random(st['f3'].shape)

    return st


class TestStructured(FitsTestCase):
    def test_structured(self):
        fname = self.data('stddata.fits')

        data1, h1 = fits.getdata(fname, ext=1, header=True)
        data2, h2 = fits.getdata(fname, ext=2, header=True)

        st = get_test_data()

        outfile = self.temp('test.fits')
        fits.writeto(outfile, data1, overwrite=True)
        fits.append(outfile, data2)

        fits.append(outfile, st)
        assert st.dtype.isnative
        assert np.all(st['f1'] == [1, 3, 5])

        data1check, h1check = fits.getdata(outfile, ext=1, header=True)
        data2check, h2check = fits.getdata(outfile, ext=2, header=True)
        stcheck, sthcheck = fits.getdata(outfile, ext=3, header=True)

        assert compare_arrays(data1, data1check, verbose=True)
        assert compare_arrays(data2, data2check, verbose=True)
        assert compare_arrays(st, stcheck, verbose=True)

        # try reading with view
        dataviewcheck, hviewcheck = fits.getdata(outfile, ext=2, header=True,
                                                 view=np.ndarray)
        assert compare_arrays(data2, dataviewcheck, verbose=True)
