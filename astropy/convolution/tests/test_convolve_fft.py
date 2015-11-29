# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import itertools

import numpy as np
from numpy.testing import assert_array_almost_equal_nulp, assert_allclose

from ..convolve import convolve_fft
from ...tests.helper import pytest


VALID_DTYPES = []
for dtype_array in ['>f4', '<f4', '>f8', '<f8']:
    for dtype_kernel in ['>f4', '<f4', '>f8', '<f8']:
        VALID_DTYPES.append((dtype_array, dtype_kernel))

BOUNDARY_OPTIONS = [None, 'fill', 'wrap']

"""
What does convolution mean?  We use the 'same size' assumption here (i.e.,
you expect an array of the exact same size as the one you put in)
Convolving any array with a kernel that is [1] should result in the same array returned
Working example array: [1, 2, 3, 4, 5]
Convolved with [1] = [1, 2, 3, 4, 5]
Convolved with [1, 1] = [1, 3, 5, 7, 9] THIS IS NOT CONSISTENT!
Convolved with [1, 0] = [1, 2, 3, 4, 5]
Convolved with [0, 1] = [0, 1, 2, 3, 4]
"""

# NOTE: use_numpy_fft is redundant if you don't have FFTW installed
option_names = ('boundary', 'interpolate_nan', 'normalize_kernel',
                'ignore_edge_zeros')
options = list(itertools.product(BOUNDARY_OPTIONS, (True, False),
                                 (True, False), (True, False)))


class TestConvolve1D(object):

    @pytest.mark.parametrize(option_names, options)
    def test_unity_1_none(self, boundary, interpolate_nan, normalize_kernel,
                          ignore_edge_zeros):
        '''
        Test that a unit kernel with a single element returns the same array
        '''

        x = np.array([1., 2., 3.], dtype='float64')

        y = np.array([1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        assert_array_almost_equal_nulp(z, x, 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3(self, boundary, interpolate_nan, normalize_kernel,
                     ignore_edge_zeros):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None).
        '''

        x = np.array([1., 2., 3.], dtype='float64')

        y = np.array([0., 1., 0.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        assert_array_almost_equal_nulp(z, x, 10)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3(self, boundary, interpolate_nan, normalize_kernel,
                       ignore_edge_zeros):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements
        '''

        x = np.array([1., 0., 3.], dtype='float64')

        y = np.array([1., 1., 1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        answer_dict = {
            'sum': np.array([1., 4., 3.], dtype='float64'),
            'sum_zeros': np.array([1., 4., 3.], dtype='float64'),
            'sum_nozeros': np.array([1., 4., 3.], dtype='float64'),
            'average': np.array([1 / 3., 4 / 3., 1.], dtype='float64'),
            'sum_wrap': np.array([4., 4., 4.], dtype='float64'),
            'average_wrap': np.array([4 / 3., 4 / 3., 4 / 3.], dtype='float64'),
            'average_zeros': np.array([1 / 3., 4 / 3., 1.], dtype='float64'),
            'average_nozeros': np.array([0.5, 4 / 3., 1.5], dtype='float64'),
        }

        if normalize_kernel:
            answer_key = 'average'
        else:
            answer_key = 'sum'

        if boundary == 'wrap':
            answer_key += '_wrap'
        elif ignore_edge_zeros:
            answer_key += '_nozeros'
        else:
            # average = average_zeros; sum = sum_zeros
            answer_key += '_zeros'

        assert_array_almost_equal_nulp(z, answer_dict[answer_key], 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3_withnan(self, boundary, interpolate_nan, normalize_kernel,
                             ignore_edge_zeros):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='float64')

        y = np.array([0., 1., 0.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        # for whatever reason, numpy's fft has very limited precision, and
        # the comparison fails unless you cast the float64 to a float16
        if hasattr(np, 'float16'):
            assert_array_almost_equal_nulp(
                np.asarray(z, dtype=np.float16),
                np.array([1., 0., 3.], dtype=np.float16), 10)
        # ASSERT equality to better than 16 bit but worse than 32 bit precision
        assert np.all(np.abs(z - np.array([1., 0., 3.])) < 1e-14)

    inputs = (np.array([1., np.nan, 3.], dtype='float64'),
              np.array([1., np.inf, 3.], dtype='float64'))
    outputs = (np.array([1., 0., 3.], dtype='float64'),
               np.array([1., 0., 3.], dtype='float64'))
    options_unity1withnan = list(itertools.product(BOUNDARY_OPTIONS,
                                                   (True, False),
                                                   (True, False),
                                                   (True, False),
                                                   inputs, outputs))

    @pytest.mark.parametrize(option_names + ('inval', 'outval'),
                             options_unity1withnan)
    def test_unity_1_withnan(self, boundary, interpolate_nan, normalize_kernel,
                             ignore_edge_zeros, inval, outval):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        x = inval

        y = np.array([1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        # for whatever reason, numpy's fft has very limited precision, and
        # the comparison fails unless you cast the float64 to a float16
        if hasattr(np, 'float16'):
            assert_array_almost_equal_nulp(np.asarray(z, dtype=np.float16),
                                           np.array([1., 0., 3.], dtype=np.float16), 10)
        assert np.all(np.abs(z - outval) < 1e-14)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3_withnan(self, boundary, interpolate_nan,
                               normalize_kernel, ignore_edge_zeros):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements. This version includes a NaN
        value in the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='float64')

        y = np.array([1., 1., 1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        answer_dict = {
            'sum': np.array([1., 4., 3.], dtype='float64'),
            'sum_nozeros': np.array([1., 4., 3.], dtype='float64'),
            'sum_zeros': np.array([1., 4., 3.], dtype='float64'),
            'sum_zeros_noignan': np.array([1., 4., 3.], dtype='float64'),
            'sum_nozeros_noignan': np.array([1., 4., 3.], dtype='float64'),
            'average': np.array([1., 2., 3.], dtype='float64'),
            'sum_wrap': np.array([4., 4., 4.], dtype='float64'),
            'sum_wrap_noignan': np.array([4., 4., 4.], dtype='float64'),
            'average_wrap': np.array([(1 + 3) / 2., 2., 2.], dtype='float64'),
            'average_wrap_noignan': np.array([4 / 3., 4 / 3., 4 / 3.], dtype='float64'),
            'average_nozeros': np.array([1, 2, 3], dtype='float64'),
            'average_nozeros_noignan': np.array([1 / 2., 4 / 3., 3 / 2.], dtype='float64'),
            'average_zeros': np.array([1 / 2., 4 / 2., 3 / 2.], dtype='float64'),
            'average_zeros_noignan': np.array([1 / 3., 4 / 3., 3 / 3.], dtype='float64'),
        }

        if normalize_kernel:
            answer_key = 'average'
        else:
            answer_key = 'sum'

        if boundary == 'wrap':
            answer_key += '_wrap'
        elif ignore_edge_zeros:
            answer_key += '_nozeros'
        else:
            # average = average_zeros; sum = sum_zeros
            answer_key += '_zeros'

        if not interpolate_nan:
            answer_key += '_noignan'

        assert_array_almost_equal_nulp(z, answer_dict[answer_key], 10)

    def test_masked_array(self):
        """
        Check whether convolve_fft works with masked arrays.
        """
        # Test masked array
        array = np.array([1., np.nan, 3.], dtype='float64')
        kernel = np.array([1, 1, 1])
        masked_array = np.ma.masked_array(array, mask=[0, 1, 0])
        result = convolve_fft(masked_array, kernel)
        assert_allclose(result, [1, 4, 3], atol=1E-10)

        # Test masked kernel
        array = np.array([1., np.nan, 3.], dtype='float64')
        kernel = np.array([1, 1, 1])
        masked_array = np.ma.masked_array(array, mask=[0, 1, 0])
        result = convolve_fft(kernel, masked_array)
        assert_allclose(result, [1, 4, 3], atol=1E-10)

    def test_normalize_function(self):
        """
        Check if convolve_fft works when passing a normalize function.
        """
        array = [1, 2, 3]
        kernel = [3, 3, 3]
        result = convolve_fft(array, kernel, normalize_kernel=np.max)
        assert_allclose(result, [3, 6, 5], atol=1E-10)


class TestConvolve2D(object):

    @pytest.mark.parametrize(option_names, options)
    def test_unity_1x1_none(self, boundary, interpolate_nan, normalize_kernel,
                            ignore_edge_zeros):
        '''
        Test that a 1x1 unit kernel returns the same array
        '''

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype='float64')

        y = np.array([[1.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        assert_array_almost_equal_nulp(z, x, 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3x3(self, boundary, interpolate_nan, normalize_kernel,
                       ignore_edge_zeros):
        '''
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None).
        '''

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype='float64')

        y = np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 0.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        assert_array_almost_equal_nulp(z, x, 10)
        # assert np.all( np.abs(z-x) < np.spacing(np.where(z>x, z, x))*2 )

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3x3(self, boundary, interpolate_nan, normalize_kernel,
                         ignore_edge_zeros):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel.
        '''

        x = np.array([[0., 0., 3.],
                      [1., 0., 0.],
                      [0., 2., 0.]], dtype='float64')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros)

        w = np.array([[4., 6., 4.],
                      [6., 9., 6.],
                      [4., 6., 4.]], dtype='float64')
        answer_dict = {
            'sum': np.array([[1., 4., 3.],
                             [3., 6., 5.],
                             [3., 3., 2.]], dtype='float64'),
            'sum_wrap': np.array([[6., 6., 6.],
                                  [6., 6., 6.],
                                  [6., 6., 6.]], dtype='float64'),
        }
        answer_dict['average'] = answer_dict['sum'] / w
        answer_dict['average_wrap'] = answer_dict['sum_wrap'] / 9.
        answer_dict['average_withzeros'] = answer_dict['sum'] / 9.
        answer_dict['sum_withzeros'] = answer_dict['sum']

        if normalize_kernel:
            answer_key = 'average'
        else:
            answer_key = 'sum'

        if boundary == 'wrap':
            answer_key += '_wrap'
        elif not ignore_edge_zeros:
            answer_key += '_withzeros'

        a = answer_dict[answer_key]
        # for reasons unknown, the Windows FFT returns an answer for the [0, 0]
        # component that is EXACTLY 10*np.spacing
        assert np.all(np.abs(z - a) <= np.spacing(np.where(z > a, z, a)) * 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3x3_withnan(self, boundary, interpolate_nan,
                               normalize_kernel, ignore_edge_zeros):
        '''
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None). This version includes a NaN value in the original
        array.
        '''

        x = np.array([[1., 2., 3.],
                      [4., np.nan, 6.],
                      [7., 8., 9.]], dtype='float64')

        y = np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 0.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros,
                         )

        a = x
        a[1, 1] = 0

        # for whatever reason, numpy's fft has very limited precision, and
        # the comparison fails unless you cast the float64 to a float16
        if hasattr(np, 'float16'):
            assert_array_almost_equal_nulp(np.asarray(z, dtype=np.float16),
                                           np.asarray(a, dtype=np.float16), 10)
        assert np.all(np.abs(z - a) < 1e-14)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3x3_withnan(self, boundary, interpolate_nan,
                                 normalize_kernel, ignore_edge_zeros):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        x = np.array([[0., 0., 3.],
                      [1., np.nan, 0.],
                      [0., 2., 0.]], dtype='float64')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                         interpolate_nan=interpolate_nan,
                         normalize_kernel=normalize_kernel,
                         ignore_edge_zeros=ignore_edge_zeros,
                         )

        w_n = np.array([[3., 5., 3.],
                        [5., 8., 5.],
                        [3., 5., 3.]], dtype='float64')
        w_z = np.array([[4., 6., 4.],
                        [6., 9., 6.],
                        [4., 6., 4.]], dtype='float64')
        answer_dict = {
            'sum': np.array([[1., 4., 3.],
                             [3., 6., 5.],
                             [3., 3., 2.]], dtype='float64'),
            'sum_wrap': np.array([[6., 6., 6.],
                                  [6., 6., 6.],
                                  [6., 6., 6.]], dtype='float64'),
        }
        answer_dict['average'] = answer_dict['sum'] / w_z
        answer_dict['average_ignan'] = answer_dict['sum'] / w_n
        answer_dict['average_wrap_ignan'] = answer_dict['sum_wrap'] / 8.
        answer_dict['average_wrap'] = answer_dict['sum_wrap'] / 9.
        answer_dict['average_withzeros'] = answer_dict['sum'] / 9.
        answer_dict['average_withzeros_ignan'] = answer_dict['sum'] / 8.
        answer_dict['sum_withzeros'] = answer_dict['sum']
        answer_dict['sum_ignan'] = answer_dict['sum']
        answer_dict['sum_withzeros_ignan'] = answer_dict['sum']
        answer_dict['sum_wrap_ignan'] = answer_dict['sum_wrap']

        if normalize_kernel:
            answer_key = 'average'
        else:
            answer_key = 'sum'

        if boundary == 'wrap':
            answer_key += '_wrap'
        elif not ignore_edge_zeros:
            answer_key += '_withzeros'

        if interpolate_nan:
            answer_key += '_ignan'

        a = answer_dict[answer_key]
        # for reasons unknown, the Windows FFT returns an answer for the [0, 0]
        # component that is EXACTLY 10*np.spacing
        assert np.all(np.abs(z - a) <= np.spacing(np.where(z > a, z, a)) * 10)

    def test_big_fail(self):
        """ Test that convolve_fft raises an exception if a too-large array is passed in """

        with pytest.raises((ValueError, MemoryError)):
            # while a good idea, this approach did not work; it actually writes to disk
            #arr = np.memmap('file.np', mode='w+', shape=(512, 512, 512), dtype=np.complex)
            # this just allocates the memory but never touches it; it's better:
            arr = np.empty([512, 512, 512], dtype=np.complex)
            # note 512**3 * 16 bytes = 2.0 GB
            convolve_fft(arr, arr)
