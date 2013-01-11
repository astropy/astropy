# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from ....tests.helper import pytest
from ....config import ConfigurationItem

from ..convolve import convolve_fft

from numpy.testing import assert_array_almost_equal_nulp

import itertools

try:
    import scipy.fftpack
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import fftw3
    HAS_FFTW = True
except ImportError:
    HAS_FFTW = False

VALID_DTYPES = []
for dtype_array in ['>f4', '<f4', '>f8', '<f8']:
    for dtype_kernel in ['>f4', '<f4', '>f8', '<f8']:
        VALID_DTYPES.append((dtype_array, dtype_kernel))

BOUNDARY_OPTIONS = [None, 'fill', 'wrap']

"""
What does convolution mean?  We use the 'same size' assumption here (i.e.,
you expect an array of the exact same size as the one you put in)
Convolving any array with a kernel that is [1] should result in the same array returned
Working example array: [1,2,3,4,5]
Convolved with [1] = [1,2,3,4,5]
Convolved with [1,1] = [1, 3, 5, 7, 9] THIS IS NOT CONSISTENT!
Convolved with [1,0] = [1, 2, 3, 4, 5]
Convolved with [0,1] = [0, 1, 2, 3, 4]
"""

# NOTE: use_numpy_fft is redundant if you don't have FFTW installed
option_names = ('boundary', 'interpolate_nan', 'normalize_kernel', 'ignore_edge_zeros', 'fft_type')
options = list(itertools.product(BOUNDARY_OPTIONS, (True, False), (True, False), (True, False), ['numpy', 'scipy', 'fftw']))


class TestConvolve1D(object):

    @pytest.mark.parametrize(option_names, options)
    def test_unity_1_none(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that a unit kernel with a single element returns the same array
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([1., 2., 3.], dtype='float64')

        y = np.array([1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        assert_array_almost_equal_nulp(z, x, 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None).
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([1., 2., 3.], dtype='float64')

        y = np.array([0., 1., 0.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        assert_array_almost_equal_nulp(z, x, 10)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([1., 0., 3.], dtype='float64')

        y = np.array([1., 1., 1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        answer_dict = {
                'sum': np.array([1., 4., 3.], dtype='float64'),
                'sum_zeros': np.array([1., 4., 3.], dtype='float64'),
                'sum_nozeros': np.array([1., 4., 3.], dtype='float64'),
                'average': np.array([1 / 3., 4 / 3., 1.], dtype='float64'),
                'sum_wrap':  np.array([4., 4., 4.], dtype='float64'),
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

        print answer_key
        assert_array_almost_equal_nulp(z, answer_dict[answer_key], 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3_withnan(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([1., np.nan, 3.], dtype='float64')

        y = np.array([0., 1., 0.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        if fft_type == 'numpy':
            # for whatever reason, numpy's fft has very limited precision, and
            # the comparison fails unless you cast the float64 to a float16
            # REMOVED because of numpy 1.4 incompatibility assert_array_almost_equal_nulp(np.asarray(z, dtype=np.float16), np.array([1.,0.,3.], dtype=np.float16), 10)
            # ASSERT equality to better than 16 bit but worse than 32 bit precision
            assert np.all(np.abs(z - np.array([1., 0., 3.])) < 1e-14)
        else:
            assert_array_almost_equal_nulp(z, np.array([1., 0., 3.], dtype='float64'), 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_1_withnan(self, boundary, interpolate_nan, normalize_kernel,
            ignore_edge_zeros, fft_type):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([1., np.nan, 3.], dtype='float64')

        y = np.array([1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        if fft_type == 'numpy':
            # for whatever reason, numpy's fft has very limited precision, and
            # the comparison fails unless you cast the float64 to a float16
            # np1.4 incompatible assert_array_almost_equal_nulp(np.asarray(z, dtype=np.float16), np.array([1.,0.,3.], dtype=np.float16), 10)
            assert np.all(np.abs(z - np.array([1., 0., 3.])) < 1e-14)
        else:
            assert_array_almost_equal_nulp(z, np.array([1., 0., 3.], dtype='float64'), 10)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3_withnan(self, boundary, interpolate_nan,
            normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements. This version includes a NaN
        value in the original array.
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([1., np.nan, 3.], dtype='float64')

        y = np.array([1., 1., 1.], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        answer_dict = {
                'sum': np.array([1., 4., 3.], dtype='float64'),
                'sum_nozeros': np.array([1., 4., 3.], dtype='float64'),
                'sum_zeros': np.array([1., 4., 3.], dtype='float64'),
                'sum_zeros_noignan': np.array([1., 4., 3.], dtype='float64'),
                'sum_nozeros_noignan': np.array([1., 4., 3.], dtype='float64'),
                'average': np.array([1., 2., 3.], dtype='float64'),
                'sum_wrap':  np.array([4., 4., 4.], dtype='float64'),
                'sum_wrap_noignan':  np.array([4., 4., 4.], dtype='float64'),
                'average_wrap': np.array([(1 + 3) / 2., 2., 2.], dtype='float64'),
                'average_wrap_noignan': np.array([4 / 3., 4 / 3., 4 / 3.], dtype='float64'),
                'average_nozeros': np.array([1, 2, 3], dtype='float64'),
                'average_nozeros_noignan': np.array([1 / 2., 4 / 3., 3 / 2.], dtype='float64'),
                'average_zeros': np.array([1 / 2., 4 / 2., 3 /2.], dtype='float64'),
                'average_zeros_noignan': np.array([1/3., 4/3., 3/3.], dtype='float64'),
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

        print boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, answer_key
        assert_array_almost_equal_nulp(z, answer_dict[answer_key], 10)


class TestConvolve2D(object):

    @pytest.mark.parametrize(option_names, options)
    def test_unity_1x1_none(self, boundary, interpolate_nan, normalize_kernel,
            ignore_edge_zeros, fft_type):
        '''
        Test that a 1x1 unit kernel returns the same array
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype='float64')

        y = np.array([[1.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        assert_array_almost_equal_nulp(z, x, 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3x3(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None).
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype='float64')

        y = np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 0.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

        assert_array_almost_equal_nulp(z, x, 10)
        #assert np.all( np.abs(z-x) < np.spacing(np.where(z>x,z,x))*2 )

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3x3(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel.
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

        x = np.array([[0., 0., 3.],
                      [1., 0., 0.],
                      [0., 2., 0.]], dtype='float64')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='float64')

        z = convolve_fft(x, y, boundary=boundary,
                interpolate_nan=interpolate_nan,
                normalize_kernel=normalize_kernel,
                ignore_edge_zeros=ignore_edge_zeros,
                fft_type=fft_type)

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
        print boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, answer_key
        print "z: ", z
        print "answer: ", a
        print "ratio: ", z / a
        # for reasons unknown, the Windows FFT returns an answer for the [0,0]
        # component that is EXACTLY 10*np.spacing
        assert np.all(np.abs(z - a) <= np.spacing(np.where(z > a, z, a)) * 10)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3x3_withnan(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None). This version includes a NaN value in the original
        array.
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

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
                fft_type=fft_type)

        a = x
        a[1, 1] = 0

        if fft_type == 'numpy':
            # for whatever reason, numpy's fft has very limited precision, and
            # the comparison fails unless you cast the float64 to a float16
            # np1.4 incompatible assert_array_almost_equal_nulp(np.asarray(z, dtype=np.float16), np.asarray(a, dtype=np.float16), 10)
            assert np.all(np.abs(z - a) < 1e-14)
        else:
            assert_array_almost_equal_nulp(z, a, 10)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3x3_withnan(self, boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, fft_type):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        if fft_type == 'fftw' and not HAS_FFTW:
            pytest.skip('fftw3 is not installed')
        elif fft_type == 'scipy' and not HAS_SCIPY:
            pytest.skip('scipy is not installed')

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
                fft_type=fft_type)

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
        print boundary, interpolate_nan, normalize_kernel, ignore_edge_zeros, answer_key
        print "z: ", z
        print "answer: ", a
        print "ratio: ", z / a
        # for reasons unknown, the Windows FFT returns an answer for the [0,0]
        # component that is EXACTLY 10*np.spacing
        assert np.all(np.abs(z - a) <= np.spacing(np.where(z > a, z, a)) * 10)
