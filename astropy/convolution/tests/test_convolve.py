# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from ..convolve import convolve, convolve_fft

from numpy.testing import assert_array_almost_equal_nulp, assert_array_almost_equal

import itertools

VALID_DTYPES = []
for dtype_array in ['>f4', '<f4', '>f8', '<f8']:
    for dtype_kernel in ['>f4', '<f4', '>f8', '<f8']:
        VALID_DTYPES.append((dtype_array, dtype_kernel))

BOUNDARY_OPTIONS = [None, 'fill', 'wrap', 'extend']
NANHANDLING_OPTIONS = ['interpolate', 'fill']
NORMALIZE_OPTIONS = [True, False]
PRESERVE_NAN_OPTIONS = [True, False]

BOUNDARIES_AND_CONVOLUTIONS = (list(zip(itertools.cycle((convolve,)),
                                        BOUNDARY_OPTIONS)) + [(convolve_fft,
                                                               'wrap'),
                                                              (convolve_fft,
                                                               'fill')])
HAS_SCIPY = True
try:
    import scipy
except ImportError:
    HAS_SCIPY = False


class TestConvolve1D:
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """

        x = [1, 4, 5, 6, 5, 7, 8]
        y = [0.2, 0.6, 0.2]
        z = convolve(x, y, boundary=None)
        assert_array_almost_equal_nulp(z,
            np.array([0., 3.6, 5., 5.6, 5.6, 6.8, 0.]), 10)

    def test_input_unmodified(self):
        """
        Test that convolve works correctly when inputs are lists
        """

        inlist = [1, 4, 5, 6, 5, 7, 8]
        x = np.array(inlist)
        y = [0.2, 0.6, 0.2]
        z = convolve(x, y, boundary=None)

        assert np.all(np.array(inlist) == x)

    @pytest.mark.parametrize(('dtype_array', 'dtype_kernel'), VALID_DTYPES)
    def test_dtype(self, dtype_array, dtype_kernel):
        '''
        Test that 32- and 64-bit floats are correctly handled
        '''

        x = np.array([1., 2., 3.], dtype=dtype_array)

        y = np.array([0., 1., 0.], dtype=dtype_kernel)

        z = convolve(x, y)

        assert x.dtype == z.dtype

    @pytest.mark.parametrize(('convfunc', 'boundary',), BOUNDARIES_AND_CONVOLUTIONS)
    def test_unity_1_none(self, boundary, convfunc):
        '''
        Test that a unit kernel with a single element returns the same array
        '''

        x = np.array([1., 2., 3.], dtype='>f8')

        y = np.array([1.], dtype='>f8')

        z = convfunc(x, y, boundary=boundary)

        np.testing.assert_allclose(z, x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_3(self, boundary):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None).
        '''

        x = np.array([1., 2., 3.], dtype='>f8')

        y = np.array([0., 1., 0.], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(z == np.array([0., 2., 0.], dtype='>f8'))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements
        '''

        x = np.array([1., 0., 3.], dtype='>f8')

        y = np.array([1., 1., 1.], dtype='>f8')

        z = convolve(x, y, boundary=boundary, normalize_kernel=False)

        if boundary is None:
            assert np.all(z == np.array([0., 4., 0.], dtype='>f8'))
        elif boundary == 'fill':
            assert np.all(z == np.array([1., 4., 3.], dtype='>f8'))
        elif boundary == 'wrap':
            assert np.all(z == np.array([4., 4., 4.], dtype='>f8'))
        else:
            assert np.all(z == np.array([2., 4., 6.], dtype='>f8'))

    @pytest.mark.parametrize(('boundary', 'nan_treatment',
                              'normalize_kernel', 'preserve_nan'),
                             itertools.product(BOUNDARY_OPTIONS,
                                               NANHANDLING_OPTIONS,
                                               NORMALIZE_OPTIONS,
                                               PRESERVE_NAN_OPTIONS))
    def test_unity_3_withnan(self, boundary, nan_treatment,
                             normalize_kernel, preserve_nan):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='>f8')

        y = np.array([0., 1., 0.], dtype='>f8')

        z = convolve(x, y, boundary=boundary, nan_treatment=nan_treatment,
                     normalize_kernel=normalize_kernel,
                     preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[1])

        x = np.nan_to_num(z)
        z = np.nan_to_num(z)

        if boundary is None:
            assert np.all(z == np.array([0., 0., 0.], dtype='>f8'))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(('boundary', 'nan_treatment',
                              'normalize_kernel', 'preserve_nan'),
                             itertools.product(BOUNDARY_OPTIONS,
                                               NANHANDLING_OPTIONS,
                                               NORMALIZE_OPTIONS,
                                               PRESERVE_NAN_OPTIONS))
    def test_uniform_3_withnan(self, boundary, nan_treatment, normalize_kernel,
                               preserve_nan):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements. This version includes a NaN
        value in the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='>f8')

        y = np.array([1., 1., 1.], dtype='>f8')

        z = convolve(x, y, boundary=boundary, nan_treatment=nan_treatment,
                     normalize_kernel=normalize_kernel,
                     preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[1])

        z = np.nan_to_num(z)

        # boundary, nan_treatment, normalize_kernel
        rslt = {
                (None, 'interpolate', True): [0, 2, 0],
                (None, 'interpolate', False): [0, 6, 0],
                (None, 'fill', True): [0, 4/3., 0],
                (None, 'fill', False): [0, 4, 0],
                ('fill', 'interpolate', True): [1/2., 2, 3/2.],
                ('fill', 'interpolate', False): [3/2., 6, 9/2.],
                ('fill', 'fill', True): [1/3., 4/3., 3/3.],
                ('fill', 'fill', False): [1, 4, 3],
                ('wrap', 'interpolate', True): [2, 2, 2],
                ('wrap', 'interpolate', False): [6, 6, 6],
                ('wrap', 'fill', True): [4/3., 4/3., 4/3.],
                ('wrap', 'fill', False): [4, 4, 4],
                ('extend', 'interpolate', True): [1, 2, 3],
                ('extend', 'interpolate', False): [3, 6, 9],
                ('extend', 'fill', True): [2/3., 4/3., 6/3.],
                ('extend', 'fill', False): [2, 4, 6],
               }[boundary, nan_treatment, normalize_kernel]
        if preserve_nan:
            rslt[1] = 0

        assert_array_almost_equal_nulp(z, np.array(rslt, dtype='>f8'), 10)


class TestConvolve2D:
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """
        x = [[1, 1, 1],
             [1, 1, 1],
             [1, 1, 1]]

        z = convolve(x, x, boundary='fill', fill_value=1, normalize_kernel=True)
        assert_array_almost_equal_nulp(z, x, 10)
        z = convolve(x, x, boundary='fill', fill_value=1, normalize_kernel=False)
        assert_array_almost_equal_nulp(z, np.array(x, float)*9, 10)

    @pytest.mark.parametrize(('dtype_array', 'dtype_kernel'), VALID_DTYPES)
    def test_dtype(self, dtype_array, dtype_kernel):
        '''
        Test that 32- and 64-bit floats are correctly handled
        '''

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype=dtype_array)

        y = np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 0.]], dtype=dtype_kernel)

        z = convolve(x, y)

        assert x.dtype == z.dtype

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_1x1_none(self, boundary):
        '''
        Test that a 1x1 unit kernel returns the same array
        '''

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype='>f8')

        y = np.array([[1.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_3x3(self, boundary):
        '''
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None).
        '''

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype='>f8')

        y = np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 0.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(z == np.array([[0., 0., 0.],
                                         [0., 5., 0.],
                                         [0., 0., 0.]], dtype='>f8'))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3x3(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel.
        '''

        x = np.array([[0., 0., 3.],
                      [1., 0., 0.],
                      [0., 2., 0.]], dtype='>f8')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary, normalize_kernel=False)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[0., 0., 0.],
                                                        [0., 6., 0.],
                                                        [0., 0., 0.]], dtype='>f8'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[1., 4., 3.],
                                                        [3., 6., 5.],
                                                        [3., 3., 2.]], dtype='>f8'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([[6., 6., 6.],
                                                        [6., 6., 6.],
                                                        [6., 6., 6.]], dtype='>f8'), 10)
        else:
            assert_array_almost_equal_nulp(z, np.array([[2., 7., 12.],
                                                        [4., 6., 8.],
                                                        [6., 5., 4.]], dtype='>f8'), 10)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_3x3_withnan(self, boundary):
        '''
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None). This version includes a NaN value in the original
        array.
        '''

        x = np.array([[1., 2., 3.],
                      [4., np.nan, 6.],
                      [7., 8., 9.]], dtype='>f8')

        y = np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 0.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary, nan_treatment='fill',
                     preserve_nan=True)

        assert np.isnan(z[1, 1])
        x = np.nan_to_num(z)
        z = np.nan_to_num(z)

        if boundary is None:
            assert np.all(z == np.array([[0., 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]], dtype='>f8'))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3x3_withnanfilled(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        x = np.array([[0., 0., 4.],
                      [1., np.nan, 0.],
                      [0., 3., 0.]], dtype='>f8')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary, nan_treatment='fill',
                     normalize_kernel=False)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[0., 0., 0.],
                                                        [0., 8., 0.],
                                                        [0., 0., 0.]], dtype='>f8'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[1., 5., 4.],
                                                        [4., 8., 7.],
                                                        [4., 4., 3.]], dtype='>f8'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([[8., 8., 8.],
                                                        [8., 8., 8.],
                                                        [8., 8., 8.]], dtype='>f8'), 10)
        elif boundary == 'extend':
            assert_array_almost_equal_nulp(z, np.array([[2., 9., 16.],
                                                        [5., 8., 11.],
                                                        [8., 7., 6.]], dtype='>f8'), 10)
        else:
            raise ValueError("Invalid boundary specification")

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3x3_withnaninterped(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        x = np.array([[0., 0., 4.],
                      [1., np.nan, 0.],
                      [0., 3., 0.]], dtype='>f8')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary, nan_treatment='interpolate',
                     normalize_kernel=True)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[0., 0., 0.],
                                                        [0., 1., 0.],
                                                        [0., 0., 0.]], dtype='>f8'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[1./8, 5./8, 4./8],
                                                        [4./8, 8./8, 7./8],
                                                        [4./8, 4./8, 3./8]], dtype='>f8'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([[1., 1., 1.],
                                                        [1., 1., 1.],
                                                        [1., 1., 1.]], dtype='>f8'), 10)
        elif boundary == 'extend':
            assert_array_almost_equal_nulp(z, np.array([[2./8, 9./8, 16./8],
                                                        [5./8, 8./8, 11./8],
                                                        [8./8, 7./8, 6./8]], dtype='>f8'), 10)
        else:
            raise ValueError("Invalid boundary specification")

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_non_normalized_kernel_2D(self, boundary):

        x = np.array([[0., 0., 4.],
                      [1., 2., 0.],
                      [0., 3., 0.]], dtype='float')

        y = np.array([[1., -1., 1.],
                      [-1., 0., -1.],
                      [1., -1., 1.]], dtype='float')

        z = convolve(x, y, boundary=boundary, nan_treatment='fill',
                     normalize_kernel=False)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[0., 0., 0.],
                                                        [0., 0., 0.],
                                                        [0., 0., 0.]], dtype='float'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[1., -5., 2.],
                                                        [1., 0., -3.],
                                                        [-2., -1., -1.]], dtype='float'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([[0., -8., 6.],
                                                        [5., 0., -4.],
                                                        [2., 3., -4.]], dtype='float'), 10)
        elif boundary == 'extend':
            assert_array_almost_equal_nulp(z, np.array([[2., -1., -2.],
                                                        [0., 0., 1.],
                                                        [2., -4., 2.]], dtype='float'), 10)
        else:
            raise ValueError("Invalid boundary specification")


class TestConvolve3D:
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """
        x = [[[1, 1, 1],
              [1, 1, 1],
              [1, 1, 1]],
             [[1, 1, 1],
              [1, 1, 1],
              [1, 1, 1]],
             [[1, 1, 1],
              [1, 1, 1],
              [1, 1, 1]]]

        z = convolve(x, x, boundary='fill', fill_value=1, normalize_kernel=False)
        assert_array_almost_equal_nulp(z / 27, x, 10)

    @pytest.mark.parametrize(('dtype_array', 'dtype_kernel'), VALID_DTYPES)
    def test_dtype(self, dtype_array, dtype_kernel):
        '''
        Test that 32- and 64-bit floats are correctly handled
        '''

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype=dtype_array)

        y = np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 0.]], dtype=dtype_kernel)

        z = convolve(x, y)

        assert x.dtype == z.dtype

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_1x1x1_none(self, boundary):
        '''
        Test that a 1x1x1 unit kernel returns the same array
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., 0., 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.array([[[1.]]], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_3x3x3(self, boundary):
        '''
        Test that a 3x3x3 unit kernel returns the same array (except when
        boundary is None).
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., 3., 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.zeros((3, 3, 3), dtype='>f8')
        y[1, 1, 1] = 1.

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(z == np.array([[[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                                         [[0., 0., 0.], [0., 3., 0.], [0., 0., 0.]],
                                         [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]], dtype='>f8'))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3x3x3(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel.
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., 3., 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.ones((3, 3, 3), dtype='>f8')

        z = convolve(x, y, boundary=boundary, normalize_kernel=False)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                                                       [[0., 0., 0.], [0., 81., 0.], [0., 0., 0.]],
                                                       [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]], dtype='>f8'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[[23., 28., 16.], [35., 46., 25.], [25., 34., 18.]],
                                                       [[40., 50., 23.], [63., 81., 36.], [46., 60., 27.]],
                                                       [[32., 40., 16.], [50., 61., 22.], [36., 44., 16.]]], dtype='>f8'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([[[81., 81., 81.], [81., 81., 81.], [81., 81., 81.]],
                                                       [[81., 81., 81.], [81., 81., 81.], [81., 81., 81.]],
                                                       [[81., 81., 81.], [81., 81., 81.], [81., 81., 81.]]], dtype='>f8'), 10)
        else:
            assert_array_almost_equal_nulp(z, np.array([[[65., 54., 43.], [75., 66., 57.], [85., 78., 71.]],
                                                       [[96., 71., 46.], [108., 81., 54.], [120., 91., 62.]],
                                                       [[127., 88., 49.], [141., 96., 51.], [155., 104., 53.]]], dtype='>f8'), 10)

    @pytest.mark.parametrize(('boundary', 'nan_treatment'),
                             itertools.product(BOUNDARY_OPTIONS,
                                               NANHANDLING_OPTIONS))
    def test_unity_3x3x3_withnan(self, boundary, nan_treatment):
        '''
        Test that a 3x3x3 unit kernel returns the same array (except when
        boundary is None). This version includes a NaN value in the original
        array.
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., np.nan, 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.zeros((3, 3, 3), dtype='>f8')
        y[1, 1, 1] = 1.

        z = convolve(x, y, boundary=boundary, nan_treatment=nan_treatment,
                     preserve_nan=True)

        assert np.isnan(z[1, 1, 1])
        x = np.nan_to_num(z)
        z = np.nan_to_num(z)

        if boundary is None:
            assert np.all(z == np.array([[[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                                         [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                                         [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]], dtype='>f8'))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3x3x3_withnan_filled(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., np.nan, 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.ones((3, 3, 3), dtype='>f8')

        z = convolve(x, y, boundary=boundary, nan_treatment='fill',
                     normalize_kernel=False)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                                                        [[0., 0., 0.], [0., 78., 0.], [0., 0., 0.]],
                                                        [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]], dtype='>f8'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[[20., 25., 13.],
                                                         [32., 43., 22.],
                                                         [22., 31., 15.]],
                                                        [[37., 47., 20.],
                                                         [60., 78., 33.],
                                                         [43., 57., 24.]],
                                                        [[29., 37., 13.],
                                                         [47., 58., 19.],
                                                         [33., 41., 13.]]], dtype='>f8'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([[[78., 78., 78.], [78., 78., 78.], [78., 78., 78.]],
                                                        [[78., 78., 78.], [78., 78., 78.], [78., 78., 78.]],
                                                        [[78., 78., 78.], [78., 78., 78.], [78., 78., 78.]]], dtype='>f8'), 10)
        elif boundary == 'extend':
            assert_array_almost_equal_nulp(z, np.array([[[62., 51., 40.],
                                                         [72., 63., 54.],
                                                         [82., 75., 68.]],
                                                        [[93., 68., 43.],
                                                         [105., 78., 51.],
                                                         [117., 88., 59.]],
                                                        [[124., 85., 46.],
                                                         [138., 93., 48.],
                                                         [152., 101., 50.]]],
                                                       dtype='>f8'), 10)
        else:
            raise ValueError("Invalid Boundary Option")

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3x3x3_withnan_interped(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., np.nan, 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.ones((3, 3, 3), dtype='>f8')

        z = convolve(x, y, boundary=boundary, nan_treatment='interpolate',
                     normalize_kernel=True)

        kernsum = y.sum() - 1  # one nan is missing
        mid = x[np.isfinite(x)].sum() / kernsum

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]],
                                                        [[0., 0., 0.], [0., 78., 0.], [0., 0., 0.]],
                                                        [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]],
                                                       dtype='>f8')/kernsum, 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[[20., 25., 13.],
                                                         [32., 43., 22.],
                                                         [22., 31., 15.]],
                                                        [[37., 47., 20.],
                                                         [60., 78., 33.],
                                                         [43., 57., 24.]],
                                                        [[29., 37., 13.],
                                                         [47., 58., 19.],
                                                         [33., 41., 13.]]],
                                                       dtype='>f8')/kernsum, 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.tile(mid.astype('>f8'), [3, 3, 3]), 10)
        elif boundary == 'extend':
            assert_array_almost_equal_nulp(z, np.array([[[62., 51., 40.],
                                                         [72., 63., 54.],
                                                         [82., 75., 68.]],
                                                        [[93., 68., 43.],
                                                         [105., 78., 51.],
                                                         [117., 88., 59.]],
                                                        [[124., 85., 46.],
                                                         [138., 93., 48.],
                                                         [152., 101., 50.]]],
                                                       dtype='>f8')/kernsum, 10)
        else:
            raise ValueError("Invalid Boundary Option")


@pytest.mark.parametrize(('convfunc', 'boundary'), BOUNDARIES_AND_CONVOLUTIONS)
def test_asymmetric_kernel(boundary, convfunc):
    '''
    Regression test for #6264: make sure that asymmetric convolution
    functions go the right direction
    '''

    x = np.array([3., 0., 1.], dtype='>f8')

    y = np.array([1, 2, 3], dtype='>f8')

    z = convolve(x, y, boundary=boundary, normalize_kernel=False)

    if boundary == 'fill':
        assert_array_almost_equal_nulp(z, np.array([6., 10., 2.], dtype='float'), 10)
    elif boundary is None:
        assert_array_almost_equal_nulp(z, np.array([0., 10., 0.], dtype='float'), 10)
    elif boundary == 'extend':
        assert_array_almost_equal_nulp(z, np.array([15., 10., 3.], dtype='float'), 10)
    elif boundary == 'wrap':
        assert_array_almost_equal_nulp(z, np.array([9., 10., 5.], dtype='float'), 10)


@pytest.mark.parametrize('ndims', (1, 2, 3))
def test_convolution_consistency(ndims):

    np.random.seed(0)
    array = np.random.randn(*([3]*ndims))
    np.random.seed(0)
    kernel = np.random.rand(*([3]*ndims))

    conv_f = convolve_fft(array, kernel, boundary='fill')
    conv_d = convolve(array, kernel, boundary='fill')

    assert_array_almost_equal_nulp(conv_f, conv_d, 30)


def test_astropy_convolution_against_numpy():
    x = np.array([1, 2, 3])
    y = np.array([5, 4, 3, 2, 1])

    assert_array_almost_equal(np.convolve(y, x, 'same'),
                              convolve(y, x, normalize_kernel=False))
    assert_array_almost_equal(np.convolve(y, x, 'same'),
                              convolve_fft(y, x, normalize_kernel=False))


@pytest.mark.skipif('not HAS_SCIPY')
def test_astropy_convolution_against_scipy():
    from scipy.signal import fftconvolve
    x = np.array([1, 2, 3])
    y = np.array([5, 4, 3, 2, 1])

    assert_array_almost_equal(fftconvolve(y, x, 'same'),
                              convolve(y, x, normalize_kernel=False))
    assert_array_almost_equal(fftconvolve(y, x, 'same'),
                              convolve_fft(y, x, normalize_kernel=False))
