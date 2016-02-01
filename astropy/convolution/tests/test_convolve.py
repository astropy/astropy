# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ...tests.helper import pytest

from ..convolve import convolve

from numpy.testing import assert_array_almost_equal_nulp

VALID_DTYPES = []
for dtype_array in ['>f4', '<f4', '>f8', '<f8']:
    for dtype_kernel in ['>f4', '<f4', '>f8', '<f8']:
        VALID_DTYPES.append((dtype_array, dtype_kernel))

BOUNDARY_OPTIONS = [None, 'fill', 'wrap', 'extend']


class TestConvolve1D(object):
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """

        x = [1, 4, 5, 6, 5, 7, 8]
        y = [0.2, 0.6, 0.2]
        z = convolve(x, y, boundary=None)
        assert_array_almost_equal_nulp(z,
            np.array([0., 3.6, 5., 5.6, 5.6, 6.8, 0.]), 10)

    @pytest.mark.parametrize(('dtype_array', 'dtype_kernel'), VALID_DTYPES)
    def test_dtype(self, dtype_array, dtype_kernel):
        '''
        Test that 32- and 64-bit floats are correctly handled
        '''

        x = np.array([1., 2., 3.], dtype=dtype_array)

        y = np.array([0., 1., 0.], dtype=dtype_kernel)

        z = convolve(x, y)

        assert x.dtype == z.dtype

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_1_none(self, boundary):
        '''
        Test that a unit kernel with a single element returns the same array
        '''

        x = np.array([1., 2., 3.], dtype='>f8')

        y = np.array([1.], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        assert np.all(z == x)

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

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(z == np.array([0., 4., 0.], dtype='>f8'))
        elif boundary == 'fill':
            assert np.all(z == np.array([1., 4., 3.], dtype='>f8'))
        elif boundary == 'wrap':
            assert np.all(z == np.array([4., 4., 4.], dtype='>f8'))
        else:
            assert np.all(z == np.array([2., 4., 6.], dtype='>f8'))

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_3_withnan(self, boundary):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='>f8')

        y = np.array([0., 1., 0.], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        assert np.isnan(z[1])
        x = np.nan_to_num(z)
        z = np.nan_to_num(z)

        if boundary is None:
            assert np.all(z == np.array([[0., 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]], dtype='>f8'))
        else:
            assert np.all(z == x)

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_uniform_3_withnan(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements. This version includes a NaN
        value in the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='>f8')

        y = np.array([1., 1., 1.], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([0., 6., 0.], dtype='>f8'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([3., 6., 5.], dtype='>f8'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([6., 6., 6.], dtype='>f8'), 10)
        else:
            assert_array_almost_equal_nulp(z, np.array([4., 6., 8.], dtype='>f8'), 10)


class TestConvolve2D(object):
    def test_list(self):
        """
        Test that convolve works correctly when inputs are lists
        """
        x = [[1, 1, 1],
             [1, 1, 1],
             [1, 1, 1]]

        z = convolve(x, x, boundary='fill', fill_value=1)
        assert_array_almost_equal_nulp(z / 9, x, 10)

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

        z = convolve(x, y, boundary=boundary)

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

        z = convolve(x, y, boundary=boundary)

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
    def test_uniform_3x3_withnan(self, boundary):
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

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert_array_almost_equal_nulp(z, np.array([[0., 0., 0.],
                                                       [0., 9., 0.],
                                                       [0., 0., 0.]], dtype='>f8'), 10)
        elif boundary == 'fill':
            assert_array_almost_equal_nulp(z, np.array([[2., 6., 5.],
                                                        [5., 9., 8.],
                                                        [5., 5., 4.]], dtype='>f8'), 10)
        elif boundary == 'wrap':
            assert_array_almost_equal_nulp(z, np.array([[9., 9., 9.],
                                                        [9., 9., 9.],
                                                        [9., 9., 9.]], dtype='>f8'), 10)
        else:
            assert_array_almost_equal_nulp(z, np.array([[3., 10., 17.],
                                                        [6., 9., 12.],
                                                        [9., 8., 7.]], dtype='>f8'), 10)


class TestConvolve3D(object):
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

        z = convolve(x, x, boundary='fill', fill_value=1)
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

        z = convolve(x, y, boundary=boundary)

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

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_unity_3x3x3_withnan(self, boundary):
        '''
        Test that a 3x3 unit kernel returns the same array (except when
        boundary is None). This version includes a NaN value in the original
        array.
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., np.nan, 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.zeros((3, 3, 3), dtype='>f8')
        y[1, 1, 1] = 1.

        z = convolve(x, y, boundary=boundary)

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
    def test_uniform_3x3x3_withnan(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        x = np.array([[[1., 2., 1.], [2., 3., 1.], [3., 2., 5.]],
                      [[4., 3., 1.], [5., np.nan, 2.], [6., 1., 1.]],
                      [[7., 0., 2.], [8., 2., 3.], [9., 2., 2.]]], dtype='>f8')

        y = np.ones((3, 3, 3), dtype='>f8')

        z = convolve(x, y, boundary=boundary)

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
