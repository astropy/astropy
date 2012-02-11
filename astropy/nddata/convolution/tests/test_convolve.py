import numpy as np

from ....tests.helper import pytest

from ..convolve import convolve


VALID_DTYPES = []
for dtype_array in ['>f4', '<f4', '>f8', '<f8']:
    for dtype_kernel in ['>f4', '<f4', '>f8', '<f8']:
        VALID_DTYPES.append((dtype_array, dtype_kernel))

BOUNDARY_OPTIONS = [None, 'fill', 'wrap', 'extend']


class TestConvolve2D(object):

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

        x = np.array([[0., 0., 27.],
                      [9., 0., 0.],
                      [0., 18., 0.]], dtype='>f8')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(z == np.array([[0., 0., 0.],
                                         [0., 0., 0.],
                                         [0., 0., 0.]], dtype='>f8'))
        elif boundary == 'fill':
            assert np.all(z == np.array([[1., 4., 3.],
                                         [3., 6., 5.],
                                         [3., 3., 2.]], dtype='>f8'))
        elif boundary == 'wrap':
            assert np.all(z == np.array([[6., 6., 6.],
                                         [6., 6., 6.],
                                         [6., 6., 6.]], dtype='>f8'))
        else:
            assert np.all(z == np.array([[2., 7., 12.],
                                         [4., 6., 8.],
                                         [6., 5., 4.]], dtype='>f8'))

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
    def test_uniform_3x3(self, boundary):
        '''
        Test that the different modes are producing the correct results using
        a 3x3 uniform kernel. This version includes a NaN value in the
        original array.
        '''

        x = np.array([[0., 0., 27.],
                      [9., np.nan, 0.],
                      [0., 18., 0.]], dtype='>f8')

        y = np.array([[1., 1., 1.],
                      [1., 1., 1.],
                      [1., 1., 1.]], dtype='>f8')

        z = convolve(x, y, boundary=boundary)

        if boundary is None:
            assert np.all(z == np.array([[0.00, 0.00, 0.00],
                                         [0.00, 0.75, 0.00],
                                         [0.00, 0.00, 0.00]], dtype='>f8'))
        elif boundary == 'fill':
            assert np.all(z == np.array([[1.75, 4.75, 3.75],
                                         [3.75, 6.75, 5.75],
                                         [3.75, 3.75, 2.75]], dtype='>f8'))
        elif boundary == 'wrap':
            assert np.all(z == np.array([[6.75, 6.75, 6.75],
                                         [6.75, 6.75, 6.75],
                                         [6.75, 6.75, 6.75]], dtype='>f8'))
        else:
            assert np.all(z == np.array([[2.75, 7.75, 12.75],
                                         [4.75, 6.75, 8.75],
                                         [6.75, 5.75, 4.75]], dtype='>f8'))
