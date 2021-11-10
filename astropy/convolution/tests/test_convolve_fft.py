# Licensed under a 3-clause BSD style license - see LICENSE.rst

import itertools
from contextlib import nullcontext

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_array_almost_equal_nulp, assert_array_equal

from astropy import units as u
from astropy.convolution.convolve import convolve, convolve_fft
from astropy.utils.exceptions import AstropyUserWarning

VALID_DTYPES = ('>f4', '<f4', '>f8', '<f8')
VALID_DTYPE_MATRIX = list(itertools.product(VALID_DTYPES, VALID_DTYPES))

BOUNDARY_OPTIONS = [None, 'fill', 'wrap']
NANTREATMENT_OPTIONS = ('interpolate', 'fill')
NORMALIZE_OPTIONS = [True, False]
PRESERVE_NAN_OPTIONS = [True, False]

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
option_names = ('boundary', 'nan_treatment', 'normalize_kernel', 'dealias')
options = list(itertools.product(BOUNDARY_OPTIONS,
                                 NANTREATMENT_OPTIONS,
                                 (True, False),
                                 (True, False)
                                 ))
option_names_preserve_nan = ('boundary', 'nan_treatment',
                             'normalize_kernel', 'preserve_nan')
options_preserve_nan = list(itertools.product(BOUNDARY_OPTIONS,
                                              NANTREATMENT_OPTIONS,
                                              (True, False),
                                              (True, False)))


def expected_boundary_warning(boundary=None):
    # Helper that returns the appropriate context manager for the boundary=None
    # warning depending on the value of boundary.
    if boundary is None:
        ctx = pytest.warns(AstropyUserWarning,
                           match='The convolve_fft version of boundary=None '
                                 'is equivalent to the convolve boundary=\'fill\'')
    else:
        ctx = nullcontext()
    return ctx


def expected_dealias_error(boundary=None, dealias=False):
    # Helper that returns the appropriate context manager for the boundary=None
    # warning depending on the value of boundary.

    if dealias and boundary == 'wrap':
        ctx = pytest.raises(ValueError)
    else:
        ctx = nullcontext()
    return ctx


def assert_floatclose(x, y):
    """Assert arrays are close to within expected floating point rounding.

    Check that the result is correct at the precision expected for 64 bit
    numbers, taking account that the tolerance has to reflect that all powers
    in the FFTs enter our values.
    """
    # The number used is set by the fact that the Windows FFT sometimes
    # returns an answer that is EXACTLY 10*np.spacing.
    assert_allclose(x, y, atol=10*np.spacing(x.max()), rtol=0.)


class TestConvolve1D:

    @pytest.mark.parametrize(option_names, options)
    def test_quantity(self, boundary, nan_treatment, normalize_kernel, dealias):
        """
        Test that convolve_fft works correctly when input array is a Quantity
        """

        x = np.array([1., 4., 5., 6., 5., 7., 8.], dtype='float64') * u.ph
        y = np.array([0.2, 0.6, 0.2], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

                assert x.unit == z.unit

    @pytest.mark.parametrize(option_names, options)
    def test_unity_1_none(self, boundary, nan_treatment, normalize_kernel, dealias):
        '''
        Test that a unit kernel with a single element returns the same array
        '''

        x = np.array([1., 2., 3.], dtype='float64')

        y = np.array([1.], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

                assert_floatclose(z, x)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3(self, boundary, nan_treatment, normalize_kernel, dealias):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None).
        '''

        x = np.array([1., 2., 3.], dtype='float64')

        y = np.array([0., 1., 0.], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

                assert_floatclose(z, x)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3(self, boundary, nan_treatment, normalize_kernel, dealias):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements
        '''

        x = np.array([1., 0., 3.], dtype='float64')

        y = np.array([1., 1., 1.], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

                answer_key = (boundary, nan_treatment, normalize_kernel)

                answer_dict = {
                    'sum_fill_zeros': np.array([1., 4., 3.], dtype='float64'),
                    'average_fill_zeros': np.array([1 / 3., 4 / 3., 1.], dtype='float64'),
                    'sum_wrap': np.array([4., 4., 4.], dtype='float64'),
                    'average_wrap': np.array([4 / 3., 4 / 3., 4 / 3.], dtype='float64'),
                }

                result_dict = {
                    # boundary, nan_treatment, normalize_kernel
                    ('fill', 'interpolate', True): answer_dict['average_fill_zeros'],
                    ('wrap', 'interpolate', True): answer_dict['average_wrap'],
                    ('fill', 'interpolate', False): answer_dict['sum_fill_zeros'],
                    ('wrap', 'interpolate', False): answer_dict['sum_wrap'],
                }
                for k in list(result_dict.keys()):
                    result_dict[(k[0], 'fill', k[2])] = result_dict[k]
                for k in list(result_dict.keys()):
                    if k[0] == 'fill':
                        result_dict[(None, k[1], k[2])] = result_dict[k]

                assert_floatclose(z, result_dict[answer_key])

    @pytest.mark.parametrize(option_names, options)
    def test_halfity_3(self, boundary, nan_treatment, normalize_kernel, dealias):
        '''
        Test that the different modes are producing the correct results using
        a uniform, non-unity kernel with three elements
        '''

        x = np.array([1., 0., 3.], dtype='float64')

        y = np.array([0.5, 0.5, 0.5], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

                answer_dict = {
                    'sum': np.array([0.5, 2.0, 1.5], dtype='float64'),
                    'sum_zeros': np.array([0.5, 2., 1.5], dtype='float64'),
                    'sum_nozeros': np.array([0.5, 2., 1.5], dtype='float64'),
                    'average': np.array([1 / 3., 4 / 3., 1.], dtype='float64'),
                    'sum_wrap': np.array([2., 2., 2.], dtype='float64'),
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
                else:
                    # average = average_zeros; sum = sum_zeros
                    answer_key += '_zeros'

                assert_floatclose(z, answer_dict[answer_key])

    @pytest.mark.parametrize(option_names_preserve_nan, options_preserve_nan)
    def test_unity_3_withnan(self, boundary, nan_treatment, normalize_kernel,
                             preserve_nan):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='float64')

        y = np.array([0., 1., 0.], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            z = convolve_fft(x, y, boundary=boundary,
                             nan_treatment=nan_treatment,
                             normalize_kernel=normalize_kernel,
                             preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[1])

        z = np.nan_to_num(z)

        assert_floatclose(z, [1., 0., 3.])

    inputs = (np.array([1., np.nan, 3.], dtype='float64'),
              np.array([1., np.inf, 3.], dtype='float64'))
    outputs = (np.array([1., 0., 3.], dtype='float64'),
               np.array([1., 0., 3.], dtype='float64'))
    options_unity1withnan = list(itertools.product(BOUNDARY_OPTIONS,
                                                   NANTREATMENT_OPTIONS,
                                                   (True, False),
                                                   (True, False),
                                                   inputs, outputs))

    @pytest.mark.parametrize(option_names_preserve_nan + ('inval', 'outval'),
                             options_unity1withnan)
    def test_unity_1_withnan(self, boundary, nan_treatment, normalize_kernel,
                             preserve_nan, inval, outval):
        '''
        Test that a unit kernel with three elements returns the same array
        (except when boundary is None). This version includes a NaN value in
        the original array.
        '''

        x = inval

        y = np.array([1.], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            z = convolve_fft(x, y, boundary=boundary,
                             nan_treatment=nan_treatment,
                             normalize_kernel=normalize_kernel,
                             preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[1])

        z = np.nan_to_num(z)

        assert_floatclose(z, outval)

    @pytest.mark.parametrize(option_names_preserve_nan, options_preserve_nan)
    def test_uniform_3_withnan(self, boundary, nan_treatment,
                               normalize_kernel, preserve_nan):
        '''
        Test that the different modes are producing the correct results using
        a uniform kernel with three elements. This version includes a NaN
        value in the original array.
        '''

        x = np.array([1., np.nan, 3.], dtype='float64')

        y = np.array([1., 1., 1.], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            z = convolve_fft(x, y, boundary=boundary,
                             nan_treatment=nan_treatment,
                             normalize_kernel=normalize_kernel,
                             preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[1])

        answer_dict = {
            'sum': np.array([1., 4., 3.], dtype='float64'),
            'sum_nozeros': np.array([1., 4., 3.], dtype='float64'),
            'sum_zeros': np.array([1., 4., 3.], dtype='float64'),
            'sum_nozeros_interpnan': np.array([1., 4., 3.], dtype='float64'),
            'average': np.array([1., 2., 3.], dtype='float64'),
            'sum_wrap': np.array([4., 4., 4.], dtype='float64'),
            'average_wrap': np.array([4/3., 4/3., 4/3.], dtype='float64'),
            'average_wrap_interpnan': np.array([2, 2, 2], dtype='float64'),
            'average_nozeros': np.array([1/2., 4/3., 3/2.], dtype='float64'),
            'average_nozeros_interpnan': np.array([1., 2., 3.], dtype='float64'),
            'average_zeros': np.array([1 / 3., 4 / 3., 3 / 3.], dtype='float64'),
            'average_zeros_interpnan': np.array([1 / 2., 4 / 2., 3 / 2.], dtype='float64'),
        }

        for key in list(answer_dict.keys()):
            if 'sum' in key:
                answer_dict[key+"_interpnan"] = answer_dict[key] * 3./2.

        if normalize_kernel:
            answer_key = 'average'
        else:
            answer_key = 'sum'

        if boundary == 'wrap':
            answer_key += '_wrap'
        else:
            # average = average_zeros; sum = sum_zeros
            answer_key += '_zeros'

        if nan_treatment == 'interpolate':
            answer_key += '_interpnan'

        posns = np.isfinite(z)

        answer = answer_dict[answer_key][posns]

        # check that fill is set and that the 1'th position that was originally
        # NaN is included in the check
        if (nan_treatment == 'fill') and posns[1]:
            # we fill the center with the sum of the input array divided by
            # three, since we've now pre-filled the center value with zero
            answer[1] = 4 / (3. if normalize_kernel else 1.)

        assert_floatclose(z[posns], answer)

    def test_nan_interpolate(self):

        # Test masked array
        array = np.array([1., np.nan, 3.], dtype='float64')
        kernel = np.array([1, 1, 1])
        masked_array = np.ma.masked_array(array, mask=[0, 1, 0])

        result = convolve_fft(masked_array, kernel, boundary='fill',
                              nan_treatment='interpolate',
                              fill_value=np.nan)

        assert_floatclose(result, [1, 2, 3])

    def test_nan_fill(self):
        # regression for #8121

        # Test masked array
        array = np.array([1., np.nan, 3.], dtype='float64')
        kernel = np.array([1, 1, 1])

        result = convolve_fft(array, kernel, boundary='fill',
                              nan_treatment='fill',
                              fill_value=0)

        # note that, because fill_value also affects boundary='fill', the edge
        # pixels are treated as zero rather than being ignored.
        assert_floatclose(result, [1/3., 4/3., 1.])

    def test_nan_fill_two(self):
        # regression for #8121

        # Test masked array
        array = np.array([1., np.nan, 3.], dtype='float64')
        kernel = np.array([1, 1, 1])

        result = convolve_fft(array, kernel, boundary='fill',
                              nan_treatment='fill',
                              fill_value=1)

        # note that, because fill_value also affects boundary='fill', the edge
        # pixels are treated as fill_value=1 rather than being ignored.
        assert_floatclose(result, [1., 5/3., 5/3.])

    def test_masked_array(self):
        """
        Check whether convolve_fft works with masked arrays.
        """

        # Test masked array
        array = np.array([1., 2., 3.], dtype='float64')
        kernel = np.array([1, 1, 1])
        masked_array = np.ma.masked_array(array, mask=[0, 1, 0])
        result = convolve_fft(masked_array, kernel, boundary='fill',
                              fill_value=0.)
        assert_floatclose(result, [1./2, 2, 3./2])

        # Now test against convolve()
        convolve_result = convolve(masked_array, kernel, boundary='fill',
                              fill_value=0.)
        assert_floatclose(convolve_result, result)

        # Test masked kernel
        array = np.array([1., 2., 3.], dtype='float64')
        kernel = np.array([1, 1, 1])
        masked_kernel = np.ma.masked_array(kernel, mask=[0, 1, 0])
        result = convolve_fft(array, masked_kernel, boundary='fill',
                              fill_value=0.)
        assert_floatclose(result, [1, 2, 1])

        # Now test against convolve()
        convolve_result = convolve(array, masked_kernel, boundary='fill',
                              fill_value=0.)
        assert_floatclose(convolve_result, result)

    def test_normalize_function(self):
        """
        Check if convolve_fft works when passing a normalize function.
        """
        array = [1, 2, 3]
        kernel = [3, 3, 3]
        result = convolve_fft(array, kernel, normalize_kernel=np.max)
        assert_floatclose(result, [3, 6, 5])

    @pytest.mark.parametrize(option_names, options)
    def test_normalization_is_respected(self, boundary,
                                        nan_treatment,
                                        normalize_kernel,
                                        dealias):
        """
        Check that if normalize_kernel is False then the normalization
        tolerance is respected.
        """
        array = np.array([1, 2, 3])
        # A simple identity kernel to which a non-zero normalization is added.
        base_kernel = np.array([1.0])

        # Use the same normalization error tolerance in all cases.
        normalization_rtol = 1e-4

        # Add the error below to the kernel.
        norm_error = [normalization_rtol / 10, normalization_rtol * 10]

        for err in norm_error:
            kernel = base_kernel + err
            result = convolve_fft(array, kernel,
                                  normalize_kernel=normalize_kernel,
                                  nan_treatment=nan_treatment,
                                  normalization_zero_tol=normalization_rtol)
            if normalize_kernel:
                # Kernel has been normalized to 1.
                assert_floatclose(result, array)
            else:
                # Kernel should not have been normalized...
                assert_floatclose(result, array * kernel)


class TestConvolve2D:

    @pytest.mark.parametrize(option_names, options)
    def test_unity_1x1_none(self, boundary, nan_treatment, normalize_kernel, dealias):
        '''
        Test that a 1x1 unit kernel returns the same array
        '''

        x = np.array([[1., 2., 3.],
                      [4., 5., 6.],
                      [7., 8., 9.]], dtype='float64')

        y = np.array([[1.]], dtype='float64')

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

                assert_floatclose(z, x)

    @pytest.mark.parametrize(option_names, options)
    def test_unity_3x3(self, boundary, nan_treatment, normalize_kernel, dealias):
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

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

                assert_floatclose(z, x)

    @pytest.mark.parametrize(option_names, options)
    def test_uniform_3x3(self, boundary, nan_treatment, normalize_kernel, dealias):
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

        with expected_boundary_warning(boundary=boundary):
            with expected_dealias_error(boundary=boundary, dealias=dealias):
                z = convolve_fft(x, y, boundary=boundary,
                                 nan_treatment=nan_treatment,
                                 fill_value=np.nan if normalize_kernel else 0,
                                 normalize_kernel=normalize_kernel,
                                 dealias=dealias)

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
                elif nan_treatment == 'fill':
                    answer_key += '_withzeros'

                a = answer_dict[answer_key]
                assert_floatclose(z, a)

    @pytest.mark.parametrize(option_names_preserve_nan, options_preserve_nan)
    def test_unity_3x3_withnan(self, boundary, nan_treatment,
                               normalize_kernel, preserve_nan):
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

        with expected_boundary_warning(boundary=boundary):
            z = convolve_fft(x, y, boundary=boundary,
                             nan_treatment=nan_treatment,
                             normalize_kernel=normalize_kernel,
                             preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[1, 1])
            z = np.nan_to_num(z)

        x = np.nan_to_num(x)

        assert_floatclose(z, x)

    @pytest.mark.parametrize(option_names_preserve_nan, options_preserve_nan)
    def test_uniform_3x3_withnan(self, boundary, nan_treatment,
                                 normalize_kernel, preserve_nan):
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

        # commented out: allow unnormalized nan-ignoring convolution
        # # kernel is not normalized, so this situation -> exception
        # if nan_treatment and not normalize_kernel:
        #     with pytest.raises(ValueError):
        #         z = convolve_fft(x, y, boundary=boundary,
        #                          nan_treatment=nan_treatment,
        #                          normalize_kernel=normalize_kernel,
        #                          ignore_edge_zeros=ignore_edge_zeros,
        #                          )
        #     return

        with expected_boundary_warning(boundary=boundary):
            z = convolve_fft(x, y, boundary=boundary,
                             nan_treatment=nan_treatment,
                             # you cannot fill w/nan, you can only interpolate over it
                             fill_value=np.nan if normalize_kernel and nan_treatment=='interpolate' else 0,
                             normalize_kernel=normalize_kernel,
                             preserve_nan=preserve_nan)

        if preserve_nan:
            assert np.isnan(z[1, 1])

        # weights
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
        answer_dict['average_interpnan'] = answer_dict['sum'] / w_n
        answer_dict['average_wrap_interpnan'] = answer_dict['sum_wrap'] / 8.
        answer_dict['average_wrap'] = answer_dict['sum_wrap'] / 9.
        answer_dict['average_withzeros'] = answer_dict['sum'] / 9.
        answer_dict['average_withzeros_interpnan'] = answer_dict['sum'] / 8.
        answer_dict['sum_withzeros'] = answer_dict['sum']
        answer_dict['sum_interpnan'] = answer_dict['sum'] * 9/8.
        answer_dict['sum_withzeros_interpnan'] = answer_dict['sum']
        answer_dict['sum_wrap_interpnan'] = answer_dict['sum_wrap'] * 9/8.

        if normalize_kernel:
            answer_key = 'average'
        else:
            answer_key = 'sum'

        if boundary == 'wrap':
            answer_key += '_wrap'
        elif nan_treatment == 'fill':
            answer_key += '_withzeros'

        if nan_treatment == 'interpolate':
            answer_key += '_interpnan'

        answer_dict[answer_key]

        # Skip the NaN at [1, 1] when preserve_nan=True
        posns = np.where(np.isfinite(z))

        # for reasons unknown, the Windows FFT returns an answer for the [0, 0]
        # component that is EXACTLY 10*np.spacing
        assert_floatclose(z[posns], z[posns])

    def test_big_fail(self):
        """ Test that convolve_fft raises an exception if a too-large array is passed in."""

        with pytest.raises((ValueError, MemoryError)):
            # while a good idea, this approach did not work; it actually writes to disk
            # arr = np.memmap('file.np', mode='w+', shape=(512, 512, 512), dtype=complex)
            # this just allocates the memory but never touches it; it's better:
            arr = np.empty([512, 512, 512], dtype=complex)
            # note 512**3 * 16 bytes = 2.0 GB
            convolve_fft(arr, arr)

    def test_padding(self):
        """
        Test that convolve_fft pads to _next_fast_lengths and does not expand all dimensions
        to length of longest side (#11242/#10047).
        """

        # old implementation expanded this to up to 2048**3
        shape = (1, 1226, 518)
        img = np.zeros(shape, dtype='float64')
        img[0, 600:610, 300:304] = 1.0
        kernel = np.zeros((1, 7, 7), dtype='float64')
        kernel[0, 3, 3] = 1.0

        with pytest.warns(AstropyUserWarning,
                          match="psf_pad was set to False, which overrides the boundary='fill'"):
            img_fft = convolve_fft(img, kernel, return_fft=True, psf_pad=False, fft_pad=False)
            assert_array_equal(img_fft.shape, shape)
            img_fft = convolve_fft(img, kernel, return_fft=True, psf_pad=False, fft_pad=True)
            # should be from either hardcoded _good_sizes[] or scipy.fft.next_fast_len()
            assert img_fft.shape in ((1, 1250, 540), (1, 1232, 525))

        img_fft = convolve_fft(img, kernel, return_fft=True, psf_pad=True, fft_pad=False)
        assert_array_equal(img_fft.shape, np.array(shape) + np.array(kernel.shape))
        img_fft = convolve_fft(img, kernel, return_fft=True, psf_pad=True, fft_pad=True)
        assert img_fft.shape in ((2, 1250, 540), (2, 1250, 525))

    @pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
    def test_non_normalized_kernel(self, boundary):

        x = np.array([[0., 0., 4.],
                      [1., 2., 0.],
                      [0., 3., 0.]], dtype='float')

        y = np.array([[1., -1., 1.],
                      [-1., 0., -1.],
                      [1., -1., 1.]], dtype='float')

        with expected_boundary_warning(boundary=boundary):
            z = convolve_fft(x, y, boundary=boundary, nan_treatment='fill',
                             normalize_kernel=False)

        if boundary in (None, 'fill'):
            assert_floatclose(z, np.array([[1., -5., 2.],
                                           [1., 0., -3.],
                                           [-2., -1., -1.]], dtype='float'))
        elif boundary == 'wrap':
            assert_floatclose(z, np.array([[0., -8., 6.],
                                           [5., 0., -4.],
                                           [2., 3., -4.]], dtype='float'))
        else:
            raise ValueError("Invalid boundary specification")


@pytest.mark.parametrize(('boundary'), BOUNDARY_OPTIONS)
def test_asymmetric_kernel(boundary):
    '''
    Make sure that asymmetric convolution
    functions go the right direction
    '''

    x = np.array([3., 0., 1.], dtype='>f8')

    y = np.array([1, 2, 3], dtype='>f8')

    with expected_boundary_warning(boundary=boundary):
        z = convolve_fft(x, y, boundary=boundary, normalize_kernel=False)

    if boundary in (None, 'fill'):
        assert_array_almost_equal_nulp(z, np.array([6., 10., 2.], dtype='float'), 10)
    elif boundary == 'wrap':
        assert_array_almost_equal_nulp(z, np.array([9., 10., 5.], dtype='float'), 10)


@pytest.mark.parametrize(('boundary', 'nan_treatment',
                          'normalize_kernel', 'preserve_nan', 'dtype'),
                         itertools.product(BOUNDARY_OPTIONS,
                                           NANTREATMENT_OPTIONS,
                                           NORMALIZE_OPTIONS,
                                           PRESERVE_NAN_OPTIONS,
                                           VALID_DTYPES))
def test_input_unmodified(boundary, nan_treatment,
                          normalize_kernel, preserve_nan, dtype):
    """
    Test that convolve_fft works correctly when inputs are lists
    """

    array = [1., 4., 5., 6., 5., 7., 8.]
    kernel = [0.2, 0.6, 0.2]
    x = np.array(array, dtype=dtype)
    y = np.array(kernel, dtype=dtype)

    # Make pseudoimmutable
    x.flags.writeable = False
    y.flags.writeable = False

    with expected_boundary_warning(boundary=boundary):
        z = convolve_fft(x, y, boundary=boundary, nan_treatment=nan_treatment,
                        normalize_kernel=normalize_kernel, preserve_nan=preserve_nan)

    assert np.all(np.array(array, dtype=dtype) == x)
    assert np.all(np.array(kernel, dtype=dtype) == y)


@pytest.mark.parametrize(('boundary', 'nan_treatment',
                          'normalize_kernel', 'preserve_nan', 'dtype'),
                         itertools.product(BOUNDARY_OPTIONS,
                                           NANTREATMENT_OPTIONS,
                                           NORMALIZE_OPTIONS,
                                           PRESERVE_NAN_OPTIONS,
                                           VALID_DTYPES))
def test_input_unmodified_with_nan(boundary, nan_treatment,
                                   normalize_kernel, preserve_nan, dtype):
    """
    Test that convolve_fft doesn't modify the input data
    """

    array = [1., 4., 5., np.nan, 5., 7., 8.]
    kernel = [0.2, 0.6, 0.2]
    x = np.array(array, dtype=dtype)
    y = np.array(kernel, dtype=dtype)

    # Make pseudoimmutable
    x.flags.writeable = False
    y.flags.writeable = False

    # make copies for post call comparison
    x_copy = x.copy()
    y_copy = y.copy()

    with expected_boundary_warning(boundary=boundary):
        z = convolve_fft(x, y, boundary=boundary, nan_treatment=nan_treatment,
                        normalize_kernel=normalize_kernel, preserve_nan=preserve_nan)

    # ( NaN == NaN ) = False
    # Only compare non NaN values for canonical equivalence
    # and then check NaN explicitly with np.isnan()
    array_is_nan = np.isnan(array)
    kernel_is_nan = np.isnan(kernel)
    array_not_nan = ~array_is_nan
    kernel_not_nan = ~kernel_is_nan
    assert np.all(x_copy[array_not_nan] == x[array_not_nan])
    assert np.all(y_copy[kernel_not_nan] == y[kernel_not_nan])
    assert np.all(np.isnan(x[array_is_nan]))
    assert np.all(np.isnan(y[kernel_is_nan]))


@pytest.mark.parametrize('error_kwarg', [{'psf_pad': True}, {'fft_pad': True}, {'dealias': True}])
def test_convolve_fft_boundary_wrap_error(error_kwarg):
    x = np.array([[1., 2., 3.],
                  [4., 5., 6.],
                  [7., 8., 9.]], dtype='>f8')
    y = np.array([[1.]], dtype='>f8')
    assert (convolve_fft(x, y, boundary='wrap') == x).all()

    with pytest.raises(ValueError) as err:
        convolve_fft(x, y, boundary='wrap', **error_kwarg)
    assert str(err.value) == \
        f"With boundary='wrap', {list(error_kwarg.keys())[0]} cannot be enabled."


def test_convolve_fft_boundary_extend_error():
    x = np.array([[1., 2., 3.],
                  [4., 5., 6.],
                  [7., 8., 9.]], dtype='>f8')
    y = np.array([[1.]], dtype='>f8')

    with pytest.raises(NotImplementedError) as err:
        convolve_fft(x, y, boundary='extend')
    assert str(err.value) == \
        "The 'extend' option is not implemented for fft-based convolution"
