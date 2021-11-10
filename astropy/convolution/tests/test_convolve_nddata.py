# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import numpy as np

from astropy.convolution.convolve import convolve, convolve_fft
from astropy.convolution.kernels import Gaussian2DKernel
from astropy.nddata import NDData


def test_basic_nddata():
    arr = np.zeros((11, 11))
    arr[5, 5] = 1
    ndd = NDData(arr)
    test_kernel = Gaussian2DKernel(1)

    result = convolve(ndd, test_kernel)

    x, y = np.mgrid[:11, :11]
    expected = result[5, 5] * np.exp(-0.5 * ((x - 5)**2 + (y - 5)**2))

    np.testing.assert_allclose(result, expected, atol=1e-6)

    resultf = convolve_fft(ndd, test_kernel)
    np.testing.assert_allclose(resultf, expected, atol=1e-6)


@pytest.mark.parametrize('convfunc',
   [lambda *args: convolve(*args, nan_treatment='interpolate', normalize_kernel=True),
    lambda *args: convolve_fft(*args, nan_treatment='interpolate', normalize_kernel=True)])
def test_masked_nddata(convfunc):
    arr = np.zeros((11, 11))
    arr[4, 5] = arr[6, 5] = arr[5, 4] = arr[5, 6] = 0.2
    arr[5, 5] = 1.5
    ndd_base = NDData(arr)

    mask = arr < 0  # this is all False
    mask[5, 5] = True
    ndd_mask = NDData(arr, mask=mask)

    arrnan = arr.copy()
    arrnan[5, 5] = np.nan
    ndd_nan = NDData(arrnan)

    test_kernel = Gaussian2DKernel(1)

    result_base = convfunc(ndd_base, test_kernel)
    result_nan = convfunc(ndd_nan, test_kernel)
    result_mask = convfunc(ndd_mask, test_kernel)

    assert np.allclose(result_nan, result_mask)
    assert not np.allclose(result_base, result_mask)
    assert not np.allclose(result_base, result_nan)

    # check to make sure the mask run doesn't talk back to the initial array
    assert np.sum(np.isnan(ndd_base.data)) != np.sum(np.isnan(ndd_nan.data))
