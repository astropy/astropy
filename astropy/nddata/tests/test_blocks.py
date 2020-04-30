# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy.nddata import reshape_as_blocks, block_reduce, block_replicate


class TestReshapeAsBlocks:
    def test_1d(self):
        data = np.arange(16)
        reshaped = reshape_as_blocks(data, 2)
        assert reshaped.shape == (8, 2)
        reshaped = reshape_as_blocks(data, 4)
        assert reshaped.shape == (4, 4)
        reshaped = reshape_as_blocks(data, 8)
        assert reshaped.shape == (2, 8)

    def test_2d(self):
        data = np.arange(16).reshape(4, 4)
        reshaped = reshape_as_blocks(data, (2, 2))
        assert reshaped.shape == (2, 2, 2, 2)

        data = np.arange(64).reshape(8, 8)
        reshaped = reshape_as_blocks(data, (2, 2))
        assert reshaped.shape == (4, 4, 2, 2)
        reshaped = reshape_as_blocks(data, (4, 4))
        assert reshaped.shape == (2, 2, 4, 4)

    def test_3d(self):
        data = np.arange(64).reshape(4, 4, 4)
        reshaped = reshape_as_blocks(data, (2, 2, 2))
        assert reshaped.shape == (2, 2, 2, 2, 2, 2)

        data = np.arange(2*3*4).reshape(2, 3, 4)
        reshaped = reshape_as_blocks(data, (2, 1, 2))
        assert reshaped.shape == (1, 3, 2, 2, 1, 2)

    def test_view(self):
        data = np.arange(16).reshape(4, 4)
        reshaped = reshape_as_blocks(data, (2, 2))
        data[0, 0] = 100
        assert reshaped[0, 0, 0, 0] == 100

    def test_invalid_block_dim(self):
        data = np.arange(64).reshape(4, 4, 4)
        match = ('block_size must be a scalar or have the same '
                 'length as the number of data dimensions')
        with pytest.raises(ValueError, match=match):
            reshape_as_blocks(data, (2, 2))

    def test_invalid_block_size(self):
        data = np.arange(16).reshape(4, 4)
        match = ('Each dimension of block_size must divide evenly '
                 'into the corresponding dimension of data')
        with pytest.raises(ValueError, match=match):
            reshape_as_blocks(data, (2, 3))

    def test_invalid_block_value(self):
        data = np.arange(16).reshape(4, 4)
        match = 'block_size elements must be integers'
        with pytest.raises(ValueError, match=match):
            reshape_as_blocks(data, (2.1, 2))

        match = 'block_size elements must be strictly positive'
        with pytest.raises(ValueError, match=match):
            reshape_as_blocks(data, (-1, 0))


class TestBlockReduce:
    def test_1d(self):
        """Test 1D array."""
        data = np.arange(4)
        expected = np.array([1, 5])
        result = block_reduce(data, 2)
        assert np.all(result == expected)

    def test_1d_mean(self):
        """Test 1D array with func=np.mean."""
        data = np.arange(4)
        block_size = 2.
        expected = block_reduce(data, block_size, func=np.sum) / block_size
        result_mean = block_reduce(data, block_size, func=np.mean)
        assert np.all(result_mean == expected)

    def test_2d(self):
        """Test 2D array."""
        data = np.arange(4).reshape(2, 2)
        expected = np.array([[6]])
        result = block_reduce(data, 2)
        assert np.all(result == expected)

    def test_2d_mean(self):
        """Test 2D array with func=np.mean."""
        data = np.arange(4).reshape(2, 2)
        block_size = 2.
        expected = (block_reduce(data, block_size, func=np.sum) /
                    block_size**2)
        result = block_reduce(data, block_size, func=np.mean)
        assert np.all(result == expected)

    def test_2d_trim(self):
        """
        Test trimming of 2D array when size is not perfectly divisible
        by block_size.
        """

        data1 = np.arange(15).reshape(5, 3)
        result1 = block_reduce(data1, 2)
        data2 = data1[0:4, 0:2]
        result2 = block_reduce(data2, 2)
        assert np.all(result1 == result2)

    def test_block_size_broadcasting(self):
        """Test scalar block_size broadcasting."""
        data = np.arange(16).reshape(4, 4)
        result1 = block_reduce(data, 2)
        result2 = block_reduce(data, (2, 2))
        assert np.all(result1 == result2)

    def test_block_size_len(self):
        """Test block_size length."""
        data = np.ones((2, 2))
        with pytest.raises(ValueError):
            block_reduce(data, (2, 2, 2))


class TestBlockReplicate:
    def test_1d(self):
        """Test 1D array."""
        data = np.arange(2)
        expected = np.array([0, 0, 0.5, 0.5])
        result = block_replicate(data, 2)
        assert np.all(result == expected)

    def test_1d_conserve_sum(self):
        """Test 1D array with conserve_sum=False."""
        data = np.arange(2)
        block_size = 2.
        expected = block_replicate(data, block_size) * block_size
        result = block_replicate(data, block_size, conserve_sum=False)
        assert np.all(result == expected)

    def test_2d(self):
        """Test 2D array."""
        data = np.arange(2).reshape(2, 1)
        expected = np.array([[0, 0], [0, 0], [0.25, 0.25], [0.25, 0.25]])
        result = block_replicate(data, 2)
        assert np.all(result == expected)

    def test_2d_conserve_sum(self):
        """Test 2D array with conserve_sum=False."""
        data = np.arange(6).reshape(2, 3)
        block_size = 2.
        expected = block_replicate(data, block_size) * block_size**2
        result = block_replicate(data, block_size, conserve_sum=False)
        assert np.all(result == expected)

    def test_block_size_broadcasting(self):
        """Test scalar block_size broadcasting."""
        data = np.arange(4).reshape(2, 2)
        result1 = block_replicate(data, 2)
        result2 = block_replicate(data, (2, 2))
        assert np.all(result1 == result2)

    def test_block_size_len(self):
        """Test block_size length."""
        data = np.arange(5)
        with pytest.raises(ValueError):
            block_replicate(data, (2, 2))
