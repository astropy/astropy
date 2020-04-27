# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy.nddata import block_reduce, block_replicate


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
