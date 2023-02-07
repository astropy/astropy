import pytest

import numpy as np
from numpy.testing import assert_equal

from astropy.utils.slicing import simplify_basic_index

from hypothesis import given
from hypothesis.extra.numpy import basic_indices

TEST_SHAPE = (13, 16, 4, 90)


class TestSanitizeBasicIndex:

    @pytest.fixture(autouse=True)
    def setup_class(self, tmp_path):
        self.shape = TEST_SHAPE
        self.data = np.random.random(TEST_SHAPE)

    @given(basic_indices(TEST_SHAPE))
    def test_indexing(self, index):
        print('before', index)
        new_index = simplify_basic_index(index, shape=self.shape)
        print('after', index)
        print('new_index', new_index)
        assert_equal(self.data[index], self.data[new_index])
        assert isinstance(new_index, tuple)
        assert len(new_index) == len(self.shape)
        for idim, idx in enumerate(new_index):
            if isinstance(idx, slice):
                assert isinstance(idx.start, int)
                assert idx.start >= 0
                assert idx.start < TEST_SHAPE[idim]
                if idx.step is None or idx.step > 0:
                    assert isinstance(idx.stop, int)
                    assert idx.stop >= 0
                    assert idx.stop <= TEST_SHAPE[idim]
                assert isinstance(idx.step, int)
            elif isinstance(idx, int):
                assert idx >= 0
            else:
                raise TypeError('Expected slice or int')
