# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np
from numpy.testing import assert_allclose
from ..sampling import downsample, upsample

DATA1 = np.array([[0., 1.], [2., 3.]])
DATA2 = np.array([[0., 0., 0.25, 0.25],
                  [0., 0., 0.25, 0.25],
                  [0.5, 0.5, 0.75, 0.75],
                  [0.5, 0.5, 0.75, 0.75]])
DATA3 = np.array([[0., 0., 0., 1., 1., 1.],
                  [0., 0., 0., 1., 1., 1.],
                  [0., 0., 0., 1., 1., 1.],
                  [2., 2., 2., 3., 3., 3.],
                  [2., 2., 2., 3., 3., 3.],
                  [2., 2., 2., 3., 3., 3.]])
DATA4 = np.array([[0., 0., 0.25, 0.25, 1.],
                  [0., 0., 0.25, 0.25, 2.],
                  [0.5, 0.5, 0.75, 0.75, 3.],
                  [0.5, 0.5, 0.75, 0.75, 4.],
                  [5.5, 5.5, 5.75, 5.75, 5.]])


class TestUpsample(object):
    """Test upsampling of images."""
    def test_factor1(self):
        """Test upsampling of an image by a factor of 1."""
        img = upsample(DATA1, 1)
        assert_allclose(img, DATA1, rtol=0, atol=1e-6)

    def test_factor2(self):
        """Test upsampling of an image by a factor of 2."""
        img = upsample(DATA1, 2)
        assert_allclose(img, DATA2, rtol=0, atol=1e-6)

    def test_factor3(self):
        """Test upsampling of an image by a factor of 3."""
        img = upsample(DATA1*9., 3)
        assert_allclose(img, DATA3, rtol=0, atol=1e-6)


class TestDownsample(object):
    """Test downsampling of images."""
    def test_factor1(self):
        """Test downsampling of an image by a factor of 1."""
        img = downsample(DATA1, 1)
        assert_allclose(img, DATA1, rtol=0, atol=1e-6)

    def test_factor2_2x2(self):
        """Test downsampling of a 2x2 image by a factor of 2."""
        img = downsample(DATA1, 2)
        assert_allclose(img, np.array([[6.]]), rtol=0, atol=1e-6)

    def test_factor2(self):
        """Test downsampling of an image by a factor of 2."""
        img = downsample(DATA2, 2)
        assert_allclose(img, DATA1, rtol=0, atol=1e-6)

    def test_factor3(self):
        """Test downsampling of an image by a factor of 3."""
        img = downsample(DATA3/9., 3)
        assert_allclose(img, DATA1, rtol=0, atol=1e-6)

    def test_factor2_roundtrip(self):
        """
        Test roundtrip of upsampling then downsampling of an image by a
        factor of 2.
        """
        img = downsample(upsample(DATA1, 2), 2)
        assert_allclose(img, DATA1, rtol=0, atol=1e-6)

    def test_factor3_roundtrip(self):
        """
        Test roundtrip of upsampling then downsampling of an image by a
        factor of 3.
        """
        img = downsample(upsample(DATA1, 3), 3)
        assert_allclose(img, DATA1, rtol=0, atol=1e-6)

    def test_factor3_4x4size(self):
        """
        Test downsampling when image dimensions are not a whole multiple
        of factor.  A 4x4 image downsampled by a factor of 3.
        """
        img = downsample(DATA2, 3)
        assert_allclose(img, np.array([[2.25]]), rtol=0, atol=1e-6)

    def test_factor2_5x5size(self):
        """
        Test downsampling when image dimensions are not a whole multiple
        of factor.  A 5x5 image downsampled by a factor of 2.
        """
        img = downsample(DATA4, 2)
        assert_allclose(img, DATA1, rtol=0, atol=1e-6)

    def test_factor4_6x6size(self):
        """
        Test downsampling when image dimensions are not a whole multiple
        of factor.  A 6x6 image downsampled by a factor for 4.
        """
        img = downsample(DATA3, 4)
        assert_allclose(img, np.array([[12.]]), rtol=0, atol=1e-6)
