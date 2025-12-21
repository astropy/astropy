"""ASV-style benchmark for separable-kernel convolve performance.

This file is intended to be copied into astropy-benchmarks/benchmarks.
"""

from __future__ import annotations

import contextlib
import importlib

import numpy as np

from astropy.convolution import convolve

IMAGE_SIZES = [128, 256, 512]
KERNEL_SIZES = [11, 21, 31]
BOUNDARIES = ["fill"]


def _make_gaussian_1d(size: int, sigma: float | None = None) -> np.ndarray:
    if size % 2 == 0:
        raise ValueError("Kernel size must be odd.")
    if sigma is None:
        sigma = size / 6.0
    x = np.arange(size, dtype=float) - size // 2
    kernel = np.exp(-(x * x) / (2.0 * sigma * sigma))
    kernel /= kernel.sum()
    return kernel


@contextlib.contextmanager
def _disable_separable_fast_path(convolve_mod):
    original = convolve_mod._factor_separable_kernel_2d
    convolve_mod._factor_separable_kernel_2d = lambda _kernel: None
    try:
        yield
    finally:
        convolve_mod._factor_separable_kernel_2d = original


class ConvolveSeparable2D:
    params = (IMAGE_SIZES, KERNEL_SIZES, BOUNDARIES)
    param_names = ["image_size", "kernel_size", "boundary"]

    def setup(self, image_size, kernel_size, boundary):
        rng = np.random.default_rng(12345)
        self.array = rng.random((image_size, image_size))
        k1d = _make_gaussian_1d(kernel_size)
        self.kernel = np.outer(k1d, k1d)
        self.boundary = boundary
        self.convolve_mod = importlib.import_module("astropy.convolution.convolve")

    def time_convolve_separable_fast(self, image_size, kernel_size, boundary):
        convolve(self.array, self.kernel, boundary=boundary)

    def time_convolve_separable_baseline(self, image_size, kernel_size, boundary):
        with _disable_separable_fast_path(self.convolve_mod):
            convolve(self.array, self.kernel, boundary=boundary)
