"""Benchmark separable-kernel convolution speedups.

This script compares baseline (forced non-separable path) against the
separable fast path for a grid of image/kernel sizes and prints a
Markdown table.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib
import time
from collections.abc import Iterable

import numpy as np

from astropy.convolution import convolve


@contextlib.contextmanager
def _disable_separable_fast_path():
    convolve_mod = importlib.import_module("astropy.convolution.convolve")
    original = convolve_mod._factor_separable_kernel_2d
    convolve_mod._factor_separable_kernel_2d = lambda _kernel: None
    try:
        yield
    finally:
        convolve_mod._factor_separable_kernel_2d = original


def _make_gaussian_1d(size: int, sigma: float | None = None) -> np.ndarray:
    if size % 2 == 0:
        raise ValueError("Kernel size must be odd.")
    if sigma is None:
        sigma = size / 6.0
    x = np.arange(size, dtype=float) - size // 2
    kernel = np.exp(-(x * x) / (2.0 * sigma * sigma))
    kernel /= kernel.sum()
    return kernel


def _time_convolve(
    array: np.ndarray,
    kernel: np.ndarray,
    boundary: str,
    repeats: int,
) -> float:
    # Warm-up
    convolve(array, kernel, boundary=boundary)
    timings = []
    for _ in range(repeats):
        start = time.perf_counter()
        convolve(array, kernel, boundary=boundary)
        timings.append(time.perf_counter() - start)
    return min(timings)


def _parse_sizes(values: Iterable[str]) -> list[int]:
    sizes: list[int] = []
    for value in values:
        for chunk in value.split(","):
            chunk = chunk.strip()
            if not chunk:
                continue
            sizes.append(int(chunk))
    return sizes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--image-sizes",
        nargs="*",
        default=["128", "256", "512"],
        help="Image sizes (comma-separated or repeated flags).",
    )
    parser.add_argument(
        "--kernel-sizes",
        nargs="*",
        default=["5", "11", "21", "31"],
        help="Odd kernel sizes (comma-separated or repeated flags).",
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--boundary", default="fill")
    parser.add_argument("--seed", type=int, default=0)

    args = parser.parse_args()

    image_sizes = _parse_sizes(args.image_sizes)
    kernel_sizes = _parse_sizes(args.kernel_sizes)

    rng = np.random.default_rng(args.seed)

    header = [
        "| image | kernel | baseline (ms) | fast (ms) | speedup |",
        "|---:|---:|---:|---:|---:|",
    ]
    rows = []

    for image_size in image_sizes:
        array = rng.random((image_size, image_size))
        for kernel_size in kernel_sizes:
            if kernel_size % 2 == 0:
                raise ValueError("Kernel sizes must be odd.")
            k1d = _make_gaussian_1d(kernel_size)
            kernel = np.outer(k1d, k1d)

            with _disable_separable_fast_path():
                baseline = _time_convolve(array, kernel, args.boundary, args.repeats)

            fast = _time_convolve(array, kernel, args.boundary, args.repeats)

            speedup = baseline / fast if fast else float("inf")
            rows.append(
                f"| {image_size} | {kernel_size} | {baseline * 1000.0:.3f} | {fast * 1000.0:.3f} | {speedup:.2f}x |"
            )

    print("\n".join(header + rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
