from typing import TypeVar

import numpy as np

type _D1 = tuple[int]
type _D2 = tuple[int, int]
type _D3 = tuple[int, int, int]
_D = TypeVar("_D", _D1, _D2, _D3)

def _convolveNd_c(
    result: np.ndarray[_D, np.dtype[np.float64]],
    array_to_convolve: np.ndarray[_D, np.dtype[np.float64]],
    kernel: np.ndarray[_D, np.dtype[np.float64]],
    nan_interpolate: bool,
    embed_result_within_padded_region: bool,
    n_threads: int,
) -> None: ...
