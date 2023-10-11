import itertools

import numpy as np
import pytest

COMPRESSION_TYPES = [
    "GZIP_1",
    "GZIP_2",
    "RICE_1",
    "HCOMPRESS_1",
    "PLIO_1",
    "NOCOMPRESS",
]


def fitsio_param_to_astropy_param(param):
    # Convert fitsio kwargs to astropy kwargs
    _map = {"qlevel": "quantize_level", "qmethod": "quantize_method"}
    param = {_map[k]: v for k, v in param.items()}

    # Map quantize_level
    if param.get("quantize_level", "missing") is None:
        param["quantize_level"] = 0.0

    return param


def _expand(*params):
    """
    Expands a list of N iterables of parameters into a flat list with all
    combinations of all parameters.
    """
    expanded = []
    for ele in params:
        expanded += list(itertools.product(*ele))
    return expanded


ALL_INTEGER_DTYPES = [
    "".join(ele)
    for ele in _expand([("<", ">"), ("i",), ("2", "4")], [("<", ">"), ("u",), ("1",)])
]
ALL_FLOAT_DTYPES = ["".join(ele) for ele in _expand([("<", ">"), ("f",), ("4", "8")])]


@pytest.fixture(
    scope="session",
    ids=lambda x: " ".join(map(str, x)),
    # The params here are compression type, parameters for the compression /
    # quantise and dtype
    params=_expand(
        # Test all compression types with default compression parameters for
        # all integers
        [
            COMPRESSION_TYPES,
            ({},),
            ALL_INTEGER_DTYPES,
        ],
        # GZIP and NOCOMPRESS support lossless non-quantized floating point data
        [
            ("GZIP_1", "GZIP_2", "NOCOMPRESS"),
            ({"qlevel": None},),
            ALL_FLOAT_DTYPES,
        ],
        # All compression types can also take quantized floating point input
        # Rather than running all quantization parameters for all algorithms
        # split up the algorithms to reduce the total number of tests.
        [
            ["GZIP_1", "GZIP_2"],
            ({"qlevel": 5, "qmethod": -1},),
            ALL_FLOAT_DTYPES,
        ],
        [
            ["RICE_1"],
            ({"qlevel": 10, "qmethod": 1},),
            ALL_FLOAT_DTYPES,
        ],
        [
            ["HCOMPRESS_1"],
            (
                {"qlevel": 20, "qmethod": 2},
                {"qlevel": 10, "qmethod": 1},
            ),
            ALL_FLOAT_DTYPES,
        ],
        # Note no PLIO here as that's intended for masks, i.e. data which can't
        # be generated with quantization.
    ),
)
def comp_param_dtype(request):
    return request.param


@pytest.fixture(scope="session")
def compression_type(comp_param_dtype):
    return comp_param_dtype[0]


@pytest.fixture(scope="session")
def compression_param(comp_param_dtype):
    return comp_param_dtype[1]


@pytest.fixture(scope="session")
def dtype(comp_param_dtype):
    return comp_param_dtype[2]


@pytest.fixture(scope="session")
def numpy_rng():
    return np.random.default_rng(0)
