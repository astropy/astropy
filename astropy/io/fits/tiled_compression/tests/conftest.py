import pytest

COMPRESSION_TYPES = [
    "GZIP_1",
    "GZIP_2",
    "RICE_1",
    "HCOMPRESS_1",
    "PLIO_1",
]


def _generate_comp_type_dtype_parameters():
    comp_type_dtype_parameters = []
    for compression_type in COMPRESSION_TYPES:
        for types in {
            # "float",
            "integer",
        }:
            # io.fits doesn't seem able to compress 64-bit data, even though e.g. GZIP_?
            # and HCOMPRESS_1 should be able to handle it.
            itemsizes = {
                "integer": [1, 2, 4],
                "float": [4, 8],
            }[types]
            for itemsize in itemsizes:
                for endian in ["<", ">"]:
                    format = "f" if types == "float" else "i"
                    format = "u" if itemsize == 1 else format
                    comp_type_dtype_parameters.append(
                        (compression_type, f"{endian}{format}{itemsize}")
                    )
    return comp_type_dtype_parameters


@pytest.fixture(
    scope="session",
    params=_generate_comp_type_dtype_parameters(),
    ids=lambda x: f"{x[0]}-{x[1]}",
)
def comp_type_dtype(request):
    return request.param


@pytest.fixture(scope="session")
def compression_type(comp_type_dtype):
    return comp_type_dtype[0]


@pytest.fixture(scope="session")
def dtype(comp_type_dtype):
    return comp_type_dtype[1]
