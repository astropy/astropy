import pytest

COMPRESSION_TYPES = [
    "GZIP_1",
    "GZIP_2",
    "RICE_1",
    "HCOMPRESS_1",
    "PLIO_1",
]

compression_dtype_parameters = []
for compression_type in COMPRESSION_TYPES:
    # io.fits doesn't seem able to compress 64-bit data, even though e.g. GZIP_?
    # and HCOMPRESS_1 should be able to handle it.
    for itemsize in [1, 2, 4]:
        for endian in ["<", ">"]:
            format = "u" if itemsize == 1 else "i"
            compression_dtype_parameters.append(
                (compression_type, f"{endian}{format}{itemsize}")
            )


@pytest.fixture(
    params=compression_dtype_parameters, scope="module", ids=lambda x: f"{x[0]} {x[1]}"
)
def compression_type_dtype(request):
    """
    This fixture provides a matrix of compression type and dtype parameters.
    """
    return request.param
