from astropy.nddata import NDData, NDDataRef, NDIOMixin  # noqa: F401

# Alias NDDataAllMixins in case this will be renamed ... :-)
NDDataIO = NDDataRef


def test_simple_write_read():
    ndd = NDDataIO([1, 2, 3])
    assert hasattr(ndd, "read")
    assert hasattr(ndd, "write")
