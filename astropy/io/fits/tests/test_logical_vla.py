import numpy as np

from astropy.io import fits


def test_write_read_logical_vla_heap_and_values(tmp_path):
    """Writing a PL (logical VLA) column should store ASCII codes
    0 (NULL), 70 ('F'), 84 ('T') in the heap and read back as
    [None, False, True]."""

    fn = tmp_path / "logical_vla.fits"

    col = fits.Column(name="a", format="PL", array=[[None, False, True]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # Access raw heap bytes for inspection
        heap = data._get_heap_data()
        # heap is a numpy byte array: convert to list of ints
        heap_list = list(heap.tolist()) if hasattr(heap, "tolist") else list(heap)

        # Expect the three bytes [0, ord('F'), ord('T')] in that order
        assert heap_list == [0, ord("F"), ord("T")]

        # Data should read back to an object array with [None, False, True]
        v = data["a"][0]
        assert isinstance(v, np.ndarray)
        # Convert to Python list for equality check
        assert v.tolist() == [None, False, True]


def test_fixed_width_logical_roundtrip(tmp_path):
    """Ensure fixed-width logical columns (TFORM='L') still roundtrip
    as boolean arrays (regression test)."""

    fn = tmp_path / "logical_fixed.fits"

    col = fits.Column(name="b", format="L", array=[True, False, True])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        assert data["b"].dtype == np.dtype("bool") or data["b"].dtype == bool
        assert data["b"].tolist() == [True, False, True]
