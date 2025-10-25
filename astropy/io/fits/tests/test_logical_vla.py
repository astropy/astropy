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


def test_logical_vla_various_input_formats(tmp_path):
    """Test that logical VLAs accept various input formats like
    'T'/'F' strings, numeric 0/1, booleans, and None."""

    fn = tmp_path / "logical_vla_inputs.fits"

    # Test various input representations
    col = fits.Column(
        name="varied",
        format="PL",
        array=[
            [True, False, None],  # bool and None
            ["T", "F", "t"],  # string characters
            [1, 0, -1],  # numeric (0=False, non-zero=True)
            [b"T", b"F", b"t"],  # bytes
        ],
    )
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # Row 0: [True, False, None]
        assert data["varied"][0].tolist() == [True, False, None]
        # Row 1: 'T', 'F', 't' all map to True, False, True
        assert data["varied"][1].tolist() == [True, False, True]
        # Row 2: 1->True, 0->False, -1->True
        assert data["varied"][2].tolist() == [True, False, True]
        # Row 3: bytes 'T', 'F', 't' -> True, False, True
        assert data["varied"][3].tolist() == [True, False, True]


def test_logical_vla_multi_row(tmp_path):
    """Test logical VLAs with multiple rows and varying lengths."""

    fn = tmp_path / "logical_vla_multirow.fits"

    col = fits.Column(
        name="multi",
        format="PL",
        array=[
            [True],
            [False, True],
            [None, False, True],
            [],  # empty array
        ],
    )
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        assert data["multi"][0].tolist() == [True]
        assert data["multi"][1].tolist() == [False, True]
        assert data["multi"][2].tolist() == [None, False, True]
        assert data["multi"][3].tolist() == []


def test_logical_vla_q_format(tmp_path):
    """Test logical VLAs with Q format (64-bit descriptors)."""

    fn = tmp_path / "logical_vla_q.fits"

    col = fits.Column(name="ql", format="QL", array=[[None, False, True]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        heap = data._get_heap_data()
        heap_list = list(heap.tolist()) if hasattr(heap, "tolist") else list(heap)
        # Should have same byte codes as P format
        assert heap_list == [0, ord("F"), ord("T")]

        v = data["ql"][0]
        assert v.tolist() == [None, False, True]


def test_logical_vla_roundtrip_preserves_values(tmp_path):
    """Test that writing and reading back preserves all logical values."""

    fn = tmp_path / "logical_vla_roundtrip.fits"

    # Create file with multiple rows
    col = fits.Column(
        name="rt",
        format="PL",
        array=[
            [True, False],
            [None, True, False],
            [False, None, True, None],
        ],
    )
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    # Read back and verify all values preserved
    with fits.open(fn) as hdul:
        data = hdul[1].data
        assert data["rt"][0].tolist() == [True, False]
        assert data["rt"][1].tolist() == [None, True, False]
        assert data["rt"][2].tolist() == [False, None, True, None]


def test_logical_vla_empty_string_input(tmp_path):
    """Test that empty strings are treated as NULL (None)."""

    fn = tmp_path / "logical_vla_empty.fits"

    col = fits.Column(name="empty", format="PL", array=[["", "T", "F"]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # Empty string should map to None (0 byte)
        assert data["empty"][0][0] is None
        assert data["empty"][0][1] is True
        assert data["empty"][0][2] is False


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
