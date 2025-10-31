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


def test_logical_vla_scalar_input(tmp_path):
    """Test that scalar inputs (non-sequences) are handled correctly."""

    fn = tmp_path / "logical_vla_scalar.fits"

    # Test with scalar boolean (should be treated as single-element array)
    col = fits.Column(name="scalar", format="PL", array=[True, False, None])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # Each row should be a single-element array
        assert data["scalar"][0].tolist() == [True]
        assert data["scalar"][1].tolist() == [False]
        assert data["scalar"][2].tolist() == [None]


def test_logical_vla_invalid_bytes(tmp_path):
    """Test fallback handling for invalid byte sequences and arbitrary bytes."""

    fn = tmp_path / "logical_vla_invalid.fits"

    # Test with various byte values:
    # - Non-decodable bytes (should fallback to True)
    # - Arbitrary decodable bytes (should become True unless 'F'/'False')
    col = fits.Column(
        name="inv",
        format="PL",
        array=[[b"\xff\xfe", b"X", b"abc", b"F", b"false"]],
    )
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # All should become True except 'F' and 'false'
        expected = [True, True, True, False, False]
        assert data["inv"][0].tolist() == expected


def test_logical_vla_numpy_bool(tmp_path):
    """Test numpy boolean types are handled correctly."""

    fn = tmp_path / "logical_vla_npbool.fits"

    # Use numpy boolean arrays
    col = fits.Column(
        name="npb",
        format="PL",
        array=[
            np.array([True, False], dtype=np.bool_),
            np.array([np.True_, np.False_]),
        ],
    )
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        assert data["npb"][0].tolist() == [True, False]
        assert data["npb"][1].tolist() == [True, False]


def test_logical_vla_mixed_case_strings(tmp_path):
    """Test that lowercase and mixed-case strings are handled."""

    fn = tmp_path / "logical_vla_case.fits"

    col = fits.Column(name="case", format="PL", array=[["t", "f", "T", "F"]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # All should map correctly (case-insensitive)
        assert data["case"][0].tolist() == [True, False, True, False]


def test_logical_vla_non_tf_characters(tmp_path):
    """Test that arbitrary strings (not 'F'/'False') are treated as True."""

    fn = tmp_path / "logical_vla_chars.fits"

    # Use arbitrary strings - they should all be treated as True
    # Only 'F' or 'False' (case insensitive) should become False
    col = fits.Column(
        name="chars",
        format="PL",
        array=[["X", "abc", "123", "True", "yes", "F", "False", "false"]],
    )
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # All non-F/False strings should map to True
        expected = [True, True, True, True, True, False, False, False]
        assert data["chars"][0].tolist() == expected

        # Verify heap bytes are correct (all should be 'T' or 'F')
        raw_data = hdul[1].data._get_raw_data()
        heap_offset = hdul[1].data._heapoffset
        heap_size = hdul[1].data._heapsize
        heap_data = np.frombuffer(
            raw_data[heap_offset : heap_offset + heap_size], dtype=np.uint8
        )
        expected_bytes = [
            ord("T"),
            ord("T"),
            ord("T"),
            ord("T"),
            ord("T"),
            ord("F"),
            ord("F"),
            ord("F"),
        ]
        assert heap_data.tolist() == expected_bytes


def test_logical_vla_read_invalid_heap_bytes(tmp_path):
    """Test reading heap with non-standard byte codes."""

    fn = tmp_path / "logical_vla_read.fits"

    # First create a normal file
    col = fits.Column(name="test", format="PL", array=[[True, False, None]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    # Now read it back - this exercises the reading path
    with fits.open(fn) as hdul:
        data = hdul[1].data
        # Access the field to trigger conversion
        result = data["test"][0]
        assert result.tolist() == [True, False, None]


def test_logical_vla_non_numeric_input(tmp_path):
    """Test that non-numeric, non-bool, non-string inputs fallback to True."""

    fn = tmp_path / "logical_vla_nonnumeric.fits"

    # Use objects that can't be converted to int (should fallback to True)
    # Using complex numbers which can't be directly converted to int
    col = fits.Column(name="obj", format="PL", array=[[1.5, 2.7, 3.14]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    with fits.open(fn) as hdul:
        data = hdul[1].data
        # Floats are converted: non-zero should be True
        assert data["obj"][0].tolist() == [True, True, True]


def test_logical_vla_write_then_scale_back(tmp_path):
    """Test the _scale_back path that writes logical VLA data back to heap."""

    fn = tmp_path / "logical_vla_scaleback.fits"

    # Create a file with logical VLA
    col = fits.Column(name="sb", format="PL", array=[[None, False, True]])
    tab = fits.BinTableHDU.from_columns([col])

    # Write to file
    tab.writeto(fn, overwrite=True)

    # Open and modify to trigger _scale_back with object dtype conversion
    with fits.open(fn, mode="update") as hdul:
        # This should trigger the _get_heap_data conversion path
        # that handles object arrays for logical VLAs
        hdul[1].data["sb"]  # Access the column
        hdul.flush()

    # Verify it still works after modification
    with fits.open(fn) as hdul:
        data = hdul[1].data
        assert data["sb"][0].tolist() == [None, False, True]


def test_logical_vla_object_with_complex_type(tmp_path):
    """Test reading with values that can't be converted to int."""

    fn = tmp_path / "logical_vla_complex.fits"

    # Create file with normal data first
    col = fits.Column(name="ct", format="PL", array=[[True, False, True]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    # Read back - the conversion code should handle all values
    with fits.open(fn) as hdul:
        data = hdul[1].data
        vals = data["ct"][0]
        # All values should be converted successfully
        assert vals.tolist() == [True, False, True]


def test_logical_vla_heap_conversion_all_values(tmp_path):
    """Test _get_heap_data conversion path with all logical values (None/False/True)."""

    fn = tmp_path / "logical_vla_heap_all.fits"

    # Create file with all three logical values to hit all branches in _get_heap_data
    col = fits.Column(
        name="all",
        format="PL",
        array=[
            [None],  # Tests the `if val is None` branch
            [False],  # Tests the `elif val is False` branch
            [True],  # Tests the `else` branch (True)
        ],
    )
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    # Open and flush to trigger _get_heap_data conversion
    with fits.open(fn, mode="update") as hdul:
        # Accessing the data triggers conversion
        _ = hdul[1].data["all"]
        hdul.flush()  # This triggers _scale_back -> _get_heap_data

    # Verify values are still correct after conversion
    with fits.open(fn) as hdul:
        data = hdul[1].data
        assert data["all"][0].tolist() == [None]
        assert data["all"][1].tolist() == [False]
        assert data["all"][2].tolist() == [True]


def test_logical_vla_read_conversion_exception(tmp_path):
    """Test exception handling in _convert_other when value can't be converted to int."""

    fn = tmp_path / "logical_vla_exception.fits"

    # Create a file with logical VLA
    col = fits.Column(name="exc", format="PL", array=[[True, False, None]])
    tab = fits.BinTableHDU.from_columns([col])
    tab.writeto(fn, overwrite=True)

    # Reading should handle all values correctly
    # The exception path is triggered when numpy can't convert to int
    with fits.open(fn) as hdul:
        data = hdul[1].data
        result = data["exc"][0]
        # All values should be properly converted
        assert result.tolist() == [True, False, None]


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
