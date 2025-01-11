import numpy as np
import pyarrow as pa

from astropy.table import Table


def convert_pa_string_array_to_numpy(str_arr: pa.ChunkedArray) -> np.ndarray:
    # First get the maximum length of the strings
    # char_lengths = pc.utf8_length(string_array)
    # max_chars = pc.max(char_lengths).as_py()
    # dtype = f"<U{max_chars}"

    # Check if string_array has any nulls
    has_null = str_arr.null_count > 0

    # Replace nulls with an empty string
    str_arr_filled = str_arr.fill_null("") if has_null else str_arr

    # Convert to NumPy array with fixed-length Unicode dtype
    np_array = str_arr_filled.to_numpy().astype(str)

    if has_null:
        mask = str_arr.is_null().to_numpy()
        np_array = np.ma.array(np_array, mask=mask, copy=False)

    return np_array


def convert_pa_array_to_numpy(arr):
    if pa.types.is_string(arr.type):
        return convert_pa_string_array_to_numpy(arr)

    # Check if array has any nulls
    has_null = arr.null_count > 0

    # Replace nulls with an empty string
    arr_filled = arr.fill_null(0) if has_null else arr

    # Convert to NumPy array with fixed-length Unicode dtype
    np_array = arr_filled.to_numpy()

    if has_null:
        mask = arr.is_null().to_numpy()
        np_array = np.ma.array(np_array, mask=mask, copy=False)
    return np_array


def convert_pa_table_to_astropy_table(table_pa: pa.Table) -> Table:
    """Convert a PyArrow Table to an Astropy Table."""
    columns = {}
    for name, col in zip(table_pa.column_names, table_pa.itercolumns()):
        col_np = convert_pa_array_to_numpy(col)
        columns[name] = col_np
    out = Table(columns, copy=False)
    return out
