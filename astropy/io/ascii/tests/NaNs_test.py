import numpy as np
import pytest
from math import floor, ceil

#solving the issue 17840 which was in case the columns were all nans the program will crash (bug)
def calculate_limits(col):
    col_max = np.nanmax(col)  #first the max(col) was changed to np.nanmax(col) so it will return None in case the column is full of np.nan or None and ignore them in case it has some np.nan or None besides other numbers
    col_min = np.nanmin(col)  #same thing here
    #this handles the case where the column is full of np.nan or None
    if not np.isnan(col_min) and not np.isnan(col_max): 
        return f"[{floor(col_min * 100) / 100.0}/{ceil(col_max * 100) / 100.0}]"
    else:
        #in case all are nan or None an empty string will be returned (depending on the expected output one of the contributors)
        return ''


@pytest.mark.parametrize(
    "col, expected_output",
    [
        ([None, None, None], ""),
        ([np.nan, np.nan, np.nan], ""), 
        ([None, 3, 7, None, 10], "[3.0/10.0]"),
        ([np.nan, 2, 5, np.nan, 8], "[2.0/8.0]"), 
    ],
)
def test_calculate_limits(col, expected_output):
    col = np.array([np.nan if v is None else v for v in col], dtype=np.float64)
    assert calculate_limits(col) == expected_output