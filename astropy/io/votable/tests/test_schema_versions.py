import re
import sys

import pytest

from astropy.io.votable.xmlutil import validate_schema
from astropy.utils.data import get_pkg_data_filename

# Each test file uses a feature that is new with its version.
V_1_2_TABLE = get_pkg_data_filename("data/empty_table.xml")
V_1_3_TABLE = get_pkg_data_filename("data/binary2_masked_strings.xml")
V_1_4_TABLE = get_pkg_data_filename("data/timesys.xml")
V_1_5_TABLE = get_pkg_data_filename("data/coosys.xml")

# Schema versions to test against
V_1_2 = "1.2"
V_1_3 = "1.3"
V_1_4 = "1.4"
V_1_5 = "1.5"
V_1_6 = "1.6"  # Next (unimplemented) schema

# Starting with version 1.4, schemas were upward compatible.
test_cases = [
    (V_1_2_TABLE, V_1_2, 0, None),
    (
        V_1_2_TABLE,
        V_1_3,
        3,
        b"No matching global declaration available for the validation root",
    ),
    (V_1_3_TABLE, V_1_3, 0, None),
    (V_1_3_TABLE, V_1_4, 0, None),
    (V_1_3_TABLE, V_1_5, 0, None),
    (
        V_1_3_TABLE,
        V_1_6,
        3,
        b"No matching global declaration available for the validation root",
    ),
    (V_1_4_TABLE, V_1_3, 3, b"TIMESYS.*element is not expected"),
    (V_1_4_TABLE, V_1_4, 0, None),
    (V_1_4_TABLE, V_1_5, 0, None),
    (
        V_1_4_TABLE,
        V_1_6,
        3,
        b"No matching global declaration available for the validation root",
    ),
    (V_1_5_TABLE, V_1_3, 3, b"attribute 'refposition' is not allowed"),
    (V_1_5_TABLE, V_1_4, 3, b"attribute 'refposition' is not allowed"),
    (V_1_5_TABLE, V_1_5, 0, None),
    (
        V_1_5_TABLE,
        V_1_6,
        3,
        b"No matching global declaration available for the validation root",
    ),
]


@pytest.mark.skipif(sys.platform.startswith("win"), reason="Cannot test on Windows")
@pytest.mark.parametrize(
    "votable_file,schema_version,expected_return_code,expected_msg_re", test_cases
)
def test_schema_versions(
    votable_file, schema_version, expected_return_code, expected_msg_re
):
    """Test that xmllint gives expected results for the given file and schema version."""
    try:
        rc, stdout, stderr = validate_schema(votable_file, schema_version)
    except OSError:
        # If xmllint is not installed, we want to skip the test
        pytest.skip("xmllint is not available so will not do schema validation")

    assert rc == expected_return_code
    if rc == 3:
        # Schema validation error.  Check the error message content.
        assert re.search(expected_msg_re, stderr)
