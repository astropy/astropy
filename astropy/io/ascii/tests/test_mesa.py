# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for MESA stellar evolution code file format reader/writer.

Test Data Files
---------------
Test data files are located in astropy/io/ascii/tests/data/:

- history_mesa.data : MESA history file with 25 models and 3 columns
  (truncated from a real stellar evolution calculation with no restarts)

- history_mesa_with_restarts.data : MESA history file with restart artifacts,
  containing 31 rows that should be cleaned to 20 rows after
  removing the outdated models 10-20 from the first run before restart

- profile_mesa.data : MESA profile file with minimal zones and columns

- history_mesa_bad_format.data : File with bad format (missing metadata lines)

- history_mesa_missing_columns.data : File missing required 'model_number' or
  `zone` column

"""
import io

import numpy as np
import pytest

from astropy.io import ascii
from astropy.utils.data import get_pkg_data_filename

# Minimal MESA history file for testing
# Note: MESA format requires column indices on both metadata and data sections
MESA_HISTORY = """                           1                    2                    3
              version_number             compiler                build
                  "r24.03.1"           "gfortran"             "13.1.0"

                           1                    2                    3                    4                    5
                model_number            num_zones             star_age                 log_L                log_Teff
                           1                 1004   1.00000000000E-05   1.23456789012E+00   3.65432109876E+00
                           2                 1025   2.00000000000E-05   1.23456789013E+00   3.65432109877E+00
                           3                 1050   3.00000000000E-05   1.23456789014E+00   3.65432109878E+00
"""

# MESA file with restart artifacts
MESA_WITH_RESTART = """                           1                    2
              version_number             compiler
                  "r24.03.1"           "gfortran"

                           1                    2                    3
                model_number            num_zones             star_age
                           1                 1004   1.00000000000E-05
                           2                 1005   2.00000000000E-05
                           3                 1006   3.00000000000E-05
                           4                 1007   4.00000000000E-05
                           5                 1008   5.00000000000E-05
                           3                 2000   3.50000000000E-05
                           5                 2002   5.50000000000E-05
                           7                 2004   7.50000000000E-05
"""

# MESA profile file (same format as history)
MESA_PROFILE = """                           1                    2
              version_number             compiler
                  "r23.05.1"              "ifort"

                           1                    2                    3                    4
                        zone                 mass               radius            log_rho
                           1   1.00000000000E-03   1.00000000000E+10   5.00000000000E+00
                           2   2.00000000000E-03   2.00000000000E+10   4.90000000000E+00
                           3   3.00000000000E-03   3.00000000000E+10   4.80000000000E+00
"""


class TestMesaReader:
    """Tests for reading MESA format files."""

    def test_read_basic(self):
        """Test basic reading of MESA file."""
        table = ascii.read(MESA_HISTORY, format="mesa")
        assert len(table) == 3
        assert len(table.columns) == 5
        assert table.colnames == [
            "model_number",
            "num_zones",
            "star_age",
            "log_L",
            "log_Teff",
        ]

    def test_read_metadata(self):
        """Test metadata extraction from MESA file."""
        table = ascii.read(MESA_HISTORY, format="mesa")
        assert "header" in table.meta
        header = table.meta["header"]
        assert header["version_number"] == "r24.03.1"
        assert header["compiler"] == "gfortran"
        assert header["build"] == "13.1.0"

    def test_read_data_values(self):
        """Test that data values are read correctly."""
        table = ascii.read(MESA_HISTORY, format="mesa")
        assert table["model_number"][0] == 1
        assert table["model_number"][2] == 3
        assert table["num_zones"][0] == 1004
        assert np.isclose(table["star_age"][0], 1.0e-5)
        assert np.isclose(table["log_L"][1], 1.23456789013)

    def test_read_data_types(self):
        """Test that data types are preserved (int vs float)."""
        table = ascii.read(MESA_HISTORY, format="mesa")
        # Integer columns should be int
        assert table["model_number"].dtype.kind == "i"
        assert table["num_zones"].dtype.kind == "i"
        # Float columns should be float
        assert table["star_age"].dtype.kind == "f"
        assert table["log_L"].dtype.kind == "f"
        assert table["log_Teff"].dtype.kind == "f"

    def test_read_profile(self):
        """Test reading MESA profile file (same format)."""
        table = ascii.read(MESA_PROFILE, format="mesa")
        assert len(table) == 3
        assert table.colnames[0] == "zone"  # Profiles have 'zone' not 'model_number'
        assert table["zone"][0] == 1
        assert table.meta["header"]["version_number"] == "r23.05.1"
        assert table.meta["header"]["compiler"] == "ifort"

    def test_restart_removal_default(self):
        """Test that restart artifacts are removed by default."""
        table = ascii.read(MESA_WITH_RESTART, format="mesa")
        # Should have 7 rows (3 removed from first run)
        assert len(table) == 5
        # Model numbers should be monotonically increasing
        assert all(table["model_number"][:-1] < table["model_number"][1:])
        # Models should be [1, 2, 3, 4, 5, 6, 7]
        assert list(table["model_number"]) == [1, 2, 3, 5, 7]

    def test_restart_removal_disabled(self):
        """Test reading with restart removal disabled."""
        table = ascii.read(MESA_WITH_RESTART, format="mesa", remove_restart_rows=False)
        # Should have all 10 rows
        assert len(table) == 8
        # Model numbers will have duplicates
        assert list(table["model_number"]) == [1, 2, 3, 4, 5, 3, 5, 7]

    def test_restart_removal_algorithm(self):
        """Test restart removal algorithm with specific scenario."""
        table = ascii.read(MESA_WITH_RESTART, format="mesa")

        # The restart happened at model 3
        # Original: 1, 2, 3, 4, 5 (first run), then 3, 4, 5, 6, 7 (restart from 3)
        # Removed: models at indices 2, 3, 4 (first occurrence of 3, 4, 5)
        # Result: models 1, 2, 3, 4, 5, 6, 7 (from restart)

        assert len(table) == 5

        # Check num_zones to verify we kept the right data
        # After restart, num_zones changes from 100x to 200x
        assert table["num_zones"][0] == 1004  # model 1
        assert table["num_zones"][1] == 1005  # model 2
        assert table["num_zones"][2] == 2000  # model 3 (from restart, not first run)
        assert table["num_zones"][3] == 2002  # model 5 (from restart)
        assert table["num_zones"][4] == 2004  # model 7


class TestMesaAutoDetection:
    """Tests for MESA file format auto-detection."""

    @classmethod
    def setup_class(cls):
        """Set up test fixtures once for the entire class."""
        cls.good_history_file_obj = open(
            get_pkg_data_filename("data/history_mesa.data")
        )
        cls.good_profile_file_obj = open(
            get_pkg_data_filename("data/profile_mesa.data")
        )
        cls.bad_format_file_obj = open(
            get_pkg_data_filename("data/history_mesa_bad_format.data")
        )
        cls.bad_columns_file_obj = open(
            get_pkg_data_filename("data/history_mesa_missing_columns.data")
        )

    @classmethod
    def teardown_class(cls):
        """Clean up file handles after tests."""
        cls.good_history_file_obj.close()
        cls.good_profile_file_obj.close()
        cls.bad_format_file_obj.close()
        cls.bad_columns_file_obj.close()

    def test_autodetect_history(self):
        """Test auto-detection of history files."""
        from astropy.io.ascii.mesa import mesa_identify

        # Access using self
        assert mesa_identify("read", "history.data", self.good_history_file_obj) is True
        self.good_history_file_obj.seek(0)
        assert (
            mesa_identify("read", "history_backup.data", self.good_history_file_obj)
            is True
        )
        self.good_history_file_obj.seek(0)
        assert mesa_identify("read", "HISTORY.DATA", self.good_history_file_obj) is True
        self.good_history_file_obj.seek(0)

    def test_autodetect_profile(self):
        """Test auto-detection of profile files."""
        from astropy.io.ascii.mesa import mesa_identify

        assert (
            mesa_identify("read", "profile1.data", self.good_profile_file_obj) is True
        )
        self.good_profile_file_obj.seek(0)
        assert (
            mesa_identify("read", "PROFILE99.DATA", self.good_profile_file_obj) is True
        )
        self.good_profile_file_obj.seek(0)

    # sad path tests; bad file names, bad format, and/or missing columns
    def test_autodetect_negative(self):
        """Test that non-MESA files are not identified."""
        from astropy.io.ascii.mesa import mesa_identify

        # bad names (not starting with history/profile and ending with .data)
        assert mesa_identify("read", "data.csv", self.good_history_file_obj) is False
        self.good_history_file_obj.seek(0)
        assert mesa_identify("read", "random.txt", self.good_history_file_obj) is False
        self.good_history_file_obj.seek(0)
        assert mesa_identify("read", "mesa.dat", self.good_history_file_obj) is False
        self.good_history_file_obj.seek(0)
        assert mesa_identify("read", None, self.good_history_file_obj) is False
        self.good_history_file_obj.seek(0)

        # good name, but bad format
        assert mesa_identify("read", "history.data", self.bad_format_file_obj) is False
        self.bad_format_file_obj.seek(0)

        # good name, good format, but doesn't have zone or model_number columns
        assert mesa_identify("read", "history.data", self.bad_columns_file_obj) is False
        self.bad_columns_file_obj.seek(0)

        # a few more cases that look like MESA file names but are not
        assert mesa_identify("read", "data.csv", None) is False

        text_lt_7_lines = io.StringIO("             1                    2\n")
        text_lt_7_lines.write("version_number             compiler\n")
        text_lt_7_lines.seek(0)
        assert mesa_identify("read", "data.csv", text_lt_7_lines) is False


class TestMesaEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_file(self):
        """Test handling of empty or too-short files."""
        with pytest.raises((ascii.InconsistentTableError, IndexError)):
            ascii.read("   \n  \n  ", format="mesa")

    def test_single_row(self):
        """Test reading file with only one data row."""
        single_row = """                           1
              version_number
                  "r24.03.1"

                           1                    2
                model_number             star_age
                           1   1.00000000000E-05
"""
        table = ascii.read(single_row, format="mesa")
        assert len(table) == 1
        assert table["model_number"][0] == 1


class TestMesaIntegration:
    """Integration tests using actual test data files."""

    def test_read_mesa_history(self):
        """Test reading the mesa_history.data test file."""
        test_file = get_pkg_data_filename("data/history_mesa.data")

        table = ascii.read(test_file, format="mesa")

        # Check basic properties
        assert len(table) == 25
        assert len(table.columns) == 3
        assert table.colnames[0] == "model_number"

        # Check metadata
        assert table.meta["header"]["version_number"] == "r24.03.1"
        assert table.meta["header"]["compiler"] == "gfortran"

        # Check monotonicity
        model_nums = table["model_number"]
        assert all(model_nums[:-1] < model_nums[1:])

        # Check first and last values
        assert table["model_number"][0] == 1
        assert table["model_number"][-1] == 25

    def test_read_history_mesa_with_restarts(self):
        """Test reading history_mesa_with_restarts.data and verifying restart removal."""
        test_file = get_pkg_data_filename("data/history_mesa_with_restarts.data")

        # Read without restart removal
        table_full = ascii.read(test_file, format="mesa", remove_restart_rows=False)

        # Read with restart removal (default)
        table_cleaned = ascii.read(test_file, format="mesa")

        # Full table should have 31 rows
        assert len(table_full) == 31

        # Cleaned table should have 30 rows (11 removed: duplicate/outdated
        # models 10-20)
        assert len(table_cleaned) == 20

        # Cleaned data should be monotonic
        assert all(
            table_cleaned["model_number"][:-1] < table_cleaned["model_number"][1:]
        )

        # Check specific pattern: restart from model 10, should have 1-10, then
        # 12, 14, ..., 30
        assert list(table_cleaned["model_number"]) == list(range(1, 11)) + list(
            range(12, 31, 2)
        )
