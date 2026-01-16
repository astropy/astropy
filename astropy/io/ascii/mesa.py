# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""mesa.py:
  Classes to read MESA stellar evolution code output files (history and profile).

MESA (Modules for Experiments in Stellar Astrophysics) is a widely-used
stellar evolution code. This reader handles the standard MESA history and
profile file formats.

See: https://docs.mesastar.org/en/stable/using_mesa/output.html


:Author: Bill Wolf
"""

import re
from pathlib import Path

import numpy as np

from . import core


def mesa_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if a file is in MESA format.

    MESA files are typically named 'history.data' or 'profileXX.data'.
    This function checks for those patterns, in addition to checking
    the structure of the first seven lines.

    Parameters
    ----------
    origin : str
        Origin of the file
    filepath : str or None
        Path to the file
    fileobj : file-like or None
        File object
    args, kwargs : tuple, dict
        Additional arguments

    Returns
    -------
    bool
        True if the file matches MESA naming patterns
    """
    if filepath is None:
        return False

    # Extract just the filename from the path

    filename = Path(filepath).name

    # Check for MESA history or profile patterns
    # history.data, history_something.data, profileXX.data, profile.data, etc.
    history_pattern = r"^history.*\.data$"
    profile_pattern = r"^profile.*\.data$"

    file_name_check = (
        re.match(history_pattern, filename, re.IGNORECASE) is not None
        or re.match(profile_pattern, filename, re.IGNORECASE) is not None
    )

    if not file_name_check:
        return False

    # file name seems legitimate, so check the first few lines for MESA structure
    if fileobj is None:
        return False

    # Read the first 7 lines to check for MESA structure
    lines = []
    for _ in range(7):
        line = fileobj.readline()
        lines.append(line)

    # Check if the first 7 lines match the expected MESA format
    if len(lines) < 7:
        return False

    # Check that line 0 (0-indexed) contains integers separated by whitespace
    meta_col_indices_line = lines[0]
    if not re.match(r"^\s*\d+\s+\d+", meta_col_indices_line):
        return False

    # count number of metadata columns and ensure that line 1 (0-indexed) has same number of names
    num_meta_cols = len(meta_col_indices_line.split())
    meta_col_names_line = lines[1]
    if len(meta_col_names_line.split()) != num_meta_cols:
        return False

    # Check that line 2 (0-indexed) contains metadata values
    meta_values_line = lines[2]
    if len(meta_values_line.split()) != num_meta_cols:
        return False

    # Similar checks for data table columns
    data_col_indices_line = lines[4]
    if not re.match(r"^\s*\d+\s+\d+", data_col_indices_line):
        return False
    num_data_cols = len(data_col_indices_line.split())
    data_col_names_line = lines[5]
    if len(data_col_names_line.split()) != num_data_cols:
        return False

    # check first line of table data (line 6, 0-indexed) has same number of entries as data columns
    data_row_line = lines[6]
    if len(data_row_line.split()) != num_data_cols:
        return False

    # at this point, file names look like a MESA file, and the first 7 lines have the expected structure
    # of a MESA file, so it is likely a MESA file. For final check, check to see if "model_number" or
    # "zone" is in the data column names, as these are very likely to be in a history or profile file, respectively.
    data_col_names = data_col_names_line.split()
    if "model_number" not in data_col_names and "zone" not in data_col_names:
        return False

    return True


class MesaHeader(core.BaseHeader):
    """Reader for MESA file headers.

    MESA files have a distinctive format:
    - Line 0: Column indices for metadata
    - Line 1: Metadata column names
    - Line 2: Metadata values
    - Line 3: Blank line
    - Line 4: Column indices for data table
    - Line 5: Column names for data table
    - Line 6+: Data rows
    """

    start_line = 5  # Column names are on line 5 (0-indexed, in original file)
    comment = None  # No comment character in MESA files

    def get_cols(self, lines):
        """Initialize the header Column objects from MESA file lines.

        Note: get_cols receives the ORIGINAL lines (not processed), so we use
        the original line indices:
        - Line 0: metadata column indices
        - Line 1: metadata column names
        - Line 2: metadata values
        - Line 3: blank line
        - Line 4: data column indices
        - Line 5: data column names <- we want this

        Parameters
        ----------
        lines : list
            List of table lines (original, not processed)
        """
        # Convert lines to a list so we can index it
        lines_list = list(lines)

        if len(lines_list) < 6:
            raise core.InconsistentTableError(
                "MESA file must have at least 6 lines (metadata + header + data)"
            )

        # Line 5 (0-indexed) contains the column names
        col_names_line = lines_list[5]

        # Split the line to get column names
        self.names = col_names_line.split()

        self._set_cols_from_names()

    def update_meta(self, lines, meta):
        """Extract MESA metadata from the file header.

        MESA files contain metadata in the first three lines:
        - Line 1: Metadata column names
        - Line 2: Metadata values

        These are stored in meta['table']['mesa_metadata'] as a dictionary.
        """
        # Get lines as a list to access by index
        lines_list = list(lines)

        if len(lines_list) < 3:
            return

        # Parse metadata from lines 1 and 2
        # Line 1: metadata column names
        # Line 2: metadata values
        meta_names_line = lines_list[1]
        meta_values_line = lines_list[2]

        # Split the lines using whitespace
        meta_names = meta_names_line.split()
        meta_values = meta_values_line.split()

        # Create metadata dictionary
        mesa_meta = {}
        for name, value in zip(meta_names, meta_values):
            # Strip quotes from string values
            if value.startswith('"') and value.endswith('"'):
                value = value.strip('"')
            # Try to convert to number if possible
            else:
                try:
                    # Try integer first
                    value = int(value)
                except ValueError:
                    try:
                        # Try float
                        value = float(value)
                    except ValueError:
                        # Keep as string
                        pass

            mesa_meta[name] = value

        # Store in table metadata
        meta.setdefault("table", {})["header"] = mesa_meta


class MesaData(core.BaseData):
    """Reader for MESA data section.

    Data starts at line 5 (0-indexed) after blank lines are removed.
    After process_lines removes the blank line 3, the structure is:
    - Line 0: metadata column indices
    - Line 1: metadata column names
    - Line 2: metadata values
    - Line 3: data column indices (line 4 in original file)
    - Line 4: data column names (line 5 in original file)
    - Line 5: first data row (line 6 in original file)
    """

    start_line = 5
    comment = None  # No comment character
    delimiter = None  # Whitespace-delimited


class Mesa(core.BaseReader):
    """MESA stellar evolution code output format.

    MESA (Modules for Experiments in Stellar Astrophysics) is a widely-used
    one-dimensional stellar evolution code. This reader handles both history
    files (evolution over time) and profile files (spatial profiles at a
    single timestep). It does NOT handle model (.mod) files or output from
    associated packages like GYRE, adipls, or Stella.

    The MESA format consists of:

    - Lines 0-2: Metadata section with column indices, names, and values
    - Line 3: Blank line
    - Line 4: Column indices for main data table
    - Line 5: Column names for main data table
    - Line 6+: Data rows (whitespace-delimited)

    The metadata from lines 1-2 is stored in the table's ``meta['header']``
    dictionary.

    For history files, this reader automatically detects and removes restart
    artifacts. When a MESA run is restarted from an earlier model, the history
    file will contain duplicate model numbers. This reader removes the earlier
    instances so that the resulting table has monotonically increasing model
    numbers suitable for analysis and plotting.

    Example::

                                           1                    2
                              version_number             compiler
                                  "r24.03.1"           "gfortran"

                                           1                    2                    3
                                model_number            num_zones             star_age
                                           1                 1004   1.00000000000E-05
                                           2                 1025   2.20000000000E-05

    Example usage::

    >>> from astropy.table import Table
    >>> from astropy.utils.data import get_pkg_data_filename
    >>> # Find the location of an example file included in astropy
    >>> mesa_file = get_pkg_data_filename('data/history_mesa.data',
    ...                                      package='astropy.io.ascii.tests')
    >>> t = Table.read(mesa_file, format='ascii.mesa')
    >>> print(t.meta['header']['version_number'])
    r24.03.1
    >>> t['model_number', 'star_age'][:3]
    <Table length=3>
    model_number        star_age
       int64            float64
    ------------ ----------------------
               1                  1e-05
               2                2.2e-05
               3 3.6400000000000004e-05
    """

    _format_name = "mesa"
    _description = "MESA stellar evolution code output"
    _io_registry_can_write = False

    header_class = MesaHeader
    data_class = MesaData

    def __init__(self, remove_restart_rows=True):
        """Initialize MESA reader.

        Parameters
        ----------
        remove_restart_rows : bool, optional
            If True (default), automatically remove restart artifacts from history
            files by detecting non-monotonic model numbers. Set to False to keep
            all rows including restart artifacts.
        """
        super().__init__()
        self.remove_restart_rows = remove_restart_rows

    def read(self, table):
        """Read input data into a Table and return the result.

        Parameters
        ----------
        table : str, file-like, list
            Input table data

        Returns
        -------
        out : `~astropy.table.Table`
            Output table with restart artifacts removed (if applicable)
        """
        out = super().read(table)

        # Check if this is a history file and remove restarts if requested
        if self.remove_restart_rows and "model_number" in out.colnames:
            out = self._remove_restart_rows(out)

        return out

    def _remove_restart_rows(self, table):
        """
        Remove rows from history files that correspond to restarts.

        When MESA restarts from an earlier model, the history file contains
        duplicate model numbers. We iterate backwards through the model_number
        column and remove rows where the model number is >= the minimum seen
        so far (going backwards). This identifies the "old" evolution that was
        restarted.

        Example:
            Model numbers: [1, 2, 3, 4, 5, 3, 4, 5, 6, 7]
            Going backwards:
            - i=9: model=7, min_seen=7, keep
            - i=8: model=6, min_seen=6, keep
            - i=7: model=5, min_seen=5, keep
            - i=6: model=4, min_seen=4, keep
            - i=5: model=3, min_seen=3, keep
            - i=4: model=5, 5 >= 3, REMOVE
            - i=3: model=4, 4 >= 3, REMOVE
            - i=2: model=3, 3 >= 3, REMOVE
            - i=1: model=2, 2 < 3, min_seen=2, keep
            - i=0: model=1, 1 < 2, min_seen=1, keep
            Result: [1, 2, 3, 4, 5, 6, 7]

        Parameters
        ----------
        table : `~astropy.table.Table`
            Input table with potential restart artifacts

        Returns
        -------
        table : `~astropy.table.Table`
            Table with restart artifacts removed
        """
        model_number = table["model_number"]
        # are looking at the reversed array to get the last occurrence
        smallest_seen = np.minimum.accumulate(model_number[::-1])
        _, reverse_index = np.unique(smallest_seen[:], return_index=True)
        # Transform the index to work for the original, non-reversed table
        keep = len(table) - reverse_index - 1
        return table[keep]


# Register the custom identifier for MESA files
# This allows auto-detection of files named history*.data or profile*.data
try:
    from astropy.table import Table

    from . import connect

    connect.io_registry.register_identifier("ascii.mesa", Table, mesa_identify)
except Exception:
    # If registration fails (e.g., during import), just skip it
    pass
