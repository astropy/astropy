# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Classes to read AAS MRT table format.

Ref: https://journals.aas.org/mrt-standards

:Copyright: Smithsonian Astrophysical Observatory (2021)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu), \
         Suyog Garg (suyog7130@gmail.com)
"""

import re
import textwrap
from string import Template

import numpy as np

from . import cds, fixedwidth

__doctest_skip__ = ["*"]


MRT_TEMPLATE = [
    "$title",
    "$authors",
    "$table",
    "================================================================================",
    "$bytebybyte",
    "$notes",
    "--------------------------------------------------------------------------------",
]


class MrtSplitter(fixedwidth.FixedWidthSplitter):
    """
    Contains the join function to left align the MRT columns
    when writing to a file.
    """

    def join(self, vals, widths):
        vals = [val + " " * (width - len(val)) for val, width in zip(vals, widths)]
        return self.delimiter.join(vals)


class MrtHeader(cds.CdsHeader):
    _subfmt = "MRT"

    def _split_float_format(self, value):
        """
        Splits a Float string into different parts to find number
        of digits after decimal and check if the value is in Scientific
        notation.

        Parameters
        ----------
        value : str
            String containing the float value to split.

        Returns
        -------
        fmt: (int, int, int, bool, bool)
            List of values describing the Float string.
            (size, dec, ent, sign, exp)
            size, length of the given string.
            ent, number of digits before decimal point.
            dec, number of digits after decimal point.
            sign, whether or not given value signed.
            exp, is value in Scientific notation?
        """
        regfloat = re.compile(
            r"""(?P<sign> [+-]*)
                (?P<ent> [^eE.]+)
                (?P<deciPt> [.]*)
                (?P<decimals> [0-9]*)
                (?P<exp> [eE]*-*)[0-9]*""",
            re.VERBOSE,
        )
        mo = regfloat.match(value)

        if mo is None:
            raise Exception(f"{value} is not a float number")
        return (
            len(value),
            len(mo.group("ent")),
            len(mo.group("decimals")),
            mo.group("sign") != "",
            mo.group("exp") != "",
        )

    def _set_column_val_limits(self, col):
        """
        Sets the ``col.min`` and ``col.max`` column attributes,
        taking into account columns with Null values.
        """
        col.max = max(col)
        col.min = min(col)
        if col.max is np.ma.core.MaskedConstant:
            col.max = None
        if col.min is np.ma.core.MaskedConstant:
            col.min = None

    def _extract_metadata_lines(self, lines, meta_pattern, allow_gaps=True):
        """Extract metadata from a list of lines

        The first line of a metadata field is identified based on meta_pattern.
        Subsequent lines are appended to the field if they start with the minimum
        MRT indentation of 2 spaces and are non-empty.

        Parameters
        ----------
        lines : list[str]
            List of lines from an input MRT file.
        meta_pattern : str
            Metadata pattern with exactly two capture groups for keys and values.
            For example '^([A-Za-z]) = (.+)$' will match '{key} = {value}'
        allow_gaps: bool
            Allow gaps between metadata fields.
            Could be separators, other table segments, etc.

        Returns
        -------
        dict[str, list[str]]
            Dictionary mapping metadata keys to a list of lines from the input.
            For the first line of each field, only the second group of meta_pattern
            is included in the list.
        """
        pattern = re.compile(meta_pattern)
        if pattern.groups != 2:
            raise ValueError(
                "meta_pattern must have exactly two capture groups."
                f" Pattern {meta_pattern} has {pattern.groups} groups."
            )
        meta_lines = {}
        in_meta = False
        for line in lines:
            match = re.match(meta_pattern, line)
            if match:
                key, val = match.groups()
                meta_lines[key] = [val]
                in_meta = True
            elif line.startswith("  ") and not line.isspace() and in_meta:
                meta_lines[key].append(line)
            elif in_meta:
                if not allow_gaps:
                    raise ValueError("Unexpected gap in metadata fields")
                in_meta = False
        return meta_lines

    def _find_map_note(self, meta_key, meta_lines):
        map_note_pattern = "^([A-Za-z]) = (.+)$"
        if (
            meta_key.startswith("Note")
            and meta_lines[0].strip() == ""
            and re.match(map_note_pattern, meta_lines[1].strip())
        ):
            n_indent = len(meta_lines[1]) - len(meta_lines[1].lstrip())
            stripped_lines = [line[n_indent:] for line in meta_lines]
            try:
                mapped_note = self._extract_metadata_lines(
                    stripped_lines, map_note_pattern, allow_gaps=False
                )
            except ValueError as e:
                if str(e) == "Unexpected gap in metadata fields":
                    return None
                else:
                    raise e
            for map_key, map_val in mapped_note.items():
                mapped_note[map_key] = "\n".join([mline.strip() for mline in map_val])
            return mapped_note

    def update_meta(self, lines, meta):
        """
        Extract any table-level metadata, e.g. keywords, comments, column metadata, from
        the table ``lines`` and update the dict ``meta`` in place.
        For MRT tables, the extracted metadata includes article title, author list,
        table name, and notes.
        """
        # Not really necessary but avoids "catching" data lines erroneously
        head_lines = [line for line in lines if line not in self.data.data_lines]

        # Store all metadata fields (top part of header and notes) in a dictionary.
        # Keep lines separated in a list for each field.
        top_keys = ["Title", "Authors", "Table"]
        meta_key_patterns = top_keys + [r"Note \(\d+\)"]
        # Match empty strings on RHS, added as any other line to meta_dict for now
        meta_pattern = f"({'|'.join(meta_key_patterns)}):(.*)$"
        meta_dict = self._extract_metadata_lines(head_lines, meta_pattern)

        # Store mapping notes in dictionaries, store all others as strings.
        # Then split authors into a Python list.
        for meta_key, meta_lines in meta_dict.items():
            mapped_note = self._find_map_note(meta_key, meta_lines)
            if mapped_note is not None:
                meta_dict[meta_key] = mapped_note
            else:
                join_char = "\n" if meta_key.startswith("Note") else " "
                meta_dict[meta_key] = join_char.join(
                    [mline.strip() for mline in meta_lines]
                )
        meta_dict["Authors"] = meta_dict["Authors"].split(", ")

        top_meta = {}
        notes = []
        for key, val in meta_dict.items():
            if key in top_keys:
                top_meta[key] = val
            elif key.startswith("Note"):
                notes.append(val)
            else:
                raise ValueError(f"Unrecognized metadata key: {key}")
        meta.setdefault("table", {})["top"] = top_meta
        meta["table"]["notes"] = notes

    def column_float_formatter(self, col):
        """
        String formatter function for a column containing Float values.
        Checks if the values in the given column are in Scientific notation,
        by splitting the value string. It is assumed that the column either has
        float values or Scientific notation.

        A ``col.formatted_width`` attribute is added to the column. It is not added
        if such an attribute is already present, say when the ``formats`` argument
        is passed to the writer. A properly formatted format string is also added as
        the ``col.format`` attribute.

        Parameters
        ----------
        col : A ``Table.Column`` object.
        """
        # maxsize: maximum length of string containing the float value.
        # maxent: maximum number of digits places before decimal point.
        # maxdec: maximum number of digits places after decimal point.
        # maxprec: maximum precision of the column values, sum of maxent and maxdec.
        maxsize, maxprec, maxent, maxdec = 1, 0, 1, 0
        sign = False
        fformat = "F"

        # Find maximum sized value in the col
        for val in col.str_vals:
            # Skip null values
            if val is None or val == "":
                continue

            # Find format of the Float string
            fmt = self._split_float_format(val)
            # If value is in Scientific notation
            if fmt[4] is True:
                # if the previous column value was in normal Float format
                # set maxsize, maxprec and maxdec to default.
                if fformat == "F":
                    maxsize, maxprec, maxdec = 1, 0, 0
                # Designate the column to be in Scientific notation.
                fformat = "E"
            else:
                # Move to next column value if
                # current value is not in Scientific notation
                # but the column is designated as such because
                # one of the previous values was.
                if fformat == "E":
                    continue

            maxsize = max(maxsize, fmt[0])
            maxent = max(maxent, fmt[1])
            maxdec = max(maxdec, fmt[2])
            if fmt[3]:
                sign = True

            maxprec = max(maxprec, fmt[1] + fmt[2])

        if fformat == "E":
            # If ``formats`` not passed.
            if getattr(col, "formatted_width", None) is None:
                col.formatted_width = maxsize
                if sign:
                    col.formatted_width += 1
            # Number of digits after decimal is replaced by the precision
            # for values in Scientific notation, when writing that Format.
            col.fortran_format = fformat + str(col.formatted_width) + "." + str(maxprec)
            col.format = str(col.formatted_width) + "." + str(maxdec) + "e"
        else:
            lead = ""
            if (
                getattr(col, "formatted_width", None) is None
            ):  # If ``formats`` not passed.
                col.formatted_width = maxent + maxdec + 1
                if sign:
                    col.formatted_width += 1
            elif col.format.startswith("0"):
                # Keep leading zero, if already set in format - primarily for `seconds` columns
                # in coordinates; may need extra case if this is to be also supported with `sign`.
                lead = "0"
            col.fortran_format = fformat + str(col.formatted_width) + "." + str(maxdec)
            col.format = lead + col.fortran_format[1:] + "f"

    def write(self, lines):
        """
        Writes the Header of the MRT table, aka ReadMe, which
        also contains the Byte-By-Byte description of the table.
        """
        byte_by_byte = self.get_byte_by_byte()

        # Fill up the full ReadMe
        default_indent = "    "
        top_meta = self.table_meta.get("top", {})
        notes = self.table_meta.get("notes", {})
        from functools import partial

        wrap_meta = partial(
            textwrap.wrap,
            subsequent_indent=default_indent,
            width=cds.MAX_SIZE_README_LINE,
            break_long_words=False,
            break_on_hyphens=False,
        )
        # Indent title and authors
        top_out = {}
        top_out["Title"] = "\n".join(wrap_meta("Title: " + top_meta.get("Title", "")))
        top_out["Table"] = "\n".join(wrap_meta("Table: " + top_meta.get("Table", "")))
        # HACK: Use numbers to avoid breaking author names, then revert to spaces
        authors_hack = "Authors: " + "6".join(top_meta.get("Authors", [])).replace(
            " ", "7"
        ).replace("6", ", ")
        top_out["Authors"] = "\n".join(wrap_meta(authors_hack)).replace("7", " ")
        notes_str = []
        for i, note in enumerate(notes):
            note_prefix = f"Note ({i+1})"
            if isinstance(note, dict):
                map_notes = []
                for key, val in note.items():
                    map_str = f"{key} = {val}"
                    # Don't use default_indent: want len("k = "), so 4 spaces
                    map_str = textwrap.indent(
                        map_str,
                        "    ",
                        predicate=lambda x: not x.startswith(f"{key} = "),
                    )
                    map_notes.append(map_str)
                note = textwrap.indent("\n".join(map_notes), default_indent)
                notes_str.append(f"{note_prefix}:\n{note}")
            elif not isinstance(note, str):
                raise TypeError(
                    f"Unexpected type {type(note)} for note {note}. Expected str or dict."
                )
            else:
                note_str = textwrap.indent(
                    f"{note_prefix}: {note}",
                    default_indent,
                    predicate=lambda x: not x.startswith("Note ("),
                )
                notes_str.append(note_str)
        notes_str = "\n".join(notes_str)
        if notes_str == "":
            notes_str = "Notes:"

        rm_template = Template("\n".join(MRT_TEMPLATE))
        readme_filled = rm_template.substitute(
            {
                "bytebybyte": byte_by_byte,
                "title": top_out["Title"],
                "authors": top_out["Authors"],
                "table": top_out["Table"],
                "notes": notes_str,
            }
        )
        lines.append(readme_filled)


class MrtData(cds.CdsData):
    """MRT table data reader."""

    _subfmt = "MRT"
    splitter_class = MrtSplitter

    def write(self, lines):
        self.splitter.delimiter = " "
        fixedwidth.FixedWidthData.write(self, lines)


class Mrt(cds.Cds):
    """AAS MRT (Machine-Readable Table) format table.

    **Reading**
    ::

      >>> from astropy.io import ascii
      >>> table = ascii.read('data.mrt', format='mrt')

    **Writing**

    Use ``ascii.write(table, 'data.mrt', format='mrt')`` to  write tables to
    Machine Readable Table (MRT) format.

    Note that the metadata of the table, apart from units, column names and
    description, will not be written. These have to be filled in by hand later.

    See also: :ref:`cds_mrt_format`.

    Caveats:

    * The Units and Explanations are available in the column ``unit`` and
      ``description`` attributes, respectively.
    * The other metadata defined by this format is not available in the output table.
    """

    _format_name = "mrt"
    _io_registry_format_aliases = ["mrt"]
    _io_registry_can_write = True
    _description = "MRT format table"

    data_class = MrtData
    header_class = MrtHeader
