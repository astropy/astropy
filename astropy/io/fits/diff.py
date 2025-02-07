# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Facilities for diffing two FITS files.  Includes objects for diffing entire
FITS files, individual HDUs, FITS headers, or just FITS data.

Used to implement the fitsdiff program.
"""

import fnmatch
import glob
import io
import operator
import os
import os.path
import textwrap
from collections import defaultdict
from inspect import signature
from itertools import islice

import numpy as np

from astropy import __version__
from astropy.utils.diff import diff_values, report_diff_values, where_not_allclose
from astropy.utils.misc import NOT_OVERWRITING_MSG

from .card import BLANK_CARD, Card

# HDUList is used in one of the doctests
from .hdu.hdulist import HDUList, fitsopen  # pylint: disable=W0611
from .hdu.table import _TableLikeHDU
from .header import Header
from .util import path_like

__all__ = [
    "FITSDiff",
    "HDUDiff",
    "HeaderDiff",
    "ImageDataDiff",
    "RawDataDiff",
    "TableDataDiff",
]

# Column attributes of interest for comparison
_COL_ATTRS = [
    ("unit", "units"),
    ("null", "null values"),
    ("bscale", "bscales"),
    ("bzero", "bzeros"),
    ("disp", "display formats"),
    ("dim", "dimensions"),
]


def _get_differences(a, b):
    relative = abs(b - a) / abs(b)
    absolute = float(abs(b - a))
    return relative, absolute


class _BaseDiff:
    """
    Base class for all FITS diff objects.

    When instantiating a FITS diff object, the first two arguments are always
    the two objects to diff (two FITS files, two FITS headers, etc.).
    Instantiating a ``_BaseDiff`` also causes the diff itself to be executed.
    The returned ``_BaseDiff`` instance has a number of attribute that describe
    the results of the diff operation.

    The most basic attribute, present on all ``_BaseDiff`` instances, is
    ``.identical`` which is `True` if the two objects being compared are
    identical according to the diff method for objects of that type.
    """

    def __init__(self, a, b):
        """
        The ``_BaseDiff`` class does not implement a ``_diff`` method and
        should not be instantiated directly. Instead instantiate the
        appropriate subclass of ``_BaseDiff`` for the objects being compared
        (for example, use `HeaderDiff` to compare two `Header` objects.
        """
        self.a = a
        self.b = b

        # For internal use in report output
        self._fileobj = None
        self._indent = 0

        self._diff()

    def __bool__(self):
        """
        A ``_BaseDiff`` object acts as `True` in a boolean context if the two
        objects compared are different.  Otherwise it acts as `False`.
        """
        return not self.identical

    @classmethod
    def fromdiff(cls, other, a, b):
        """
        Returns a new Diff object of a specific subclass from an existing diff
        object, passing on the values for any arguments they share in common
        (such as ignore_keywords).

        For example::

            >>> from astropy.io import fits
            >>> hdul1, hdul2 = fits.HDUList(), fits.HDUList()
            >>> headera, headerb = fits.Header(), fits.Header()
            >>> fd = fits.FITSDiff(hdul1, hdul2, ignore_keywords=['*'])
            >>> hd = fits.HeaderDiff.fromdiff(fd, headera, headerb)
            >>> list(hd.ignore_keywords)
            ['*']
        """
        sig = signature(cls.__init__)
        # The first 3 arguments of any Diff initializer are self, a, and b.
        kwargs = {}
        for arg in list(sig.parameters.keys())[3:]:
            if hasattr(other, arg):
                kwargs[arg] = getattr(other, arg)

        return cls(a, b, **kwargs)

    @property
    def identical(self):
        """
        `True` if all the ``.diff_*`` attributes on this diff instance are
        empty, implying that no differences were found.

        Any subclass of ``_BaseDiff`` must have at least one ``.diff_*``
        attribute, which contains a non-empty value if and only if some
        difference was found between the two objects being compared.
        """
        return not any(
            getattr(self, attr) for attr in self.__dict__ if attr.startswith("diff_")
        )

    def report(self, fileobj=None, indent=0, overwrite=False):
        """
        Generates a text report on the differences (if any) between two
        objects, and either returns it as a string or writes it to a file-like
        object.

        Parameters
        ----------
        fileobj : file-like, string, or None, optional
            If `None`, this method returns the report as a string. Otherwise it
            returns `None` and writes the report to the given file-like object
            (which must have a ``.write()`` method at a minimum), or to a new
            file at the path specified.

        indent : int
            The number of 4 space tabs to indent the report.

        overwrite : bool, optional
            If ``True``, overwrite the output file if it exists. Raises an
            ``OSError`` if ``False`` and the output file exists. Default is
            ``False``.

        Returns
        -------
        report : str or None
        """
        return_string = False
        filepath = None

        if isinstance(fileobj, path_like):
            fileobj = os.path.expanduser(fileobj)
            if os.path.exists(fileobj) and not overwrite:
                raise OSError(NOT_OVERWRITING_MSG.format(fileobj))
            else:
                filepath = fileobj
                fileobj = open(filepath, "w")
        elif fileobj is None:
            fileobj = io.StringIO()
            return_string = True

        self._fileobj = fileobj
        self._indent = indent  # This is used internally by _writeln

        try:
            self._report()
        finally:
            if filepath:
                fileobj.close()

        if return_string:
            return fileobj.getvalue()

    def _writeln(self, text):
        self._fileobj.write(textwrap.indent(text, self._indent * "  ") + "\n")

    def _diff(self):
        raise NotImplementedError

    def _report(self):
        raise NotImplementedError


class FITSDiff(_BaseDiff):
    """Diff two FITS files by filename, or two `HDUList` objects.

    `FITSDiff` objects have the following diff attributes:

    - ``diff_hdu_count``: If the FITS files being compared have different
      numbers of HDUs, this contains a 2-tuple of the number of HDUs in each
      file.

    - ``diff_hdus``: If any HDUs with the same index are different, this
      contains a list of 2-tuples of the HDU index and the `HDUDiff` object
      representing the differences between the two HDUs.
    """

    def __init__(
        self,
        a,
        b,
        ignore_hdus=[],
        ignore_keywords=[],
        ignore_comments=[],
        ignore_fields=[],
        numdiffs=10,
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
    ):
        """
        Parameters
        ----------
        a : str or `HDUList`
            The filename of a FITS file on disk, or an `HDUList` object.

        b : str or `HDUList`
            The filename of a FITS file on disk, or an `HDUList` object to
            compare to the first file.

        ignore_hdus : sequence, optional
            HDU names to ignore when comparing two FITS files or HDU lists; the
            presence of these HDUs and their contents are ignored.  Wildcard
            strings may also be included in the list.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers; the presence
            of these keywords and their values are ignored.  Wildcard strings
            may also be included in the list.

        ignore_comments : sequence, optional
            A list of header keywords whose comments should be ignored in the
            comparison.  May contain wildcard strings as with ignore_keywords.

        ignore_fields : sequence, optional
            The (case-insensitive) names of any table columns to ignore if any
            table data is to be compared.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored
            (default: True).

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain
            whitespace (default: True).
        """
        if isinstance(a, (str, os.PathLike)):
            try:
                a = fitsopen(a)
            except Exception as exc:
                raise OSError(f"error opening file a ({a})") from exc
            close_a = True
        else:
            close_a = False

        if isinstance(b, (str, os.PathLike)):
            try:
                b = fitsopen(b)
            except Exception as exc:
                raise OSError(f"error opening file b ({b})") from exc
            close_b = True
        else:
            close_b = False

        # Normalize keywords/fields to ignore to upper case
        self.ignore_hdus = {k.upper() for k in ignore_hdus}
        self.ignore_keywords = {k.upper() for k in ignore_keywords}
        self.ignore_comments = {k.upper() for k in ignore_comments}
        self.ignore_fields = {k.upper() for k in ignore_fields}

        self.numdiffs = numdiffs
        self.rtol = rtol
        self.atol = atol

        self.ignore_blanks = ignore_blanks
        self.ignore_blank_cards = ignore_blank_cards

        # Some hdu names may be pattern wildcards.  Find them.
        self.ignore_hdu_patterns = set()
        for name in list(self.ignore_hdus):
            if name != "*" and glob.has_magic(name):
                self.ignore_hdus.remove(name)
                self.ignore_hdu_patterns.add(name)

        self.diff_hdu_count = ()
        self.diff_hdus = []

        try:
            super().__init__(a, b)
        finally:
            if close_a:
                a.close()
            if close_b:
                b.close()

    def _diff(self):
        if len(self.a) != len(self.b):
            self.diff_hdu_count = (len(self.a), len(self.b))

        # Record filenames for use later in _report
        self.filenamea = self.a.filename()
        if not self.filenamea:
            self.filenamea = f"<{self.a.__class__.__name__} object at {id(self.a):#x}>"

        self.filenameb = self.b.filename()
        if not self.filenameb:
            self.filenameb = f"<{self.b.__class__.__name__} object at {id(self.b):#x}>"

        if self.ignore_hdus:
            self.a = HDUList([h for h in self.a if h.name not in self.ignore_hdus])
            self.b = HDUList([h for h in self.b if h.name not in self.ignore_hdus])
        if self.ignore_hdu_patterns:
            a_names = [hdu.name for hdu in self.a]
            b_names = [hdu.name for hdu in self.b]
            for pattern in self.ignore_hdu_patterns:
                a_ignored = fnmatch.filter(a_names, pattern)
                self.a = HDUList([h for h in self.a if h.name not in a_ignored])
                b_ignored = fnmatch.filter(b_names, pattern)
                self.b = HDUList([h for h in self.b if h.name not in b_ignored])

        # For now, just compare the extensions one by one in order.
        # Might allow some more sophisticated types of diffing later.

        # TODO: Somehow or another simplify the passing around of diff
        # options--this will become important as the number of options grows
        for idx in range(min(len(self.a), len(self.b))):
            hdu_diff = HDUDiff.fromdiff(self, self.a[idx], self.b[idx])

            if not hdu_diff.identical:
                if (
                    self.a[idx].name == self.b[idx].name
                    and self.a[idx].ver == self.b[idx].ver
                ):
                    self.diff_hdus.append(
                        (idx, hdu_diff, self.a[idx].name, self.a[idx].ver)
                    )
                else:
                    self.diff_hdus.append((idx, hdu_diff, "", self.a[idx].ver))

    def _report(self):
        wrapper = textwrap.TextWrapper(initial_indent="  ", subsequent_indent="  ")

        self._fileobj.write("\n")
        self._writeln(f" fitsdiff: {__version__}")
        self._writeln(f" a: {self.filenamea}\n b: {self.filenameb}")

        if self.ignore_hdus:
            ignore_hdus = " ".join(sorted(self.ignore_hdus))
            self._writeln(" HDU(s) not to be compared:\n" + wrapper.fill(ignore_hdus))

        if self.ignore_hdu_patterns:
            ignore_hdu_patterns = " ".join(sorted(self.ignore_hdu_patterns))
            self._writeln(
                " HDU(s) not to be compared:\n" + wrapper.fill(ignore_hdu_patterns)
            )

        if self.ignore_keywords:
            ignore_keywords = " ".join(sorted(self.ignore_keywords))
            self._writeln(
                " Keyword(s) not to be compared:\n" + wrapper.fill(ignore_keywords)
            )

        if self.ignore_comments:
            ignore_comments = " ".join(sorted(self.ignore_comments))
            self._writeln(
                " Keyword(s) whose comments are not to be compared:\n"
                + wrapper.fill(ignore_comments)
            )

        if self.ignore_fields:
            ignore_fields = " ".join(sorted(self.ignore_fields))
            self._writeln(
                " Table column(s) not to be compared:\n" + wrapper.fill(ignore_fields)
            )

        self._writeln(
            f" Maximum number of different data values to be reported: {self.numdiffs}"
        )
        self._writeln(
            f" Relative tolerance: {self.rtol}, Absolute tolerance: {self.atol}"
        )

        if self.diff_hdu_count:
            self._fileobj.write("\n")
            self._writeln("Files contain different numbers of HDUs:")
            self._writeln(f" a: {self.diff_hdu_count[0]}")
            self._writeln(f" b: {self.diff_hdu_count[1]}")

            if not self.diff_hdus:
                self._writeln("No differences found between common HDUs.")
                return
        elif not self.diff_hdus:
            self._fileobj.write("\n")
            self._writeln("No differences found.")
            return

        for idx, hdu_diff, extname, extver in self.diff_hdus:
            # print out the extension heading
            if idx == 0:
                self._fileobj.write("\n")
                self._writeln("Primary HDU:")
            else:
                self._fileobj.write("\n")
                if extname:
                    self._writeln(f"Extension HDU {idx} ({extname}, {extver}):")
                else:
                    self._writeln(f"Extension HDU {idx}:")
            hdu_diff.report(self._fileobj, indent=self._indent + 1)


class HDUDiff(_BaseDiff):
    """
    Diff two HDU objects, including their headers and their data (but only if
    both HDUs contain the same type of data (image, table, or unknown).

    `HDUDiff` objects have the following diff attributes:

    - ``diff_extnames``: If the two HDUs have different EXTNAME values, this
      contains a 2-tuple of the different extension names.

    - ``diff_extvers``: If the two HDUS have different EXTVER values, this
      contains a 2-tuple of the different extension versions.

    - ``diff_extlevels``: If the two HDUs have different EXTLEVEL values, this
      contains a 2-tuple of the different extension levels.

    - ``diff_extension_types``: If the two HDUs have different XTENSION values,
      this contains a 2-tuple of the different extension types.

    - ``diff_headers``: Contains a `HeaderDiff` object for the headers of the
      two HDUs. This will always contain an object--it may be determined
      whether the headers are different through ``diff_headers.identical``.

    - ``diff_data``: Contains either a `ImageDataDiff`, `TableDataDiff`, or
      `RawDataDiff` as appropriate for the data in the HDUs, and only if the
      two HDUs have non-empty data of the same type (`RawDataDiff` is used for
      HDUs containing non-empty data of an indeterminate type).
    """

    def __init__(
        self,
        a,
        b,
        ignore_keywords=[],
        ignore_comments=[],
        ignore_fields=[],
        numdiffs=10,
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
    ):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers; the presence
            of these keywords and their values are ignored.  Wildcard strings
            may also be included in the list.

        ignore_comments : sequence, optional
            A list of header keywords whose comments should be ignored in the
            comparison.  May contain wildcard strings as with ignore_keywords.

        ignore_fields : sequence, optional
            The (case-insensitive) names of any table columns to ignore if any
            table data is to be compared.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored
            (default: True).

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain
            whitespace (default: True).
        """
        self.ignore_keywords = {k.upper() for k in ignore_keywords}
        self.ignore_comments = {k.upper() for k in ignore_comments}
        self.ignore_fields = {k.upper() for k in ignore_fields}

        self.rtol = rtol
        self.atol = atol

        self.numdiffs = numdiffs
        self.ignore_blanks = ignore_blanks
        self.ignore_blank_cards = ignore_blank_cards

        self.diff_extnames = ()
        self.diff_extvers = ()
        self.diff_extlevels = ()
        self.diff_extension_types = ()
        self.diff_headers = None
        self.diff_data = None

        super().__init__(a, b)

    def _diff(self):
        if self.a.name != self.b.name:
            self.diff_extnames = (self.a.name, self.b.name)

        if self.a.ver != self.b.ver:
            self.diff_extvers = (self.a.ver, self.b.ver)

        if self.a.level != self.b.level:
            self.diff_extlevels = (self.a.level, self.b.level)

        if self.a.header.get("XTENSION") != self.b.header.get("XTENSION"):
            self.diff_extension_types = (
                self.a.header.get("XTENSION"),
                self.b.header.get("XTENSION"),
            )

        self.diff_headers = HeaderDiff.fromdiff(
            self, self.a.header.copy(), self.b.header.copy()
        )

        if self.a.data is None or self.b.data is None:
            # TODO: Perhaps have some means of marking this case
            pass
        elif self.a.is_image and self.b.is_image:
            self.diff_data = ImageDataDiff.fromdiff(self, self.a.data, self.b.data)
            # Clean up references to (possibly) memmapped arrays so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None
        elif isinstance(self.a, _TableLikeHDU) and isinstance(self.b, _TableLikeHDU):
            # TODO: Replace this if/when _BaseHDU grows a .is_table property
            self.diff_data = TableDataDiff.fromdiff(self, self.a.data, self.b.data)
            # Clean up references to (possibly) memmapped arrays so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None
        elif not self.diff_extension_types:
            # Don't diff the data for unequal extension types that are not
            # recognized image or table types
            self.diff_data = RawDataDiff.fromdiff(self, self.a.data, self.b.data)
            # Clean up references to (possibly) memmapped arrays so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None

    def _report(self):
        if self.identical:
            self._writeln(" No differences found.")
        if self.diff_extension_types:
            self._writeln(
                " Extension types differ:\n  a: {}\n  b: {}".format(
                    *self.diff_extension_types
                )
            )
        if self.diff_extnames:
            self._writeln(
                " Extension names differ:\n  a: {}\n  b: {}".format(*self.diff_extnames)
            )
        if self.diff_extvers:
            self._writeln(
                " Extension versions differ:\n  a: {}\n  b: {}".format(
                    *self.diff_extvers
                )
            )

        if self.diff_extlevels:
            self._writeln(
                " Extension levels differ:\n  a: {}\n  b: {}".format(
                    *self.diff_extlevels
                )
            )

        if not self.diff_headers.identical:
            self._fileobj.write("\n")
            self._writeln(" Headers contain differences:")
            self.diff_headers.report(self._fileobj, indent=self._indent + 1)

        if self.diff_data is not None and not self.diff_data.identical:
            self._fileobj.write("\n")
            self._writeln(" Data contains differences:")
            self.diff_data.report(self._fileobj, indent=self._indent + 1)


class HeaderDiff(_BaseDiff):
    """
    Diff two `Header` objects.

    `HeaderDiff` objects have the following diff attributes:

    - ``diff_keyword_count``: If the two headers contain a different number of
      keywords, this contains a 2-tuple of the keyword count for each header.

    - ``diff_keywords``: If either header contains one or more keywords that
      don't appear at all in the other header, this contains a 2-tuple
      consisting of a list of the keywords only appearing in header a, and a
      list of the keywords only appearing in header b.

    - ``diff_duplicate_keywords``: If a keyword appears in both headers at
      least once, but contains a different number of duplicates (for example, a
      different number of HISTORY cards in each header), an item is added to
      this dict with the keyword as the key, and a 2-tuple of the different
      counts of that keyword as the value.  For example::

          {'HISTORY': (20, 19)}

      means that header a contains 20 HISTORY cards, while header b contains
      only 19 HISTORY cards.

    - ``diff_keyword_values``: If any of the common keyword between the two
      headers have different values, they appear in this dict.  It has a
      structure similar to ``diff_duplicate_keywords``, with the keyword as the
      key, and a 2-tuple of the different values as the value.  For example::

          {'NAXIS': (2, 3)}

      means that the NAXIS keyword has a value of 2 in header a, and a value of
      3 in header b.  This excludes any keywords matched by the
      ``ignore_keywords`` list.

    - ``diff_keyword_comments``: Like ``diff_keyword_values``, but contains
      differences between keyword comments.

    `HeaderDiff` objects also have a ``common_keywords`` attribute that lists
    all keywords that appear in both headers.
    """

    def __init__(
        self,
        a,
        b,
        ignore_keywords=[],
        ignore_comments=[],
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
    ):
        """
        Parameters
        ----------
        a : `~astropy.io.fits.Header` or string or bytes
            A header.

        b : `~astropy.io.fits.Header` or string or bytes
            A header to compare to the first header.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers; the presence
            of these keywords and their values are ignored.  Wildcard strings
            may also be included in the list.

        ignore_comments : sequence, optional
            A list of header keywords whose comments should be ignored in the
            comparison.  May contain wildcard strings as with ignore_keywords.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored
            (default: True).

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain
            whitespace (default: True).
        """
        self.ignore_keywords = {k.upper() for k in ignore_keywords}
        self.ignore_comments = {k.upper() for k in ignore_comments}

        self.rtol = rtol
        self.atol = atol

        self.ignore_blanks = ignore_blanks
        self.ignore_blank_cards = ignore_blank_cards

        self.ignore_keyword_patterns = set()
        self.ignore_comment_patterns = set()
        for keyword in list(self.ignore_keywords):
            keyword = keyword.upper()
            if keyword != "*" and glob.has_magic(keyword):
                self.ignore_keywords.remove(keyword)
                self.ignore_keyword_patterns.add(keyword)
        for keyword in list(self.ignore_comments):
            keyword = keyword.upper()
            if keyword != "*" and glob.has_magic(keyword):
                self.ignore_comments.remove(keyword)
                self.ignore_comment_patterns.add(keyword)

        # Keywords appearing in each header
        self.common_keywords = []

        # Set to the number of keywords in each header if the counts differ
        self.diff_keyword_count = ()

        # Set if the keywords common to each header (excluding ignore_keywords)
        # appear in different positions within the header
        # TODO: Implement this
        self.diff_keyword_positions = ()

        # Keywords unique to each header (excluding keywords in
        # ignore_keywords)
        self.diff_keywords = ()

        # Keywords that have different numbers of duplicates in each header
        # (excluding keywords in ignore_keywords)
        self.diff_duplicate_keywords = {}

        # Keywords common to each header but having different values (excluding
        # keywords in ignore_keywords)
        self.diff_keyword_values = defaultdict(list)

        # Keywords common to each header but having different comments
        # (excluding keywords in ignore_keywords or in ignore_comments)
        self.diff_keyword_comments = defaultdict(list)

        if isinstance(a, str):
            a = Header.fromstring(a)
        if isinstance(b, str):
            b = Header.fromstring(b)

        if not (isinstance(a, Header) and isinstance(b, Header)):
            raise TypeError(
                "HeaderDiff can only diff astropy.io.fits.Header "
                "objects or strings containing FITS headers."
            )

        super().__init__(a, b)

    # TODO: This doesn't pay much attention to the *order* of the keywords,
    # except in the case of duplicate keywords.  The order should be checked
    # too, or at least it should be an option.
    def _diff(self):
        if self.ignore_blank_cards:
            cardsa = [c for c in self.a.cards if str(c) != BLANK_CARD]
            cardsb = [c for c in self.b.cards if str(c) != BLANK_CARD]
        else:
            cardsa = list(self.a.cards)
            cardsb = list(self.b.cards)

        # build dictionaries of keyword values and comments
        def get_header_values_comments(cards):
            values = {}
            comments = {}
            for card in cards:
                value = card.value
                if self.ignore_blanks and isinstance(value, str):
                    value = value.rstrip()
                values.setdefault(card.keyword, []).append(value)
                comments.setdefault(card.keyword, []).append(card.comment)
            return values, comments

        valuesa, commentsa = get_header_values_comments(cardsa)
        valuesb, commentsb = get_header_values_comments(cardsb)

        keywordsa = set(valuesa)
        keywordsb = set(valuesb)

        self.common_keywords = sorted(keywordsa.intersection(keywordsb))
        if len(cardsa) != len(cardsb):
            self.diff_keyword_count = (len(cardsa), len(cardsb))

        # Any other diff attributes should exclude ignored keywords
        keywordsa = keywordsa.difference(self.ignore_keywords)
        keywordsb = keywordsb.difference(self.ignore_keywords)
        if self.ignore_keyword_patterns:
            for pattern in self.ignore_keyword_patterns:
                keywordsa = keywordsa.difference(fnmatch.filter(keywordsa, pattern))
                keywordsb = keywordsb.difference(fnmatch.filter(keywordsb, pattern))

        if "*" in self.ignore_keywords:
            # Any other differences between keywords are to be ignored
            return

        left_only_keywords = sorted(keywordsa.difference(keywordsb))
        right_only_keywords = sorted(keywordsb.difference(keywordsa))

        if left_only_keywords or right_only_keywords:
            self.diff_keywords = (left_only_keywords, right_only_keywords)

        # Compare count of each common keyword
        for keyword in self.common_keywords:
            if keyword in self.ignore_keywords:
                continue
            if self.ignore_keyword_patterns:
                skip = False
                for pattern in self.ignore_keyword_patterns:
                    if fnmatch.fnmatch(keyword, pattern):
                        skip = True
                        break
                if skip:
                    continue

            counta = len(valuesa[keyword])
            countb = len(valuesb[keyword])
            if counta != countb:
                self.diff_duplicate_keywords[keyword] = (counta, countb)

            # Compare keywords' values and comments
            for a, b in zip(valuesa[keyword], valuesb[keyword]):
                if diff_values(a, b, rtol=self.rtol, atol=self.atol):
                    self.diff_keyword_values[keyword].append((a, b))
                else:
                    # If there are duplicate keywords we need to be able to
                    # index each duplicate; if the values of a duplicate
                    # are identical use None here
                    self.diff_keyword_values[keyword].append(None)

            if not any(self.diff_keyword_values[keyword]):
                # No differences found; delete the array of Nones
                del self.diff_keyword_values[keyword]

            if "*" in self.ignore_comments or keyword in self.ignore_comments:
                continue
            if self.ignore_comment_patterns:
                skip = False
                for pattern in self.ignore_comment_patterns:
                    if fnmatch.fnmatch(keyword, pattern):
                        skip = True
                        break
                if skip:
                    continue

            for a, b in zip(commentsa[keyword], commentsb[keyword]):
                if diff_values(a, b):
                    self.diff_keyword_comments[keyword].append((a, b))
                else:
                    self.diff_keyword_comments[keyword].append(None)

            if not any(self.diff_keyword_comments[keyword]):
                del self.diff_keyword_comments[keyword]

    def _report(self):
        if self.diff_keyword_count:
            self._writeln(" Headers have different number of cards:")
            self._writeln(f"  a: {self.diff_keyword_count[0]}")
            self._writeln(f"  b: {self.diff_keyword_count[1]}")
        if self.diff_keywords:
            for keyword in self.diff_keywords[0]:
                if keyword in Card._commentary_keywords:
                    val = self.a[keyword][0]
                else:
                    val = self.a[keyword]
                self._writeln(f" Extra keyword {keyword!r:8} in a: {val!r}")
            for keyword in self.diff_keywords[1]:
                if keyword in Card._commentary_keywords:
                    val = self.b[keyword][0]
                else:
                    val = self.b[keyword]
                self._writeln(f" Extra keyword {keyword!r:8} in b: {val!r}")

        if self.diff_duplicate_keywords:
            for keyword, count in sorted(self.diff_duplicate_keywords.items()):
                self._writeln(f" Inconsistent duplicates of keyword {keyword!r:8}:")
                self._writeln(
                    "  Occurs {} time(s) in a, {} times in (b)".format(*count)
                )

        if self.diff_keyword_values or self.diff_keyword_comments:
            for keyword in self.common_keywords:
                report_diff_keyword_attr(
                    self._fileobj,
                    "values",
                    self.diff_keyword_values,
                    keyword,
                    ind=self._indent,
                )
                report_diff_keyword_attr(
                    self._fileobj,
                    "comments",
                    self.diff_keyword_comments,
                    keyword,
                    ind=self._indent,
                )


# TODO: It might be good if there was also a threshold option for percentage of
# different pixels: For example ignore if only 1% of the pixels are different
# within some threshold.  There are lots of possibilities here, but hold off
# for now until specific cases come up.


class ImageDataDiff(_BaseDiff):
    """
    Diff two image data arrays (really any array from a PRIMARY HDU or an IMAGE
    extension HDU, though the data unit is assumed to be "pixels").

    `ImageDataDiff` objects have the following diff attributes:

    - ``diff_dimensions``: If the two arrays contain either a different number
      of dimensions or different sizes in any dimension, this contains a
      2-tuple of the shapes of each array.  Currently no further comparison is
      performed on images that don't have the exact same dimensions.

    - ``diff_pixels``: If the two images contain any different pixels, this
      contains a list of 2-tuples of the array index where the difference was
      found, and another 2-tuple containing the different values.  For example,
      if the pixel at (0, 0) contains different values this would look like::

          [(0, 0), (1.1, 2.2)]

      where 1.1 and 2.2 are the values of that pixel in each array.  This
      array only contains up to ``self.numdiffs`` differences, for storage
      efficiency.

    - ``diff_total``: The total number of different pixels found between the
      arrays.  Although ``diff_pixels`` does not necessarily contain all the
      different pixel values, this can be used to get a count of the total
      number of differences found.

    - ``diff_ratio``: Contains the ratio of ``diff_total`` to the total number
      of pixels in the arrays.
    """

    def __init__(self, a, b, numdiffs=10, rtol=0.0, atol=0.0):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0
        """
        self.numdiffs = numdiffs
        self.rtol = rtol
        self.atol = atol

        self.diff_dimensions = ()
        self.diff_pixels = []
        self.diff_ratio = 0

        # self.diff_pixels only holds up to numdiffs differing pixels, but this
        # self.diff_total stores the total count of differences between
        # the images, but not the different values
        self.diff_total = 0

        super().__init__(a, b)

    def _diff(self):
        if self.a.shape != self.b.shape:
            self.diff_dimensions = (self.a.shape, self.b.shape)
            # Don't do any further comparison if the dimensions differ
            # TODO: Perhaps we could, however, diff just the intersection
            # between the two images
            return

        # Find the indices where the values are not equal
        # If neither a nor b are floating point (or complex), ignore rtol and
        # atol
        if not (
            np.issubdtype(self.a.dtype, np.inexact)
            or np.issubdtype(self.b.dtype, np.inexact)
        ):
            rtol = 0
            atol = 0
        else:
            rtol = self.rtol
            atol = self.atol

        diffs = where_not_allclose(self.a, self.b, atol=atol, rtol=rtol)

        self.diff_total = len(diffs[0])

        if self.diff_total == 0:
            # Then we're done
            return

        if self.numdiffs < 0:
            numdiffs = self.diff_total
        else:
            numdiffs = self.numdiffs

        self.diff_pixels = [
            (idx, (self.a[idx], self.b[idx]))
            for idx in islice(zip(*diffs), 0, numdiffs)
        ]
        self.diff_ratio = float(self.diff_total) / float(len(self.a.flat))

    def _report(self):
        if self.diff_dimensions:
            dimsa = " x ".join(str(d) for d in reversed(self.diff_dimensions[0]))
            dimsb = " x ".join(str(d) for d in reversed(self.diff_dimensions[1]))
            self._writeln(" Data dimensions differ:")
            self._writeln(f"  a: {dimsa}")
            self._writeln(f"  b: {dimsb}")
            # For now we don't do any further comparison if the dimensions
            # differ; though in the future it might be nice to be able to
            # compare at least where the images intersect
            self._writeln(" No further data comparison performed.")
            return

        if not self.diff_pixels:
            return

        max_relative = 0
        max_absolute = 0

        for index, values in self.diff_pixels:
            # Convert to int to avoid np.int64 in list repr.
            index = [int(x + 1) for x in reversed(index)]
            self._writeln(f" Data differs at {index}:")
            report_diff_values(
                values[0],
                values[1],
                fileobj=self._fileobj,
                indent_width=self._indent + 1,
                rtol=self.rtol,
                atol=self.atol,
            )
            rdiff, adiff = _get_differences(values[0], values[1])
            max_relative = max(max_relative, rdiff)
            max_absolute = max(max_absolute, adiff)

        if self.diff_total > self.numdiffs:
            self._writeln(" ...")
        self._writeln(
            f" {self.diff_total} different pixels found "
            f"({self.diff_ratio:.2%} different)."
        )
        self._writeln(f" Maximum relative difference: {max_relative}")
        self._writeln(f" Maximum absolute difference: {max_absolute}")


class RawDataDiff(ImageDataDiff):
    """
    `RawDataDiff` is just a special case of `ImageDataDiff` where the images
    are one-dimensional, and the data is treated as a 1-dimensional array of
    bytes instead of pixel values.  This is used to compare the data of two
    non-standard extension HDUs that were not recognized as containing image or
    table data.

    `ImageDataDiff` objects have the following diff attributes:

    - ``diff_dimensions``: Same as the ``diff_dimensions`` attribute of
      `ImageDataDiff` objects. Though the "dimension" of each array is just an
      integer representing the number of bytes in the data.

    - ``diff_bytes``: Like the ``diff_pixels`` attribute of `ImageDataDiff`
      objects, but renamed to reflect the minor semantic difference that these
      are raw bytes and not pixel values.  Also the indices are integers
      instead of tuples.

    - ``diff_total`` and ``diff_ratio``: Same as `ImageDataDiff`.
    """

    def __init__(self, a, b, numdiffs=10):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).
        """
        self.diff_dimensions = ()
        self.diff_bytes = []

        super().__init__(a, b, numdiffs=numdiffs)

    def _diff(self):
        super()._diff()
        if self.diff_dimensions:
            self.diff_dimensions = (
                self.diff_dimensions[0][0],
                self.diff_dimensions[1][0],
            )

        self.diff_bytes = [(x[0], y) for x, y in self.diff_pixels]
        del self.diff_pixels

    def _report(self):
        if self.diff_dimensions:
            self._writeln(" Data sizes differ:")
            self._writeln(f"  a: {self.diff_dimensions[0]} bytes")
            self._writeln(f"  b: {self.diff_dimensions[1]} bytes")
            # For now we don't do any further comparison if the dimensions
            # differ; though in the future it might be nice to be able to
            # compare at least where the images intersect
            self._writeln(" No further data comparison performed.")
            return

        if not self.diff_bytes:
            return

        for index, values in self.diff_bytes:
            self._writeln(f" Data differs at byte {index}:")
            report_diff_values(
                values[0],
                values[1],
                fileobj=self._fileobj,
                indent_width=self._indent + 1,
                rtol=self.rtol,
                atol=self.atol,
            )

        self._writeln(" ...")
        self._writeln(
            f" {self.diff_total} different bytes found "
            f"({self.diff_ratio:.2%} different)."
        )


class TableDataDiff(_BaseDiff):
    """
    Diff two table data arrays. It doesn't matter whether the data originally
    came from a binary or ASCII table--the data should be passed in as a
    recarray.

    `TableDataDiff` objects have the following diff attributes:

    - ``diff_column_count``: If the tables being compared have different
      numbers of columns, this contains a 2-tuple of the column count in each
      table.  Even if the tables have different column counts, an attempt is
      still made to compare any columns they have in common.

    - ``diff_columns``: If either table contains columns unique to that table,
      either in name or format, this contains a 2-tuple of lists. The first
      element is a list of columns (these are full `Column` objects) that
      appear only in table a.  The second element is a list of tables that
      appear only in table b.  This only lists columns with different column
      definitions, and has nothing to do with the data in those columns.

    - ``diff_column_names``: This is like ``diff_columns``, but lists only the
      names of columns unique to either table, rather than the full `Column`
      objects.

    - ``diff_column_attributes``: Lists columns that are in both tables but
      have different secondary attributes, such as TUNIT or TDISP.  The format
      is a list of 2-tuples: The first a tuple of the column name and the
      attribute, the second a tuple of the different values.

    - ``diff_values``: `TableDataDiff` compares the data in each table on a
      column-by-column basis.  If any different data is found, it is added to
      this list.  The format of this list is similar to the ``diff_pixels``
      attribute on `ImageDataDiff` objects, though the "index" consists of a
      (column_name, row) tuple.  For example::

          [('TARGET', 0), ('NGC1001', 'NGC1002')]

      shows that the tables contain different values in the 0-th row of the
      'TARGET' column.

    - ``diff_total`` and ``diff_ratio``: Same as `ImageDataDiff`.

    `TableDataDiff` objects also have a ``common_columns`` attribute that lists
    the `Column` objects for columns that are identical in both tables, and a
    ``common_column_names`` attribute which contains a set of the names of
    those columns.
    """

    def __init__(self, a, b, ignore_fields=[], numdiffs=10, rtol=0.0, atol=0.0):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        ignore_fields : sequence, optional
            The (case-insensitive) names of any table columns to ignore if any
            table data is to be compared.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0
        """
        self.ignore_fields = set(ignore_fields)
        self.numdiffs = numdiffs
        self.rtol = rtol
        self.atol = atol

        self.common_columns = []
        self.common_column_names = set()

        # self.diff_columns contains columns with different column definitions,
        # but not different column data. Column data is only compared in
        # columns that have the same definitions
        self.diff_rows = ()
        self.diff_column_count = ()
        self.diff_columns = ()

        # If two columns have the same name+format, but other attributes are
        # different (such as TUNIT or such) they are listed here
        self.diff_column_attributes = []

        # Like self.diff_columns, but just contains a list of the column names
        # unique to each table, and in the order they appear in the tables
        self.diff_column_names = ()
        self.diff_values = []

        self.diff_ratio = 0
        self.diff_total = 0

        super().__init__(a, b)

    def _diff(self):
        # Much of the code for comparing columns is similar to the code for
        # comparing headers--consider refactoring
        colsa = self.a.columns
        colsb = self.b.columns

        if len(colsa) != len(colsb):
            self.diff_column_count = (len(colsa), len(colsb))

        # Even if the number of columns are unequal, we still do comparison of
        # any common columns
        colsa = {c.name.lower(): c for c in colsa}
        colsb = {c.name.lower(): c for c in colsb}

        if "*" in self.ignore_fields:
            # If all columns are to be ignored, ignore any further differences
            # between the columns
            return

        # Keep the user's original ignore_fields list for reporting purposes,
        # but internally use a case-insensitive version
        ignore_fields = {f.lower() for f in self.ignore_fields}

        # It might be nice if there were a cleaner way to do this, but for now
        # it'll do
        for fieldname in ignore_fields:
            fieldname = fieldname.lower()
            if fieldname in colsa:
                del colsa[fieldname]
            if fieldname in colsb:
                del colsb[fieldname]

        colsa_set = set(colsa.values())
        colsb_set = set(colsb.values())
        self.common_columns = sorted(
            colsa_set.intersection(colsb_set), key=operator.attrgetter("name")
        )

        self.common_column_names = {col.name.lower() for col in self.common_columns}

        left_only_columns = {
            col.name.lower(): col for col in colsa_set.difference(colsb_set)
        }
        right_only_columns = {
            col.name.lower(): col for col in colsb_set.difference(colsa_set)
        }

        if left_only_columns or right_only_columns:
            self.diff_columns = (left_only_columns, right_only_columns)
            self.diff_column_names = ([], [])

        if left_only_columns:
            for col in self.a.columns:
                if col.name.lower() in left_only_columns:
                    self.diff_column_names[0].append(col.name)

        if right_only_columns:
            for col in self.b.columns:
                if col.name.lower() in right_only_columns:
                    self.diff_column_names[1].append(col.name)

        # If the tables have a different number of rows, we don't compare the
        # columns right now.
        # TODO: It might be nice to optionally compare the first n rows where n
        # is the minimum of the row counts between the two tables.
        if len(self.a) != len(self.b):
            self.diff_rows = (len(self.a), len(self.b))
            return

        # If the tables contain no rows there's no data to compare, so we're
        # done at this point. (See ticket #178)
        if len(self.a) == len(self.b) == 0:
            return

        # Like in the old fitsdiff, compare tables on a column by column basis
        # The difficulty here is that, while FITS column names are meant to be
        # case-insensitive, Astropy still allows, for the sake of flexibility,
        # two columns with the same name but different case.  When columns are
        # accessed in FITS tables, a case-sensitive is tried first, and failing
        # that a case-insensitive match is made.
        # It's conceivable that the same column could appear in both tables
        # being compared, but with different case.
        # Though it *may* lead to inconsistencies in these rare cases, this
        # just assumes that there are no duplicated column names in either
        # table, and that the column names can be treated case-insensitively.
        for col in self.common_columns:
            name_lower = col.name.lower()
            if name_lower in ignore_fields:
                continue

            cola = colsa[name_lower]
            colb = colsb[name_lower]

            for attr, _ in _COL_ATTRS:
                vala = getattr(cola, attr, None)
                valb = getattr(colb, attr, None)
                if diff_values(vala, valb):
                    self.diff_column_attributes.append(
                        ((col.name.upper(), attr), (vala, valb))
                    )

            arra = self.a[col.name]
            arrb = self.b[col.name]

            if np.issubdtype(arra.dtype, np.floating) and np.issubdtype(
                arrb.dtype, np.floating
            ):
                diffs = where_not_allclose(arra, arrb, rtol=self.rtol, atol=self.atol)
            elif "P" in col.format or "Q" in col.format:
                diffs = (
                    [
                        idx
                        for idx in range(len(arra))
                        if not np.allclose(
                            arra[idx], arrb[idx], rtol=self.rtol, atol=self.atol
                        )
                    ],
                )
            else:
                diffs = np.where(arra != arrb)

            self.diff_total += len(set(diffs[0]))

            if self.numdiffs >= 0:
                if len(self.diff_values) >= self.numdiffs:
                    # Don't save any more diff values
                    continue

                # Add no more diff'd values than this
                max_diffs = self.numdiffs - len(self.diff_values)
            else:
                max_diffs = len(diffs[0])

            last_seen_idx = None
            for idx in islice(diffs[0], 0, max_diffs):
                if idx == last_seen_idx:
                    # Skip duplicate indices, which my occur when the column
                    # data contains multi-dimensional values; we're only
                    # interested in storing row-by-row differences
                    continue
                last_seen_idx = idx
                self.diff_values.append(((col.name, idx), (arra[idx], arrb[idx])))

        total_values = len(self.a) * len(self.a.dtype.fields)
        self.diff_ratio = float(self.diff_total) / float(total_values)

    def _report(self):
        if self.diff_column_count:
            self._writeln(" Tables have different number of columns:")
            self._writeln(f"  a: {self.diff_column_count[0]}")
            self._writeln(f"  b: {self.diff_column_count[1]}")

        if self.diff_column_names:
            # Show columns with names unique to either table
            for name in self.diff_column_names[0]:
                format = self.diff_columns[0][name.lower()].format
                self._writeln(f" Extra column {name} of format {format} in a")
            for name in self.diff_column_names[1]:
                format = self.diff_columns[1][name.lower()].format
                self._writeln(f" Extra column {name} of format {format} in b")

        col_attrs = dict(_COL_ATTRS)
        # Now go through each table again and show columns with common
        # names but other property differences...
        for col_attr, vals in self.diff_column_attributes:
            name, attr = col_attr
            self._writeln(f" Column {name} has different {col_attrs[attr]}:")
            report_diff_values(
                vals[0],
                vals[1],
                fileobj=self._fileobj,
                indent_width=self._indent + 1,
                rtol=self.rtol,
                atol=self.atol,
            )

        if self.diff_rows:
            self._writeln(" Table rows differ:")
            self._writeln(f"  a: {self.diff_rows[0]}")
            self._writeln(f"  b: {self.diff_rows[1]}")
            self._writeln(" No further data comparison performed.")
            return

        if not self.diff_values:
            return

        # Finally, let's go through and report column data differences:
        for (col, row), values in self.diff_values:
            self._writeln(f" Column {col} data differs in row {row}:")
            report_diff_values(
                values[0],
                values[1],
                fileobj=self._fileobj,
                indent_width=self._indent + 1,
                rtol=self.rtol,
                atol=self.atol,
            )

        if self.diff_values and self.numdiffs < self.diff_total:
            self._writeln(
                f" ...{self.diff_total - self.numdiffs} additional difference(s) found."
            )

        if self.diff_total > self.numdiffs:
            self._writeln(" ...")

        self._writeln(
            f" {self.diff_total} different table data element(s) found "
            f"({self.diff_ratio:.2%} different)."
        )


def report_diff_keyword_attr(fileobj, attr, diffs, keyword, ind=0):
    """
    Write a diff between two header keyword values or comments to the specified
    file-like object.
    """
    if keyword in diffs:
        vals = diffs[keyword]
        for idx, val in enumerate(vals):
            if val is None:
                continue
            if idx == 0:
                dup = ""
            else:
                dup = f"[{idx + 1}]"
            fileobj.write(
                textwrap.indent(
                    f" Keyword {keyword:8}{dup} has different {attr}:\n", ind * "  "
                )
            )
            report_diff_values(val[0], val[1], fileobj=fileobj, indent_width=ind + 1)
