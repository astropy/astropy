# Licensed under a 3-clause BSD style license - see LICENSE.rst

import contextlib
import re
import warnings
from operator import itemgetter

import numpy as np

__all__ = ["IORegistryError"]


class IORegistryError(Exception):
    """Custom error for registry clashes."""


# -----------------------------------------------------------------------------


class _UnifiedIORegistryBase:
    """Base class for registries in Astropy's Unified IO.

    This base class provides identification functions and miscellaneous
    utilities. For an example how to build a registry subclass we suggest
    :class:`~astropy.io.registry.UnifiedInputRegistry`, which enables
    read-only registries. These higher-level subclasses will probably serve
    better as a baseclass, for instance
    :class:`~astropy.io.registry.UnifiedIORegistry` subclasses both
    :class:`~astropy.io.registry.UnifiedInputRegistry` and
    :class:`~astropy.io.registry.UnifiedOutputRegistry` to enable both
    reading from and writing to files.

    .. versionadded:: 5.0

    """

    def __init__(self):
        # registry of identifier functions
        self._identifiers = {}

        # what this class can do: e.g. 'read' &/or 'write'
        self._registries = {}
        self._registries["identify"] = {
            "attr": "_identifiers",
            "column": "Auto-identify",
        }
        self._registries_order = ("identify",)  # match keys in `_registries`

        # If multiple formats are added to one class the update of the docs is quite
        # expensive. Classes for which the doc update is temporarily delayed are added
        # to this set.
        self._delayed_docs_classes = set()

    @property
    def available_registries(self):
        """Available registries.

        Returns
        -------
        ``dict_keys``
        """
        return self._registries.keys()

    def get_formats(self, data_class=None, filter_on=None):
        """
        Get the list of registered formats as a `~astropy.table.Table`.

        Parameters
        ----------
        data_class : class or None, optional
            Filter readers/writer to match data class (default = all classes).
        filter_on : str or None, optional
            Which registry to show. E.g. "identify"
            If None search for both.  Default is None.

        Returns
        -------
        format_table : :class:`~astropy.table.Table`
            Table of available I/O formats.

        Raises
        ------
        ValueError
            If ``filter_on`` is not None nor a registry name.
        """
        from astropy.table import Table

        # set up the column names
        colnames = (
            "Data class",
            "Format",
            *[self._registries[k]["column"] for k in self._registries_order],
            "Deprecated",
        )
        i_dataclass = colnames.index("Data class")
        i_format = colnames.index("Format")
        i_regstart = colnames.index(
            self._registries[self._registries_order[0]]["column"]
        )
        i_deprecated = colnames.index("Deprecated")

        # registries
        regs = set()
        for k in self._registries.keys() - {"identify"}:
            regs |= set(getattr(self, self._registries[k]["attr"]))
        format_classes = sorted(regs, key=itemgetter(0))
        # the format classes from all registries except "identify"

        rows = []
        for fmt, cls in format_classes:
            # see if can skip, else need to document in row
            if data_class is not None and not self._is_best_match(
                data_class, cls, format_classes
            ):
                continue

            # flags for each registry
            has_ = {
                k: "Yes" if (fmt, cls) in getattr(self, v["attr"]) else "No"
                for k, v in self._registries.items()
            }

            # Check if this is a short name (e.g. 'rdb') which is deprecated in
            # favor of the full 'ascii.rdb'.
            ascii_format_class = ("ascii." + fmt, cls)
            # deprecation flag
            deprecated = "Yes" if ascii_format_class in format_classes else ""

            # add to rows
            rows.append(
                (
                    cls.__name__,
                    fmt,
                    *[has_[n] for n in self._registries_order],
                    deprecated,
                )
            )

        # filter_on can be in self_registries_order or None
        if str(filter_on).lower() in self._registries_order:
            index = self._registries_order.index(str(filter_on).lower())
            rows = [row for row in rows if row[i_regstart + index] == "Yes"]
        elif filter_on is not None:
            raise ValueError(
                'unrecognized value for "filter_on": {0}.\n'
                f"Allowed are {self._registries_order} and None."
            )

        # Sorting the list of tuples is much faster than sorting it after the
        # table is created. (#5262)
        if rows:
            # Indices represent "Data Class", "Deprecated" and "Format".
            data = list(
                zip(*sorted(rows, key=itemgetter(i_dataclass, i_deprecated, i_format)))
            )
        else:
            data = None

        # make table
        # need to filter elementwise comparison failure issue
        # https://github.com/numpy/numpy/issues/6784
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=FutureWarning)

            format_table = Table(data, names=colnames)
            if not np.any(format_table["Deprecated"].data == "Yes"):
                format_table.remove_column("Deprecated")

        return format_table

    @contextlib.contextmanager
    def delay_doc_updates(self, cls):
        """Contextmanager to disable documentation updates when registering
        reader and writer. The documentation is only built once when the
        contextmanager exits.

        .. versionadded:: 1.3

        Parameters
        ----------
        cls : class
            Class for which the documentation updates should be delayed.

        Notes
        -----
        Registering multiple readers and writers can cause significant overhead
        because the documentation of the corresponding ``read`` and ``write``
        methods are build every time.

        Examples
        --------
        see for example the source code of ``astropy.table.__init__``.
        """
        self._delayed_docs_classes.add(cls)

        yield

        self._delayed_docs_classes.discard(cls)
        for method in self._registries.keys() - {"identify"}:
            self._update__doc__(cls, method)

    # =========================================================================
    # Identifier methods

    def register_identifier(self, data_format, data_class, identifier, force=False):
        """
        Associate an identifier function with a specific data type.

        Parameters
        ----------
        data_format : str
            The data format identifier. This is the string that is used to
            specify the data type when reading/writing.
        data_class : class
            The class of the object that can be written.
        identifier : function
            A function that checks the argument specified to `read` or `write` to
            determine whether the input can be interpreted as a table of type
            ``data_format``. This function should take the following arguments:

               - ``origin``: A string ``"read"`` or ``"write"`` identifying whether
                 the file is to be opened for reading or writing.
               - ``path``: The path to the file.
               - ``fileobj``: An open file object to read the file's contents, or
                 `None` if the file could not be opened.
               - ``*args``: Positional arguments for the `read` or `write`
                 function.
               - ``**kwargs``: Keyword arguments for the `read` or `write`
                 function.

            One or both of ``path`` or ``fileobj`` may be `None`.  If they are
            both `None`, the identifier will need to work from ``args[0]``.

            The function should return True if the input can be identified
            as being of format ``data_format``, and False otherwise.
        force : bool, optional
            Whether to override any existing function if already present.
            Default is ``False``.

        Examples
        --------
        To set the identifier based on extensions, for formats that take a
        filename as a first argument, you can do for example

        .. code-block:: python

            from astropy.io.registry import register_identifier
            from astropy.table import Table
            def my_identifier(*args, **kwargs):
                return isinstance(args[0], str) and args[0].endswith('.tbl')
            register_identifier('ipac', Table, my_identifier)
            unregister_identifier('ipac', Table)
        """
        if not (data_format, data_class) in self._identifiers or force:  # noqa: E713
            self._identifiers[(data_format, data_class)] = identifier
        else:
            raise IORegistryError(
                f"Identifier for format {data_format!r} and class"
                f" {data_class.__name__!r} is already defined"
            )

    def unregister_identifier(self, data_format, data_class):
        """
        Unregister an identifier function.

        Parameters
        ----------
        data_format : str
            The data format identifier.
        data_class : class
            The class of the object that can be read/written.
        """
        if (data_format, data_class) in self._identifiers:
            self._identifiers.pop((data_format, data_class))
        else:
            raise IORegistryError(
                f"No identifier defined for format {data_format!r} and class"
                f" {data_class.__name__!r}"
            )

    def identify_format(self, origin, data_class_required, path, fileobj, args, kwargs):
        """Loop through identifiers to see which formats match.

        Parameters
        ----------
        origin : str
            A string ``"read`` or ``"write"`` identifying whether the file is to be
            opened for reading or writing.
        data_class_required : object
            The specified class for the result of `read` or the class that is to be
            written.
        path : str or path-like or None
            The path to the file or None.
        fileobj : file-like or None.
            An open file object to read the file's contents, or ``None`` if the
            file could not be opened.
        args : sequence
            Positional arguments for the `read` or `write` function. Note that
            these must be provided as sequence.
        kwargs : dict-like
            Keyword arguments for the `read` or `write` function. Note that this
            parameter must be `dict`-like.

        Returns
        -------
        valid_formats : list
            List of matching formats.
        """
        valid_formats = []
        for data_format, data_class in self._identifiers:
            if self._is_best_match(data_class_required, data_class, self._identifiers):
                if self._identifiers[(data_format, data_class)](
                    origin, path, fileobj, *args, **kwargs
                ):
                    valid_formats.append(data_format)

        return valid_formats

    # =========================================================================
    # Utils

    def _get_format_table_str(self, data_class, filter_on):
        """``get_formats()``, without column "Data class", as a str."""
        format_table = self.get_formats(data_class, filter_on)
        format_table.remove_column("Data class")
        format_table_str = "\n".join(format_table.pformat(max_lines=-1))
        return format_table_str

    def _is_best_match(self, class1, class2, format_classes):
        """Determine if class2 is the "best" match for class1 in the list of classes.

        It is assumed that (class2 in classes) is True.
        class2 is the best match if:

        - ``class1`` is a subclass of ``class2`` AND
        - ``class2`` is the nearest ancestor of ``class1`` that is in classes
          (which includes the case that ``class1 is class2``)
        """
        if issubclass(class1, class2):
            classes = {cls for fmt, cls in format_classes}
            for parent in class1.__mro__:
                if parent is class2:  # class2 is closest registered ancestor
                    return True
                if parent in classes:  # class2 was superseded
                    return False
        return False

    def _get_valid_format(self, mode, cls, path, fileobj, args, kwargs):
        """
        Returns the first valid format that can be used to read/write the data in
        question.  Mode can be either 'read' or 'write'.
        """
        valid_formats = self.identify_format(mode, cls, path, fileobj, args, kwargs)

        if len(valid_formats) == 0:
            format_table_str = self._get_format_table_str(cls, mode.capitalize())
            raise IORegistryError(
                "Format could not be identified based on the"
                " file name or contents, please provide a"
                " 'format' argument.\n"
                f"The available formats are:\n{format_table_str}"
            )
        elif len(valid_formats) > 1:
            return self._get_highest_priority_format(mode, cls, valid_formats)

        return valid_formats[0]

    def _get_highest_priority_format(self, mode, cls, valid_formats):
        """
        Returns the reader or writer with the highest priority. If it is a tie,
        error.
        """
        if mode == "read":
            format_dict = self._readers
            mode_loader = "reader"
        elif mode == "write":
            format_dict = self._writers
            mode_loader = "writer"

        best_formats = []
        current_priority = -np.inf
        for format in valid_formats:
            try:
                _, priority = format_dict[(format, cls)]
            except KeyError:
                # We could throw an exception here, but get_reader/get_writer handle
                # this case better, instead maximally deprioritise the format.
                priority = -np.inf

            if priority == current_priority:
                best_formats.append(format)
            elif priority > current_priority:
                best_formats = [format]
                current_priority = priority

        if len(best_formats) > 1:
            raise IORegistryError(
                "Format is ambiguous - options are:"
                f" {', '.join(sorted(valid_formats, key=itemgetter(0)))}"
            )
        return best_formats[0]

    def _update__doc__(self, data_class, readwrite):
        """
        Update the docstring to include all the available readers / writers for
        the ``data_class.read``/``data_class.write`` functions (respectively).
        Don't update if the data_class does not have the relevant method.
        """
        # abort if method "readwrite" isn't on data_class
        if not hasattr(data_class, readwrite):
            return

        from .interface import UnifiedReadWrite

        FORMATS_TEXT = "The available built-in formats are:"

        # Get the existing read or write method and its docstring
        class_readwrite_func = getattr(data_class, readwrite)

        if not isinstance(class_readwrite_func.__doc__, str):
            # No docstring--could just be test code, or possibly code compiled
            # without docstrings
            return

        lines = class_readwrite_func.__doc__.splitlines()

        # Find the location of the existing formats table if it exists
        sep_indices = [ii for ii, line in enumerate(lines) if FORMATS_TEXT in line]
        if sep_indices:
            # Chop off the existing formats table, including the initial blank line
            chop_index = sep_indices[0]
            lines = lines[:chop_index]

        # Find the minimum indent, skipping the first line because it might be odd
        matches = [re.search(r"(\S)", line) for line in lines[1:]]
        left_indent = " " * min(match.start() for match in matches if match)

        # Get the available unified I/O formats for this class
        # Include only formats that have a reader, and drop the 'Data class' column
        format_table = self.get_formats(data_class, readwrite.capitalize())
        format_table.remove_column("Data class")

        # Get the available formats as a table, then munge the output of pformat()
        # a bit and put it into the docstring.
        new_lines = format_table.pformat(max_lines=-1, max_width=80)
        table_rst_sep = re.sub("-", "=", new_lines[1])
        new_lines[1] = table_rst_sep
        new_lines.insert(0, table_rst_sep)
        new_lines.append(table_rst_sep)

        # Check for deprecated names and include a warning at the end.
        if "Deprecated" in format_table.colnames:
            new_lines.extend(
                [
                    "",
                    "Deprecated format names like ``aastex`` will be "
                    "removed in a future version. Use the full ",
                    "name (e.g. ``ascii.aastex``) instead.",
                ]
            )

        new_lines = [FORMATS_TEXT, ""] + new_lines
        lines.extend([left_indent + line for line in new_lines])

        # Depending on Python version and whether class_readwrite_func is
        # an instancemethod or classmethod, one of the following will work.
        if isinstance(class_readwrite_func, UnifiedReadWrite):
            class_readwrite_func.__class__.__doc__ = "\n".join(lines)
        else:
            try:
                class_readwrite_func.__doc__ = "\n".join(lines)
            except AttributeError:
                class_readwrite_func.__func__.__doc__ = "\n".join(lines)
