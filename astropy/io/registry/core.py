# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
from collections import OrderedDict

from .base import IORegistryError, _UnifiedIORegistryBase

__all__ = ["UnifiedIORegistry", "UnifiedInputRegistry", "UnifiedOutputRegistry"]


PATH_TYPES = (str, os.PathLike)  # TODO! include bytes


def _expand_user_in_args(args):
    # Conservatively attempt to apply `os.path.expanduser` to the first
    # argument, which can be either a path or the contents of a table.
    if len(args) and isinstance(args[0], PATH_TYPES):
        ex_user = os.path.expanduser(args[0])
        if ex_user != args[0] and os.path.exists(os.path.dirname(ex_user)):
            args = (ex_user,) + args[1:]
    return args


# -----------------------------------------------------------------------------


class UnifiedInputRegistry(_UnifiedIORegistryBase):
    """Read-only Unified Registry.

    .. versionadded:: 5.0

    Examples
    --------
    First let's start by creating a read-only registry.

    .. code-block:: python

        >>> from astropy.io.registry import UnifiedInputRegistry
        >>> read_reg = UnifiedInputRegistry()

    There is nothing in this registry. Let's make a reader for the
    :class:`~astropy.table.Table` class::

        from astropy.table import Table

        def my_table_reader(filename, some_option=1):
            # Read in the table by any means necessary
            return table  # should be an instance of Table

    Such a function can then be registered with the I/O registry::

        read_reg.register_reader('my-table-format', Table, my_table_reader)

    Note that we CANNOT then read in a table with::

        d = Table.read('my_table_file.mtf', format='my-table-format')

    Why? because ``Table.read`` uses Astropy's default global registry and this
    is a separate registry.
    Instead we can read by the read method on the registry::

        d = read_reg.read(Table, 'my_table_file.mtf', format='my-table-format')

    """

    def __init__(self):
        super().__init__()  # set _identifiers
        self._readers = OrderedDict()
        self._registries["read"] = dict(attr="_readers", column="Read")
        self._registries_order = ("read", "identify")

    # =========================================================================
    # Read methods

    def register_reader(
        self, data_format, data_class, function, force=False, priority=0
    ):
        """
        Register a reader function.

        Parameters
        ----------
        data_format : str
            The data format identifier. This is the string that will be used to
            specify the data type when reading.
        data_class : class
            The class of the object that the reader produces.
        function : function
            The function to read in a data object.
        force : bool, optional
            Whether to override any existing function if already present.
            Default is ``False``.
        priority : int, optional
            The priority of the reader, used to compare possible formats when
            trying to determine the best reader to use. Higher priorities are
            preferred over lower priorities, with the default priority being 0
            (negative numbers are allowed though).
        """
        if (data_format, data_class) not in self._readers or force:
            self._readers[(data_format, data_class)] = function, priority
        else:
            raise IORegistryError(
                f"Reader for format '{data_format}' and class '{data_class.__name__}'"
                " is already defined"
            )

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, "read")

    def unregister_reader(self, data_format, data_class):
        """
        Unregister a reader function.

        Parameters
        ----------
        data_format : str
            The data format identifier.
        data_class : class
            The class of the object that the reader produces.
        """
        if (data_format, data_class) in self._readers:
            self._readers.pop((data_format, data_class))
        else:
            raise IORegistryError(
                f"No reader defined for format '{data_format}' and class"
                f" '{data_class.__name__}'"
            )

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, "read")

    def get_reader(self, data_format, data_class):
        """Get reader for ``data_format``.

        Parameters
        ----------
        data_format : str
            The data format identifier. This is the string that is used to
            specify the data type when reading/writing.
        data_class : class
            The class of the object that can be written.

        Returns
        -------
        reader : callable
            The registered reader function for this format and class.
        """
        readers = [(fmt, cls) for fmt, cls in self._readers if fmt == data_format]
        for reader_format, reader_class in readers:
            if self._is_best_match(data_class, reader_class, readers):
                return self._readers[(reader_format, reader_class)][0]
        else:
            format_table_str = self._get_format_table_str(data_class, "Read")
            raise IORegistryError(
                f"No reader defined for format '{data_format}' and class"
                f" '{data_class.__name__}'.\n\nThe available formats"
                f" are:\n\n{format_table_str}"
            )

    def read(self, cls, *args, format=None, cache=False, **kwargs):
        """
        Read in data.

        Parameters
        ----------
        cls : class
        *args
            The arguments passed to this method depend on the format.
        format : str or None
        cache : bool
            Whether to cache the results of reading in the data.
        **kwargs
            The arguments passed to this method depend on the format.

        Returns
        -------
        object or None
            The output of the registered reader.
        """
        ctx = None
        try:
            # Expand a tilde-prefixed path if present in args[0]
            args = _expand_user_in_args(args)

            if format is None:
                path = None
                fileobj = None

                if len(args):
                    if isinstance(args[0], PATH_TYPES) and not os.path.isdir(args[0]):
                        from astropy.utils.data import get_readable_fileobj

                        # path might be a os.PathLike object
                        if isinstance(args[0], os.PathLike):
                            args = (os.fspath(args[0]),) + args[1:]
                        path = args[0]
                        try:
                            ctx = get_readable_fileobj(
                                args[0], encoding="binary", cache=cache
                            )
                            fileobj = ctx.__enter__()
                        except OSError:
                            raise
                        except Exception:
                            fileobj = None
                        else:
                            args = [fileobj] + list(args[1:])
                    elif hasattr(args[0], "read"):
                        path = None
                        fileobj = args[0]

                format = self._get_valid_format(
                    "read", cls, path, fileobj, args, kwargs
                )

            reader = self.get_reader(format, cls)
            data = reader(*args, **kwargs)

            if not isinstance(data, cls):
                # User has read with a subclass where only the parent class is
                # registered.  This returns the parent class, so try coercing
                # to desired subclass.
                try:
                    data = cls(data)
                except Exception:
                    raise TypeError(
                        f"could not convert reader output to {cls.__name__} class."
                    )
        finally:
            if ctx is not None:
                ctx.__exit__(*sys.exc_info())

        return data


# -----------------------------------------------------------------------------


class UnifiedOutputRegistry(_UnifiedIORegistryBase):
    """Write-only Registry.

    .. versionadded:: 5.0
    """

    def __init__(self):
        super().__init__()
        self._writers = OrderedDict()
        self._registries["write"] = dict(attr="_writers", column="Write")
        self._registries_order = ("write", "identify")

    # =========================================================================
    # Write Methods

    def register_writer(
        self, data_format, data_class, function, force=False, priority=0
    ):
        """
        Register a table writer function.

        Parameters
        ----------
        data_format : str
            The data format identifier. This is the string that will be used to
            specify the data type when writing.
        data_class : class
            The class of the object that can be written.
        function : function
            The function to write out a data object.
        force : bool, optional
            Whether to override any existing function if already present.
            Default is ``False``.
        priority : int, optional
            The priority of the writer, used to compare possible formats when trying
            to determine the best writer to use. Higher priorities are preferred
            over lower priorities, with the default priority being 0 (negative
            numbers are allowed though).
        """
        if not (data_format, data_class) in self._writers or force:  # noqa: E713
            self._writers[(data_format, data_class)] = function, priority
        else:
            raise IORegistryError(
                f"Writer for format '{data_format}' and class '{data_class.__name__}'"
                " is already defined"
            )

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, "write")

    def unregister_writer(self, data_format, data_class):
        """
        Unregister a writer function.

        Parameters
        ----------
        data_format : str
            The data format identifier.
        data_class : class
            The class of the object that can be written.
        """
        if (data_format, data_class) in self._writers:
            self._writers.pop((data_format, data_class))
        else:
            raise IORegistryError(
                f"No writer defined for format '{data_format}' and class"
                f" '{data_class.__name__}'"
            )

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, "write")

    def get_writer(self, data_format, data_class):
        """Get writer for ``data_format``.

        Parameters
        ----------
        data_format : str
            The data format identifier. This is the string that is used to
            specify the data type when reading/writing.
        data_class : class
            The class of the object that can be written.

        Returns
        -------
        writer : callable
            The registered writer function for this format and class.
        """
        writers = [(fmt, cls) for fmt, cls in self._writers if fmt == data_format]
        for writer_format, writer_class in writers:
            if self._is_best_match(data_class, writer_class, writers):
                return self._writers[(writer_format, writer_class)][0]
        else:
            format_table_str = self._get_format_table_str(data_class, "Write")
            raise IORegistryError(
                f"No writer defined for format '{data_format}' and class"
                f" '{data_class.__name__}'.\n\nThe available formats"
                f" are:\n\n{format_table_str}"
            )

    def write(self, data, *args, format=None, **kwargs):
        """
        Write out data.

        Parameters
        ----------
        data : object
            The data to write.
        *args
            The arguments passed to this method depend on the format.
        format : str or None
        **kwargs
            The arguments passed to this method depend on the format.

        Returns
        -------
        object or None
            The output of the registered writer. Most often `None`.

            .. versionadded:: 4.3
        """
        # Expand a tilde-prefixed path if present in args[0]
        args = _expand_user_in_args(args)

        if format is None:
            path = None
            fileobj = None
            if len(args):
                if isinstance(args[0], PATH_TYPES):
                    # path might be a os.PathLike object
                    if isinstance(args[0], os.PathLike):
                        args = (os.fspath(args[0]),) + args[1:]
                    path = args[0]
                    fileobj = None
                elif hasattr(args[0], "read"):
                    path = None
                    fileobj = args[0]

            format = self._get_valid_format(
                "write", data.__class__, path, fileobj, args, kwargs
            )

        writer = self.get_writer(format, data.__class__)
        return writer(data, *args, **kwargs)


# -----------------------------------------------------------------------------


class UnifiedIORegistry(UnifiedInputRegistry, UnifiedOutputRegistry):
    """Unified I/O Registry.

    .. versionadded:: 5.0
    """

    def __init__(self):
        super().__init__()
        self._registries_order = ("read", "write", "identify")

    def get_formats(self, data_class=None, readwrite=None):
        """
        Get the list of registered I/O formats as a `~astropy.table.Table`.

        Parameters
        ----------
        data_class : class, optional
            Filter readers/writer to match data class (default = all classes).

        readwrite : str or None, optional
            Search only for readers (``"Read"``) or writers (``"Write"``).
            If None search for both.  Default is None.

            .. versionadded:: 1.3

        Returns
        -------
        format_table : :class:`~astropy.table.Table`
            Table of available I/O formats.
        """
        return super().get_formats(data_class, readwrite)
