# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
from collections import OrderedDict

import numpy as np

from .base import IORegistryError, UnifiedIORegistryBase

__all__ = ['UnifiedIORegistry', 'UnifiedInputRegistry', 'UnifiedOutputRegistry']


PATH_TYPES = (str, os.PathLike)  # TODO! include bytes


# -----------------------------------------------------------------------------

class UnifiedInputRegistry(UnifiedIORegistryBase):
    """Read-only Registry."""

    def __init__(self):
        super().__init__()  # set _identifiers
        self._readers = OrderedDict()
        self._registries["read"] = dict(attr="_readers", column="Read")
        self._registries_order = ("read", "identify")

    def get_formats(self, data_class=None, *args):
        return super().get_formats(data_class, filter_on="Read")

    # =========================================================================
    # Read methods

    def register_reader(self, data_format, data_class, function, force=False,
                        priority=0):
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
        if not (data_format, data_class) in self._readers or force:
            self._readers[(data_format, data_class)] = function, priority
        else:
            raise IORegistryError("Reader for format '{}' and class '{}' is "
                              'already defined'
                              ''.format(data_format, data_class.__name__))

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, 'read')

    def unregister_reader(self, data_format, data_class):
        """
        Unregister a reader function

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
            raise IORegistryError("No reader defined for format '{}' and class '{}'"
                                  ''.format(data_format, data_class.__name__))

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, 'read')

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
            format_table_str = self._get_format_table_str(data_class, 'Read')
            raise IORegistryError(
                "No reader defined for format '{}' and class '{}'.\n\nThe "
                "available formats are:\n\n{}".format(
                    data_format, data_class.__name__, format_table_str))

    def read(self, cls, *args, format=None, cache=False, **kwargs):
        """
        Read in data.

        The arguments passed to this method depend on the format.
        """

        ctx = None
        try:
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
                            ctx = get_readable_fileobj(args[0], encoding='binary', cache=cache)
                            fileobj = ctx.__enter__()
                        except OSError:
                            raise
                        except Exception:
                            fileobj = None
                        else:
                            args = [fileobj] + list(args[1:])
                    elif hasattr(args[0], 'read'):
                        path = None
                        fileobj = args[0]

                format = self._get_valid_format(
                    'read', cls, path, fileobj, args, kwargs)

            reader = self.get_reader(format, cls)
            data = reader(*args, **kwargs)

            if not isinstance(data, cls):
                # User has read with a subclass where only the parent class is
                # registered.  This returns the parent class, so try coercing
                # to desired subclass.
                try:
                    data = cls(data)
                except Exception:
                    raise TypeError('could not convert reader output to {} '
                                    'class.'.format(cls.__name__))
        finally:
            if ctx is not None:
                ctx.__exit__(*sys.exc_info())

        return data


# -----------------------------------------------------------------------------

class UnifiedOutputRegistry(UnifiedIORegistryBase):
    """Write-only Registry."""

    def __init__(self):
        super().__init__()
        self._writers = OrderedDict()
        self._registries["write"] = dict(attr="_writers", column="Write")
        self._registries_order = ("write", "identify", )

    def get_formats(self, data_class=None, *args):
        return super().get_formats(data_class, filter_on="Write")

    # =========================================================================
    # Write Methods

    def register_writer(self, data_format, data_class, function, force=False, priority=0):
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
        if not (data_format, data_class) in self._writers or force:
            self._writers[(data_format, data_class)] = function, priority
        else:
            raise IORegistryError("Writer for format '{}' and class '{}' is "
                                  'already defined'
                                  ''.format(data_format, data_class.__name__))

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, 'write')

    def unregister_writer(self, data_format, data_class):
        """
        Unregister a writer function

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
            raise IORegistryError("No writer defined for format '{}' and class '{}'"
                                  ''.format(data_format, data_class.__name__))

        if data_class not in self._delayed_docs_classes:
            self._update__doc__(data_class, 'write')

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
            format_table_str = self._get_format_table_str(data_class, 'Write')
            raise IORegistryError(
                "No writer defined for format '{}' and class '{}'.\n\nThe "
                "available formats are:\n\n{}".format(
                    data_format, data_class.__name__, format_table_str))

    def write(self, data, *args, format=None, **kwargs):
        """
        Write out data.

        The arguments passed to this method depend on the format.
        """

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
                elif hasattr(args[0], 'read'):
                    path = None
                    fileobj = args[0]

            format = self._get_valid_format(
                'write', data.__class__, path, fileobj, args, kwargs)

        writer = self.get_writer(format, data.__class__)
        return writer(data, *args, **kwargs)


# -----------------------------------------------------------------------------

class UnifiedIORegistry(UnifiedInputRegistry, UnifiedOutputRegistry):
    """Unified I/O Registry"""

    def __init__(self):
        super().__init__()
        self._registries_order = ("read", "write", "identify")

    def get_formats(self, data_class=None, readwrite=None):
        """
        Get the list of registered I/O formats as a Table.

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
        return UnifiedIORegistryBase.get_formats(self, data_class, filter_on=readwrite)
