# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.io import registry

from .info import serialize_method_as

__all__ = ["TableRead", "TableWrite"]
__doctest_skip__ = ["TableRead", "TableWrite"]


class TableRead(registry.UnifiedReadWrite):
    """Read and parse a data table and return as a Table.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily reading a file in many supported data formats
    using syntax such as::

      >>> from astropy.table import Table
      >>> dat = Table.read('table.dat', format='ascii')
      >>> events = Table.read('events.fits', format='fits')

    Get help on the available readers for ``Table`` using the``help()`` method::

      >>> Table.read.help()  # Get help reading Table and list supported formats
      >>> Table.read.help('fits')  # Get detailed help on Table FITS reader
      >>> Table.read.list_formats()  # Print list of available formats

    See also: https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data reader. If supplied the
        first argument is typically the input filename.
    format : str
        File format specifier.
    units : list, dict, optional
        List or dict of units to apply to columns
    descriptions : list, dict, optional
        List or dict of descriptions to apply to columns
    **kwargs : dict, optional
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `~astropy.table.Table`
        Table corresponding to file contents

    Notes
    -----
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "read", registry=None)
        # uses default global registry

    def __call__(self, *args, **kwargs):
        cls = self._cls
        units = kwargs.pop("units", None)
        descriptions = kwargs.pop("descriptions", None)

        out = self.registry.read(cls, *args, **kwargs)

        # For some readers (e.g., ascii.ecsv), the returned `out` class is not
        # guaranteed to be the same as the desired output `cls`.  If so,
        # try coercing to desired class without copying (io.registry.read
        # would normally do a copy).  The normal case here is swapping
        # Table <=> QTable.
        if cls is not out.__class__:
            try:
                out = cls(out, copy=False)
            except Exception:
                raise TypeError(
                    f"could not convert reader output to {cls.__name__} class."
                )

        out._set_column_attribute("unit", units)
        out._set_column_attribute("description", descriptions)

        return out


class TableWrite(registry.UnifiedReadWrite):
    """
    Write this Table object out in the specified format.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily writing a file in many supported data formats
    using syntax such as::

      >>> from astropy.table import Table
      >>> dat = Table([[1, 2], [3, 4]], names=('a', 'b'))
      >>> dat.write('table.dat', format='ascii')

    Get help on the available writers for ``Table`` using the``help()`` method::

      >>> Table.write.help()  # Get help writing Table and list supported formats
      >>> Table.write.help('fits')  # Get detailed help on Table FITS writer
      >>> Table.write.list_formats()  # Print list of available formats

    The ``serialize_method`` argument is explained in the section on
    `Table serialization methods
    <https://docs.astropy.org/en/latest/io/unified.html#table-serialization-methods>`_.

    See also: https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data writer. If supplied the
        first argument is the output filename.
    format : str
        File format specifier.
    serialize_method : str, dict, optional
        Serialization method specifier for columns.
    **kwargs : dict, optional
        Keyword arguments passed through to data writer.

    Notes
    -----
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "write", registry=None)
        # uses default global registry

    def __call__(self, *args, serialize_method=None, **kwargs):
        instance = self._instance
        with serialize_method_as(instance, serialize_method):
            self.registry.write(instance, *args, **kwargs)
