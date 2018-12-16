# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.io import registry

from .info import serialize_method_as

__doctest_skip__ = ['TableRead', 'TableWrite']


class TableRead(registry.UnifiedReadWrite):
    """
    Read and parse a data table and return as a Table.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily reading a file in many supported data formats
    using syntax such as::

      >>> from astropy.table import Table
      >>> dat = Table.read('table.dat', format='ascii')
      >>> events = Table.read('events.fits', format='fits')

    See also: http://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    format : str
        File format specifier.
    *args : tuple, optional
        Positional arguments passed through to data reader. If supplied the
        first argument is the input filename.
    **kwargs : dict, optional
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `Table`
        Table corresponding to file contents

    Notes
    -----
    """
    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'read')

    def __call__(self, *args, **kwargs):
        cls = self._cls
        out = registry.read(cls, *args, **kwargs)

        # For some readers (e.g., ascii.ecsv), the returned `out` class is not
        # guaranteed to be the same as the desired output `cls`.  If so,
        # try coercing to desired class without copying (io.registry.read
        # would normally do a copy).  The normal case here is swapping
        # Table <=> QTable.
        if cls is not out.__class__:
            try:
                out = cls(out, copy=False)
            except Exception:
                raise TypeError('could not convert reader output to {0} '
                                'class.'.format(cls.__name__))
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

    See also: http://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    format : str
        File format specifier.
    serialize_method : str, dict, optional
        Serialization method specifier for columns.
    *args : tuple, optional
        Positional arguments passed through to data writer. If supplied the
        first argument is the output filename.
    **kwargs : dict, optional
        Keyword arguments passed through to data writer.

    Notes
    -----
    """
    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'write')

    def __call__(self, *args, **kwargs):
        serialize_method = kwargs.pop('serialize_method', None)
        instance = self._instance
        with serialize_method_as(instance, serialize_method):
            registry.write(instance, *args, **kwargs)
