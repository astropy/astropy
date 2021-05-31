# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import warnings

from astropy.io import registry as io_registry
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["CosmologyRead", "CosmologyWrite"]
__doctest_skip__ = ["CosmologyRead", "CosmologyWrite"]


class CosmologyRead(io_registry.UnifiedReadWrite):
    """Read and parse data to a `~astropy.cosmology.Cosmology`.

    This function provides the Cosmology interface to the Astropy unified I/O
    layer. This allows easily reading a file in many supported data formats
    using syntax such as::

        >>> from astropy.cosmology import Cosmology
        >>> cosmo1 = Cosmology.read('cosmo1.ecsv')
        >>> cosmo2 = Cosmology.read('cosmo2.json')

    Get help on the available readers using the ``help()`` method::

      >>> Cosmology.read.help()  # Get help reading and list supported formats
      >>> Cosmology.read.help('json')  # Get detailed help on JSON reader
      >>> Cosmology.read.list_formats()  # Print list of available formats

    See also: https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : (optional)
        Positional arguments passed through to data reader. If supplied the
        first argument is typically the input filename.
        If the first argument is a Mapping all other arguments are ignored.
    format : str (optional, keyword-only)
        File format specifier.
    **kwargs
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `~astropy.cosmology.Cosmology` subclass instance
        `~astropy.cosmology.Cosmology` corresponding to file contents.
        If the table in the file holds multiple rows, the row index must be
        specified.

    Methods
    -------
    from_mapping
        Parse a mapping instance (e.g. dict) into a cosmology instance.
    from_table
        Parse a `~astropy.table.Table` instance into a cosmology instance.

    Notes
    -----
    The cosmology class must be specified by name in the metadata (if a Table)
    or as a field titled ``cosmology``.

    ``read`` may only be called from `~astropy.cosmology.Cosmology`.

    The other parameters are passed to the initialization signature. Parameters
    not in the signature will be moved to the metadata, except if the key is
    already present, in which case a `TypeError` will be raised.

    Metadata can be included on the Table or as a key "meta" for a JSON.

    Raises
    ------
    TypeError
        If read is called not from the Cosmology base class.
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "read")

        from astropy.cosmology.core import Cosmology

        if cls is not Cosmology:
            warnings.warn("``Cosmology.read()`` must be called from the "
                          "``Cosmology`` base class.",
                          category=AstropyUserWarning)

    def __call__(self, *args, index=None, **kwargs):
        from astropy.cosmology.core import Cosmology

        if self._cls is not Cosmology:
            raise TypeError("``Cosmology.read()`` must be called from the "
                            "``Cosmology`` base class.")

        cosmo = io_registry.read(self._cls, *args, index=index, **kwargs)
        return cosmo

    @staticmethod
    def from_mapping(mapping, key=None, *, move_to_meta=False):
        from astropy.cosmology.io import from_mapping

        return from_mapping(mapping, key=key, move_to_meta=move_to_meta)

    @staticmethod
    def from_table(table, index=None, *, move_to_meta=False):
        from astropy.cosmology.io import from_table

        return from_table(table, index=index, move_to_meta=move_to_meta)


class CosmologyWrite(io_registry.UnifiedReadWrite):
    """Write this Cosmology object out in the specified format.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily writing a file in many supported data formats
    using syntax such as::

      >>> from astropy.table import Table
      >>> dat = Table(dict(cosmology="FlatLambdaCDM"], H0=[67], Om0=[0.3]))
      >>> dat.write('table.ecsv', format='ascii.ecsv')

    Get help on the available writers for ``Cosmology`` using the``help()``
    method::

      >>> Table.write.help()  # Get help writing Table and list supported formats
      >>> Table.write.help('JSON')  # Get detailed help on Cosmology JSON writer
      >>> Table.write.list_formats()  # Print list of available formats

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data writer. If supplied the
        first argument is the output filename.
    format : str
        File format specifier.
    **kwargs : dict, optional
        Keyword arguments passed through to data writer.

    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "write")

    def __call__(self, *args, **kwargs):
        io_registry.write(self._instance, *args, **kwargs)

    def to_mapping(self):
        from astropy.cosmology.io import to_mapping

        return to_mapping(self._instance)

    def to_table(self):
        from astropy.cosmology.io import to_table

        return to_table(self._instance)
