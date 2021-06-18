# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import warnings

from astropy.io import registry as io_registry
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["CosmologyRead", "CosmologyWrite"]
__doctest_skip__ = ["CosmologyRead", "CosmologyWrite"]


class CosmologyRead(io_registry.UnifiedReadWrite):
    """Read and parse data to a `~astropy.cosmology.Cosmology`.

    This can *only* be called from `~astropy.cosmology.Cosmology`,
    not any subclasses.

    This function provides the Cosmology interface to the Astropy unified I/O
    layer. This allows easily reading a file in supported data formats using
    syntax such as::

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
    format : str (optional, keyword-only)
        File format specifier.
    **kwargs
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `~astropy.cosmology.Cosmology` subclass instance
        `~astropy.cosmology.Cosmology` corresponding to file contents.

    Methods
    -------
    from_mapping
        Parse a mapping instance (e.g. dict) into a cosmology instance.
    from_table
        Parse a `~astropy.table.Table` instance into a cosmology instance.

    Notes
    -----
    The cosmology class must be specified by name in the metadata (if a |Table|)
    or as a field titled ``cosmology`` if a mapping.
    Metadata can be included on the |Table| or as a key "meta" for a JSON.

    Warns
    -----
    `~astropy.utils.exceptions.AstropyUserWarning`
        If ``read`` is examined not from the Cosmology base class.
    """

    def __new__(cls, instance, cosmo_cls):
        from astropy.cosmology.core import Cosmology

        # warn that ``read`` can only be called from the base class.
        # and return None to prevent usage.
        if cosmo_cls is not Cosmology:
            warnings.warn(("``Cosmology.read()`` must be called from the "
                           "``Cosmology`` base class."),
                           category=AstropyUserWarning)
            return None  # cannot use .read on subclasses, only `Cosmology`.

        return super().__new__(cls)

    def __init__(self, instance, cosmo_cls):
        super().__init__(instance, cosmo_cls, "read")

    def __call__(self, *args, **kwargs):
        cosmo = io_registry.read(self._cls, *args, **kwargs)
        return cosmo

    @staticmethod
    def from_mapping(mapping, *, move_to_meta=False):
        """Load `~astropy.cosmology.Cosmology` from mapping object.

        Parameters
        ----------
        mapping : mapping
            Must have field "cosmology".

        move_to_meta : bool (optional, keyword-only)
            Whether to move arguments not in the initialization signature to the
            metadata. This will only have an effect if there is not variable
            keyword-only argument.

        Returns
        -------
        `~astropy.cosmology.Cosmology` subclass instance

        Raises
        ------
        TypeError
            If not called from the Cosmology base class.

        See Also
        --------
        astropy.cosmology.io.from_mapping
        """
        from astropy.cosmology.io import from_mapping

        return from_mapping(mapping, move_to_meta=move_to_meta)

    @staticmethod
    def from_table(table, index=None, *, move_to_meta=False):
        """Instantiate a `~astropy.cosmology.Cosmology` from a |QTable|.

        Parameters
        ----------
        cosmology : `~astropy.cosmology.Cosmology` class
        table : `~astropy.QTable`
        index : int or None, optional
            The row from table.

        Returns
        -------
        `~astropy.cosmology.Cosmology` subclass instance

        Raises
        ------
        TypeError
            If not called from the Cosmology base class.

        See Also
        --------
        astropy.cosmology.io.from_table
        """
        from astropy.cosmology.io import from_table

        return from_table(table, index=index, move_to_meta=move_to_meta)


class CosmologyWrite(io_registry.UnifiedReadWrite):
    """Write this Cosmology object out in the specified format.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily writing a file in supported data formats
    using syntax such as::

      >>> from astropy.cosmology import Planck18
      >>> Planck18.write('table.ecsv', format='ascii.ecsv')

    Get help on the available writers for ``Cosmology`` using the``help()``
    method::

      >>> Cosmology.write.help()  # Get help writing and list supported formats
      >>> Cosmology.write.help('JSON')  # Get detailed help on Cosmology JSON writer
      >>> Cosmology.write.list_formats()  # Print list of available formats

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
        """Serialize the Cosmology class, inputs, and metadata into a dict.

        Returns
        -------
        dict
            with key-values for the cosmology parameters and also:

            - 'cosmology' : the class
            - 'meta' : the contents of the cosmology's metadata attribute

        See Also
        --------
        astropy.cosmology.io.to_mapping
        """
        from astropy.cosmology.io import to_mapping

        return to_mapping(self._instance)

    def to_table(self):
        """Serialize the Cosmology into a `~astropy.table.QTable`.

        Returns
        -------
        `~astropy.table.QTable`
            with columns for the cosmology parameters, and metadata and
            cosmology class name in the Table's ``meta`` attribute

        See Also
        --------
        astropy.cosmology.io.to_table
        """
        from astropy.cosmology.io import to_table

        return to_table(self._instance)
