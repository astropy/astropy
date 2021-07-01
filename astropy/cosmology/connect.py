# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import warnings

from astropy.io import registry as io_registry
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["CosmologyRead", "CosmologyWrite", "CosmologyFrom", "CosmologyTo"]
__doctest_skip__ = __all__


# ==============================================================================
# Read / Write

class CosmologyRead(io_registry.UnifiedReadWrite):
    """Read and parse data to a `~astropy.cosmology.Cosmology`.

    This is (currently) *only* implemented on `~astropy.cosmology.Cosmology`,
    not any subclass.

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
    Metadata can be included on the |Table| or as a key "meta" for a JSON file.

    Warns
    -----
    `~astropy.utils.exceptions.AstropyUserWarning`
        If ``read`` is examined not from the Cosmology base class.
    """

    def __new__(cls, instance, cosmo_cls):
        from astropy.cosmology.core import Cosmology

        # warn that ``read`` is not (yet) implemented for subclasses
        if cosmo_cls is not Cosmology:
            warnings.warn(("``Cosmology.read()`` is not (yet) implemented for "
                           "``Cosmology`` subclasses."),
                           category=AstropyUserWarning)
            return NotImplemented  # TODO! implement for subclasses.

        return super().__new__(cls)

    def __init__(self, instance, cosmo_cls):
        super().__init__(instance, cosmo_cls, "read")

    def __call__(self, *args, **kwargs):
        cosmo = io_registry.read(self._cls, *args, **kwargs)
        return cosmo


class CosmologyWrite(io_registry.UnifiedReadWrite):
    """Write this Cosmology object out in the specified format.

    This function provides the Cosmology interface to the astropy unified I/O
    layer.  This allows easily writing a file in supported data formats
    using syntax such as::

      >>> from astropy.cosmology import Planck18
      >>> Planck18.write('planck18.ecsv')

    Get help on the available writers for ``Cosmology`` using the``help()``
    method::

      >>> Cosmology.write.help()  # Get help writing and list supported formats
      >>> Cosmology.write.help('JSON')  # Get details on Cosmology JSON writer
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


# ==============================================================================
# Convert Between Instances

class CosmologyFrom(io_registry.UnifiedReadWrite):
    """Transform object to a `~astropy.cosmology.Cosmology`."""

    def __new__(cls, instance, cosmo_cls):
        from astropy.cosmology.core import Cosmology

        # warn that ``read`` is not (yet) implemented for subclasses
        if cosmo_cls is not Cosmology:
            warnings.warn(("``Cosmology.read()`` is not (yet) implemented for "
                           "``Cosmology`` subclasses."),
                           category=AstropyUserWarning)
            return NotImplemented  # TODO! implement for subclasses.

        return super().__new__(cls)

    def __init__(self, instance, cosmo_cls):
        super().__init__(instance, cosmo_cls, "read")

    def __call__(self, *args, **kwargs):
        cosmo = io_registry.read(self._cls, *args, **kwargs)
        return cosmo


class CosmologyTo(io_registry.UnifiedReadWrite):
    """Convert this Cosmology object to another thing.
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "write")

    def __call__(self, *args, **kwargs):
        return io_registry.write(self._instance, *args, **kwargs)
