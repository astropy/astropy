# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the I/O mixin to the NDData class.


from astropy.io import registry

__all__ = ['NDIOMixin']
__doctest_skip__ = ['NDDataRead', 'NDDataWrite']


class NDDataRead(registry.UnifiedReadWrite):
    """Read and parse gridded N-dimensional data and return as an NDData-derived
    object.

    This function provides the NDDataBase interface to the astropy unified I/O
    layer.  This allows easily reading a file in the supported data formats,
    for example::

      >>> from astropy.nddata import CCDData
      >>> dat = CCDData.read('image.fits')

    Get help on the available readers for ``CCDData`` using the``help()`` method::

      >>> CCDData.read.help()  # Get help reading CCDData and list supported formats
      >>> CCDData.read.help('fits')  # Get detailed help on CCDData FITS reader
      >>> CCDData.read.list_formats()  # Print list of available formats

    See also:

    - https://docs.astropy.org/en/stable/nddata
    - https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data reader. If supplied the
        first argument is the input filename.
    format : str, optional
        File format specifier.
    **kwargs : dict, optional
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `NDData` subclass
        NDData-basd object corresponding to file contents

    Notes
    -----
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'read')

    def __call__(self, *args, **kwargs):
        return registry.read(self._cls, *args, **kwargs)


class NDDataWrite(registry.UnifiedReadWrite):
    """Write this CCDData object out in the specified format.

    This function provides the NDData interface to the astropy unified I/O
    layer.  This allows easily writing a file in many supported data formats
    using syntax such as::

      >>> from astropy.nddata import CCDData
      >>> dat = CCDData(np.zeros((12, 12)), unit='adu')  # 12x12 image of zeros
      >>> dat.write('zeros.fits')

    Get help on the available writers for ``CCDData`` using the``help()`` method::

      >>> CCDData.write.help()  # Get help writing CCDData and list supported formats
      >>> CCDData.write.help('fits')  # Get detailed help on CCDData FITS writer
      >>> CCDData.write.list_formats()  # Print list of available formats

    See also:

    - https://docs.astropy.org/en/stable/nddata
    - https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data writer. If supplied the
        first argument is the output filename.
    format : str, optional
        File format specifier.
    **kwargs : dict, optional
        Keyword arguments passed through to data writer.

    Notes
    -----
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, 'write')

    def __call__(self, *args, **kwargs):
        registry.write(self._instance, *args, **kwargs)


class NDIOMixin:
    """
    Mixin class to connect NDData to the astropy input/output registry.

    This mixin adds two methods to its subclasses, ``read`` and ``write``.
    """
    read = registry.UnifiedReadWriteMethod(NDDataRead)
    write = registry.UnifiedReadWriteMethod(NDDataWrite)
