*****************
Tiled Compression
*****************

.. warning::
   This module is in development (so marked as private), anything may change in future releases.
   This documentation is provided to aid in further development of this submodule and related functionality.

This module implements the compression and decompression algorithms, and associated functionality for FITS Tiled Image Compression.
The goal of this submodule is to expose a useful Python API, which different functionality can be built on for reading these files.

The functionality is roughly split up into the following sections:

1. Low level compression and decompression functions implemented in cfitsio (for all algorithms other than the GZIP ones, which use the Python stdlib).
2. The quantize and dequantize functions from cfitsio.
3. A Python C-API module which wraps all the compression and quantize cfitsio functions.
4. `numcodecs <https://numcodecs.readthedocs.io/>`__ style ``Codec`` classes for each compression algorithms.
5. `~astropy.io.fits.hdu.compressed._tiled_compression.compress_image_data` and
   `~astropy.io.fits.hdu.compressed._tiled_compression.decompress_image_data_section` functions which
   are called from `~astropy.io.fits.CompImageHDU`.


.. automodapi:: astropy.io.fits.hdu.compressed._tiled_compression

.. automodapi:: astropy.io.fits.hdu.compressed._codecs

.. automodapi:: astropy.io.fits.hdu.compressed._quantization
