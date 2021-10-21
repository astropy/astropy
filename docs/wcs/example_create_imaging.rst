.. _example_create_imaging:

First Example
^^^^^^^^^^^^^

This example, rather than starting from a FITS header, sets WCS values
programmatically, uses those settings to transform some points, and then
saves those settings to a new FITS header.

.. literalinclude:: examples/programmatic.py
   :language: python

.. note::
    The members of the WCS object correspond roughly to the key/value
    pairs in the FITS header.  However, they are adjusted and
    normalized in a number of ways that make performing the WCS
    transformation easier.  Therefore, they can not be relied upon to
    get the original values in the header.  To build up a FITS header
    directly and specifically, use `astropy.io.fits.Header` directly.
