Example Usage
=============

Loading WCS information from a FITS file
----------------------------------------

This example loads a FITS file (supplied on the commandline) and uses
the WCS cards in its primary header to transform.

.. literalinclude:: examples/from_file.py
   :language: python

Building a WCS structure programmatically
-----------------------------------------

This example, rather than starting from a FITS header, sets WCS values
programmatically, uses those settings to transform some points, and then
saves those settings to a new FITS header.

.. literalinclude:: examples/programmatic.py
   :language: python

