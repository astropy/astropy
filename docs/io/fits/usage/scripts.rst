Executable Scripts
------------------

Astropy installs a couple of useful utility programs on your system that are
built with Astropy.

fitsinfo
^^^^^^^^
.. automodule:: astropy.io.fits.scripts.fitsinfo

fitsheader
^^^^^^^^^^
.. automodule:: astropy.io.fits.scripts.fitsheader

fitscheck
^^^^^^^^^
.. automodule:: astropy.io.fits.scripts.fitscheck

With Astropy installed, please run ``fitscheck --help`` to see the full program
usage documentation.

.. _fitsdiff:

fitsdiff
^^^^^^^^

.. currentmodule:: astropy.io.fits

``fitsdiff`` provides a thin command-line wrapper around the :class:`FITSDiff`
interface--it outputs the report from a :class:`FITSDiff` of two FITS files,
and like common diff-like commands returns a 0 status code if no differences
were found, and 1 if differences were found:

With Astropy installed, please run ``fitscheck --help`` to see the full program
usage documentation.
