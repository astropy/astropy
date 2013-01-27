.. _io_registry:

I/O Registry
============

.. note:: The I/O registry is only meant to be used directly by users who want
          to define their own custom readers/writers. Users who want to find
          out more about what formats are supported by
          :class:`~astropy.table.table.Table` by default should see
          :ref:`table_io` (no formats are currently defined for
          :class:`~astropy.nddata.nddata.NDData`, but this will be added in
          future).

The I/O registry is a sub-module used to define the readers/writers available
for the :class:`~astropy.table.table.Table` and
:class:`~astropy.nddata.nddata.NDData` classes.

The following example demonstrates how to create a reader for the
:class:`~astropy.table.table.Table` class. First, we can create a highly
simplistic FITS reader which just reads the data as a structured array::

    from astropy.table import Table

    def fits_reader(filename, hdu=1):
        from astropy.io import fits
        data = fits.open(filename)[hdu].data
        return Table(data)

and then register it with `astropy.table`::

    from astropy.table import io_registry
    io_registry.register_reader('fits', fits_reader)

Reader functions can take any arguments except ``format`` (since this
is reserved for the ``Table.read`` method) and should return a
``Table`` object.

We can then read in a FITS table with::

    t = Table.read('catalog.fits', format='fits')

In practice, it would be nice to have the ``read`` method automatically
identify that this file was a FITS file, so we can construct a function that
can recognize FITS files, which we refer to here as an *identifier*
function. An identifier function should take three arguments: the first
should be a string which indicates whether the identifier is being called
from ``read`` or ``write``, and the second and third are the positional and
keyword arguments passed to ``Table.read`` respectively (and are therefore a
list and a dictionary). We can write a simplistic function that only looks
at filenames (but in practice, this function could even look at the first
few bytes of the file for example). The only requirement is that it return a
boolean indicating whether the input matches that expected for the format::

    def fits_identify(origin, args, kwargs):
        return isinstance(args[0], basestring) and \
               args[0].lower().split('.')[-1] in ['fits', 'fit']

.. note:: Identifier functions should be prepared for arbitrary input - in
          particular, the first argument may not be a filename or file
          object, so it should not assume that this is the case.

We then register this identifier function with ``astropy.table``::

    io_registry.register_identifier('fits', fits_identify)

And we can then do::

    t = Table.read('catalog.fits')

If multiple formats match the current input, then an exception is
raised, and similarly if no format matches the current input. In that
case, the format should be explicitly given with the ``format=``
keyword argument.

Similarly, it is possible to create custom writers. To go with our simplistic FITS reader above, we can write a simplistic FITS writer::

   def fits_writer(table, filename, clobber=False):
       import numpy as np
       from astropy.io import fits
       fits.writeto(filename, np.array(table), clobber=clobber)

We then register the writer::

   io_registry.register_writer('fits', fits_writer)

And we can then write the file out to a FITS file::

   t.write('catalog_new.fits', format='fits')

If we have registered the identifier as above, we can simply do::

   t.write('catalog_new.fits')

Reference/API
=============

.. automodapi:: astropy.io.registry
