.. _io_registry:

************************************
I/O Registry (`astropy.io.registry`)
************************************

.. note::

          The I/O registry is only meant to be used directly by users who want
          to define their own custom readers/writers. Users who want to find
          out more about what built-in formats are supported by
          :class:`~astropy.table.Table` by default should see
          :ref:`table_io`. No built-in formats are currently defined for
          :class:`~astropy.nddata.NDData`, but this will be added in
          future).

Introduction
============

The I/O registry is a sub-module used to define the readers/writers available
for the :class:`~astropy.table.Table` and
:class:`~astropy.nddata.NDData` classes.

Using `astropy.io.registry`
===========================

The following example demonstrates how to create a reader for the
:class:`~astropy.table.Table` class. First, we can create a highly
simplistic FITS reader which just reads the data as a structured array::

    from astropy.table import Table

    def fits_table_reader(filename, hdu=1):
        from astropy.io import fits
        data = fits.open(filename)[hdu].data
        return Table(data)

and then register it::

    from astropy.io import registry
    registry.register_reader('fits', Table, fits_table_reader)

Reader functions can take any arguments except ``format`` (since this
is reserved for :func:`~astropy.io.registry.read`) and should return an instance of the class specified as the second argument of ``register_reader`` (:class:`~astropy.table.Table` in the above case.)

We can then read in a FITS table with::

    t = Table.read('catalog.fits', format='fits')

In practice, it would be nice to have the ``read`` method automatically
identify that this file was a FITS file, so we can construct a function that
can recognize FITS files, which we refer to here as an *identifier* function.
An identifier function should take a first argument that should be a string
which indicates whether the identifier is being called from ``read`` or
``write``, and should then accept arbitrary number of positional and keyword
arguments via ``*args`` and ``**kwargs``, which are the arguments passed to
``Table.read``. We can write a simplistic function that only looks at
filenames (but in practice, this function could even look at the first few
bytes of the file for example). The only requirement is that it return a
boolean indicating whether the input matches that expected for the format::

    def fits_identify(origin, *args, **kwargs):
        return isinstance(args[0], basestring) and \
               args[0].lower().split('.')[-1] in ['fits', 'fit']

.. note:: Identifier functions should be prepared for arbitrary input - in
          particular, the first argument may not be a filename or file
          object, so it should not assume that this is the case.

We then register this identifier function::

    registry.register_identifier('fits', Table, fits_identify)

And we can then do::

    t = Table.read('catalog.fits')

If multiple formats match the current input, then an exception is
raised, and similarly if no format matches the current input. In that
case, the format should be explicitly given with the ``format=``
keyword argument.

Similarly, it is possible to create custom writers. To go with our simplistic FITS reader above, we can write a simplistic FITS writer::

   def fits_table_writer(table, filename, clobber=False):
       import numpy as np
       from astropy.io import fits
       fits.writeto(filename, np.array(table), clobber=clobber)

We then register the writer::

   io_registry.register_writer('fits', Table, fits_table_writer)

And we can then write the file out to a FITS file::

   t.write('catalog_new.fits', format='fits')

If we have registered the identifier as above, we can simply do::

   t.write('catalog_new.fits')

Reference/API
=============

.. automodapi:: astropy.io.registry
