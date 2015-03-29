I/O mixin
=========

The I/O mixin, `~astropy.nddata.NDIOMixin`, adds ``read`` and ``write``
methods that us the astropy I/O registry.

The mixin itself simply creates the read/write methods; it does not register
any readers or writers with the I/O registry. Subclasses of
`~astropy.nddata.NDDataBase` or `~astropy.nddata.NDData` need to include this
mixin, implement a reader and writer, *and* register it with the I/O
framework. See :ref:`io_registry` for details.
