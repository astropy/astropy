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

This section demonstrates how to create a custom reader/writer. A reader is
written as a function that can take any arguments except ``format`` (which is
needed when manually specifying the format - see below) and returns an
instance of the :class:`~astropy.table.Table` or
:class:`~astropy.nddata.NDData` classes (or sub-classes). For demonstration,
let's assume that we are trying to write a reader/writer for the
:class:`~astropy.table.Table` class::

    from astropy.table import Table
    def my_table_reader(filename, some_option=1):
        ...  # Read in the table by any means necessary
        return table  # should be an instance of Table

Such a function can then be registered with the I/O registry::

    from astropy.io import registry
    registry.register_reader('my-table-format', Table, my_table_reader)

where the first arguent is the name of the format, the second argument is the
class that the function returns an instance for, and the third argument is the
reader itself.

We can then read in a table with::

    d = Table.read('my_table_file.mtf', format='my-table-format')

In practice, it would be nice to have the ``read`` method automatically
identify that this file is in the ``my-table-format`` format, so we can
construct a function that can recognize these files, which we refer to here as
an *identifier* function.

An identifier function should take a first argument that is a string
which indicates whether the identifier is being called from ``read`` or
``write``, and should then accept arbitrary number of positional and keyword
arguments via ``*args`` and ``**kwargs``, which are the arguments passed to
the ``read`` method.

In the above case, we can write a simplistic function that only looks at
filenames (but in practice, this function could even look at the first few
bytes of the file for example). The only requirement for the identifier
function is that it return a boolean indicating whether the input matches that
expected for the format. In our example, we want to automatically recognize
files with filenames ending in ``.mtf`` as being in the ``my-table-format``
format::

    def identify_mtf(origin, *args, **kwargs):
        return (isinstance(args[0], basestring) and
                args[0].lower().split('.')[-1] in ['mtf'])

.. note:: Identifier functions should be prepared for arbitrary input - in
          particular, the first argument may not be a filename or file
          object, so it should not assume that this is the case.

We then register this identifier function, similarly to the reader function::

    registry.register_identifier('my-table-format', Table, identify_mtf)

Having registered this function, we can then do::

    t = Table.read('catalog.mtf')

If multiple formats match the current input, then an exception is
raised, and similarly if no format matches the current input. In that
case, the format should be explicitly given with the ``format=``
keyword argument.

It is also possible to create custom writers. To go with our custom reader
above, we can write a custom writer::

   def my_table_writer(table, filename, clobber=False):
       ...  # Write the table out to a file

Writer functons should take a dataset object (either an instance of the
:class:`~astropy.table.Table` or :class:`~astropy.nddata.NDData`
classes or sub-classes), and any number of subsequent positional and keyword
arguments - although as for the reader, the ``format`` keyword argument cannot
be used.

We then register the writer::

   registry.register_writer('my-custom-format', Table, my_table_writer)

and we can write the table out to a file::

   t.write('catalog_new.mtf', format='my-table-format')

Since we have already registered the identifier function, we can also simply
do::

   t.write('catalog_new.mtf')

Reference/API
=============

.. automodapi:: astropy.io.registry
