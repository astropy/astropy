***********************************************************
Miscellaneous: HDF5, YAML, ASDF, pickle (`astropy.io.misc`)
***********************************************************

The `astropy.io.misc` module contains miscellaneous input/output routines that
do not fit elsewhere, and are often used by other ``astropy`` sub-packages. For
example, `astropy.io.misc.hdf5` contains functions to read/write
:class:`~astropy.table.Table` objects from/to HDF5 files, but these
should not be imported directly by users. Instead, users can access this
functionality via the :class:`~astropy.table.Table` class itself (see
:ref:`table_io`). Routines that are intended to be used directly by users are
listed in the `astropy.io.misc` section.

.. automodapi:: astropy.io.misc
   :headings: =-

.. automodapi:: astropy.io.misc.hdf5
   :headings: =-

.. automodapi:: astropy.io.misc.yaml
   :headings: =-

astropy.io.misc.asdf Package
============================

The **asdf** sub-package contains code that is used to serialize ``astropy``
types so that they can be represented and stored using the Advanced Scientific
Data Format (ASDF).

If both **asdf** and **astropy** are installed, no further configuration is
required in order to process ASDF files that contain **astropy** types. The
**asdf** package has been designed to automatically detect the presence of the
tags defined by **astropy**.

For convenience, users can write `~astropy.table.Table` objects to ASDF files
using the :ref:`table_io`. See :ref:`asdf_io` below.

Documentation on the ASDF Standard can be found `here
<https://asdf-standard.readthedocs.io>`__. Documentation on the ASDF Python
module can be found `here <https://asdf.readthedocs.io>`__. Additional details
for Astropy developers can be found in :ref:`asdf_dev`.

.. _asdf_io:

Using ASDF With Table I/O
-------------------------

ASDF provides readers and writers for `~astropy.table.Table` using the
:ref:`table_io`. This makes it convenient to read and write ASDF files with
`~astropy.table.Table` data.

Basic Usage
^^^^^^^^^^^

Given a table, it is possible to write it out to an ASDF file::

    from astropy.table import Table

    # Create a simple table
    t = Table(dtype=[('a', 'f4'), ('b', 'i4'), ('c', 'S2')])
    # Write the table to an ASDF file
    t.write('table.asdf')

The I/O registry automatically selects the appropriate writer function to use
based on the ``.asdf`` extension of the output file.

Reading a file generated in this way is also possible using
`~astropy.table.Table.read`::

    t2 = Table.read('table.asdf')

The I/O registry automatically selects the appropriate reader function based on
the extension of the input file.

In the case of both reading and writing, if the file extension is not ``.asdf``
it is possible to explicitly specify the reader/writer function to be used::

    t3 = Table.read('table.zxcv', format='asdf')

Advanced Usage
^^^^^^^^^^^^^^

The fundamental ASDF data structure is the tree, which is a nested
combination of basic data structures (see `this
<https://asdf.readthedocs.io/en/latest/asdf/features.html#data-model>`_
for a more detailed description). At the top level, the tree is a `dict`.

The consequence of this is that a `~astropy.table.Table` object (or any object,
for that matter) can be stored at any arbitrary location within an ASDF tree.
The basic writer use case described above stores the given
`~astropy.table.Table` at the top of the tree using a default key. The basic
reader case assumes that a `~astropy.table.Table` is stored in the same place.

However, it may sometimes be useful for users to specify a different top-level
key to be used for storage and retrieval of a `~astropy.table.Table` from an
ASDF file. For this reason, the ASDF I/O interface provides ``data_key`` as an
optional keyword when writing and reading::

    from astropy.table import Table

    t = Table(dtype=[('a', 'f4'), ('b', 'i4'), ('c', 'S2')])
    # Write the table to an ASDF file using a non-default key
    t.write('foo.asdf', data_key='foo')

A `~astropy.table.Table` stored using a custom data key can be retrieved by
passing the same argument to `~astropy.table.Table.read`::

    foo = Table.read('foo.asdf', data_key='foo')

The ``data_key`` option only applies to `~astropy.table.Table` objects that are
stored at the top of the ASDF tree. For full generality, users may pass a
callback when writing or reading ASDF files to define precisely where the
`~astropy.table.Table` object should be placed in the tree. The option for the
write case is ``make_tree``. The function callback should accept exactly one
argument, which is the `~astropy.table.Table` object, and should return a
`dict` representing the tree to be stored::

    def make_custom_tree(table):
        # Return a nested tree where the table is stored at the second level
        return dict(foo=dict(bar=table))

    t = Table(dtype=[('a', 'f4'), ('b', 'i4'), ('c', 'S2')])
    # Write the table to an ASDF file using a non-default key
    t.write('foobar.asdf', make_tree=make_custom_tree)

Similarly, when reading an ASDF file, the user can pass a custom callback to
locate the table within the ASDF tree. The option in this case is
``find_table``. The callback should accept exactly one argument, which is an
`dict` representing the ASDF tree, and it should return a
`~astropy.table.Table` object::

    def find_table(tree):
        # This returns the Table that was stored by the example above
        return tree['foo']['bar']

    foo = Table.read('foobar.asdf', find_table=find_table)

.. _asdf_dev:

Details
-------

The **asdf** sub-package defines classes, referred to as **tags**, that
implement the logic for serialization and deserialization of ``astropy`` types.
Users should never need to refer to tag implementations directly. Their
presence should be entirely transparent when processing ASDF files.

ASDF makes use of abstract data type definitions called **schemas**. The tag
classes provided here are specific implementations of particular schemas. Some
of the tags in ``astropy`` (e.g., those related to transforms) implement schemas
that are defined by the ASDF Standard. In other cases, both the tags and
schemas are defined within ``astropy`` (e.g., those related to many of the
coordinate frames). Documentation of the individual schemas defined by
``astropy`` can be found below in the :ref:`asdf_schemas` section.

Not all ``astropy`` types are currently serializable by ASDF. Attempting to
write unsupported types to an ASDF file will lead to a ``RepresenterError``. In
order to support new types, new tags and schemas must be created. See `Writing
ASDF Extensions <https://asdf.readthedocs.io/en/latest/asdf/extensions.html>`_
for additional details.

Schemas
-------

Documentation for each of the individual ASDF schemas defined by ``astropy``
can be found below.

.. toctree::
   :maxdepth: 2

   asdf-schemas
