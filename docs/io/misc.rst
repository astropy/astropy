**************************************************************
HDF5, Parquet, PyArrow CSV, YAML (`astropy.io.misc`)
**************************************************************

The `astropy.io.misc` module contains miscellaneous input/output routines that
do not fit elsewhere, and are often used by other ``astropy`` sub-packages. For
example, `astropy.io.misc.hdf5` contains functions to read/write
:class:`~astropy.table.Table` objects from/to HDF5 files, but these
should not be imported directly by users. Instead, users can access this
functionality via the :class:`~astropy.table.Table` class itself (see
:ref:`table_io`). Routines that are intended to be used directly by users are
listed in the `astropy.io.misc` section.

Reference/API
=============

.. toctree::
   :maxdepth: 2

   misc_ref_api
