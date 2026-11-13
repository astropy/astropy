**********************************************************
ECSV, HDF5, Parquet, PyArrow CSV, YAML (`astropy.io.misc`)
**********************************************************

The `astropy.io.misc` module contains miscellaneous input/output routines that
do not fit elsewhere, and are often used by other ``astropy`` sub-packages. For
example, `astropy.io.misc.hdf5` contains functions to read/write
:class:`~astropy.table.Table` objects from/to HDF5 files, but these
should not be imported directly by users. Instead, users can access this
functionality via the :class:`~astropy.table.Table` class itself (see
:ref:`table_io`). Routines that are intended to be used directly by users are
listed in the `astropy.io.misc` section.

YAML Serialization
==================

The `astropy.io.misc.yaml` module provides :func:`~astropy.io.misc.yaml.dump`
and :func:`~astropy.io.misc.yaml.load` functions for serializing Astropy and
NumPy objects. NumPy masked arrays preserve their data (including values under
the mask), mask, fill value, and hard-mask setting. The ``numpy.ma.masked``
singleton is also supported::

    >>> import numpy as np
    >>> from astropy.io.misc import yaml
    >>> values = np.ma.array([1, 2], mask=[False, True], fill_value=-99)
    >>> restored = yaml.load(yaml.dump(values))
    >>> restored.data.tolist()
    [1, 2]
    >>> restored.mask.tolist()
    [False, True]
    >>> restored.filled().tolist()
    [1, -99]
    >>> yaml.load(yaml.dump(np.ma.masked)) is np.ma.masked
    True

These objects can also be stored in table or column metadata when writing
ECSV files. Reading their new YAML tags requires a version of Astropy that
supports masked-array serialization.

Structured NumPy arrays and scalars retain their field offsets, titles, item
size, and alignment. Previously written structured-array YAML remains
readable. However, older readers cannot reconstruct the new representation
of packed records nested inside an aligned dtype; these layouts require the
updated reader.

Reference/API
=============

.. toctree::
   :maxdepth: 2

   misc_ref_api
