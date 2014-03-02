.. include:: references.txt

.. _base_class_elements:

Base class elements
----------------------------

The key elements in :mod:`astropy.io.ascii` are:

* :class:`~astropy.io.ascii.Column`: Internal storage of column properties and data ()
* :class:`Reader <astropy.io.ascii.BaseReader>`: Base class to handle reading and writing tables.
* :class:`Inputter <astropy.io.ascii.BaseInputter>`: Get the lines from the table input.
* :class:`Splitter <astropy.io.ascii.BaseSplitter>`: Split the lines into string column values.
* :class:`Header <astropy.io.ascii.BaseHeader>`: Initialize output columns based on the table header or user input.
* :class:`Data <astropy.io.ascii.BaseData>`: Populate column data from the table.
* :class:`Outputter <astropy.io.ascii.BaseOutputter>`: Convert column data to the specified output format, e.g. `numpy` structured array.

Each of these elements is an inheritable class with attributes that control the
corresponding functionality.  In this way the large number of tweakable
parameters is modularized into managable groups.  Where it makes sense these
attributes are actually functions that make it easy to handle special cases.
