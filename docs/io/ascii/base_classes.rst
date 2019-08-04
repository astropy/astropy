.. include:: references.txt

.. _base_class_elements:

Base Class Elements
*******************

The key elements in :mod:`astropy.io.ascii` are:

* :class:`~astropy.io.ascii.Column`: internal storage of column properties and data.
* :class:`Reader <astropy.io.ascii.BaseReader>`: base class to handle reading and writing tables.
* :class:`Inputter <astropy.io.ascii.BaseInputter>`: gets the lines from the table input.
* :class:`Splitter <astropy.io.ascii.BaseSplitter>`: splits the lines into string column values.
* :class:`Header <astropy.io.ascii.BaseHeader>`: initializes output columns based on the table header or user input.
* :class:`Data <astropy.io.ascii.BaseData>`: populates column data from the table.
* :class:`Outputter <astropy.io.ascii.BaseOutputter>`: converts column data to the specified output format (e.g., ``numpy`` structured array).

Each of these elements is an inheritable class with attributes that control the
corresponding functionality. In this way, the large number of tunable
parameters are modularized into manageable groups. In certain places these
attributes are actually functions for handling special cases.
