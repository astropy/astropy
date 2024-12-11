.. _astropy-table-and-dataframes:

Astropy Table and DataFrames
============================

`Pandas <https://pandas.pydata.org/>`_ is a popular data manipulation library for Python
that provides a `~pandas.DataFrame` object which is similar to `astropy.table`. A common
question is why Astropy does not use `~pandas.DataFrame` as the base table object. The
answer stems from a number of domain-specific requirements related to astronomical data
and analysis.

Units and Quantities
--------------------

Astronomy is a physical science, and the data often have units associated with
them. The `astropy.table` package natively supports |Quantity| columns, which are a
powerful way to attach units to array data and perform unit-aware operations. In
addition, the base `~astropy.table.Column` class holds a ``unit`` attribute as
metadata to allow tracking of the units of the data for applications not using
|Quantity|.

*Pandas does not provide support for units.*

Multi-dimensional and Structured Columns
----------------------------------------

Astronomers deal with images, spectra, and other multi-dimensional data that are
commonly stored in a table. An example is a source catalog with an image thumbnail and a
spectrum for each source. Structured columns are less common, but are useful for storing
vectorized data like an `~astropy.coordinates.EarthLocation` in a table.

*Pandas is not able to natively store multi-dimensional or structured columns.*

Lossless representation of FITS and VOTable data via metadata
-------------------------------------------------------------

The `astropy.table` package strives to provide lossless representation of FITS and
VOTable data. This means that when you read a FITS or VOTable file into a table and then
write it back out, the data will be effectively identical. This is made possible by
robust support for table and column metadata which allows storing and propagating common
column information such as the unit, description, and format. For VOTable data, more
information like the UCD is maintained.

*Pandas provides limited support for metadata, but as of late-2024 it is highlighted as
"experimental" in the documentation.*

Time and Coordinates
--------------------

Time and coordinates are fundamental to astronomy, and astropy provides robust support
for them with the `~astropy.time.Time` and `~astropy.coordinates.SkyCoord` classes.
Arrays of times and coordinates can be natively stored in `astropy.table`, meaning that
the full power of these objects is available when working with them as columns within a
table.

*Pandas supports* `timeseries
<https://pandas.pydata.org/docs/user_guide/timeseries.html>`_ *data, but with key
limitations*:

- Leap seconds are not supported. In many circumstances (for instance planning an
  observation) this limitation is not acceptable.
- Pandas times are stored with 64-bit precision, which is not sufficient for some
  astronomical applications. Astropy uses 128-bit precision for time to allow
  sub-nanosecond precision over the age of the universe.
- Different :ref:`time scales <time-scale>` common in astronomy (e.g., TAI, UT1) are
  not supported.
- :ref:`Time formats <time-format>` used in astronomy such as the FITS time format are
  not supported.

*Pandas does not support sky coordinate columns.*

Responsiveness to Community Needs
---------------------------------

The `astropy.table` package is developed by the Astropy community, which is focused on
the needs of astronomers and astrophysicists. This means that the development of the
package can be responsive to the needs of this community and we can develop features
without being constrained by the potential impact to the far broader user base of
Pandas.

Interoperability
----------------

We recognize that Pandas is a popular library and that there are many users who are
familiar with it. For this reason, we have made it easy to convert between
`astropy.table` and `~pandas.DataFrame`, as documented in :ref:`pandas`. This allows
users to take advantage of the features of both packages as needed,
within the limitations stated above.

We are also committed to supporting interoperability with a more generalized concept of
the DataFrame, with packages like `polars <https://pola.rs/>`_ gaining popularity.
