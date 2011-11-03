Introduction to ``astropy.io.vo.ssa``
=====================================

.. note::

   The network functionality in the `~astropy.io.vo.conesearch`,
   `~astropy.io.vo.ssa` and `~astropy.io.vo.image` modules is
   experimental and subject to change.

`astropy.io.vo.ssa` is a Python package to make performing Simple
Spectral Access searches on popular catalogs simple.

This package uses the `VO Simple Spectral Access Protocol
<http://www.ivoa.net/Documents/REC/DAL/SSA-20080201.html>`_.

Simple Usage
------------

Some examples of simple usage follow::

  from vo import ssa
  tab = ssa.query_data(pos=10.2, size=65.4, band='J')

will return a `astropy.io.vo.tree.VOTableFile` object containing a
list of spectral data matching the query.

This list can then be further filtered, and a single row of it passed
to `~astropy.io.vo.ssa.get_data` to get the data itself::

  ssa.get_data(tab.array[0])

.. _range-list-format:

Range-list format
-----------------

Some parameters must be in range-list format as defined in Section
8.7.2 of `Simple Spectral Access Protocol
<http://www.ivoa.net/Documents/REC/DAL/SSA-20080201.html>`_.

These arguments may be a string as passed verbatim to the service
expecting a range-list, or a sequence.  If a sequence, each item must
be either:

    - a numeric value

    - a named value, such as, for example, 'J' for named spectrum
      (if the *numeric* kwarg is False)

    - a 2-tuple indicating a range

    - the last item my be a string indicating the frame of reference

.. TODO: Write a better description of the range-list format
