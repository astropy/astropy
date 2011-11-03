Introduction to ``astropy.io.vo.conesearch``
============================================

.. note::

   The network functionality in the `~astropy.io.vo.conesearch`,
   `~astropy.io.vo.ssa` and `~astropy.io.vo.image` modules is
   experimental and subject to change.

`~astropy.io.vo.conesearch` is a Python package to make performing VO
cone searches on popular catalogs simple.

This package uses the `VO conesearch protocol (v1.03)
<http://www.ivoa.net/Documents/REC/DAL/ConeSearch-20080222.html>`_.

Simple Usage
------------

Some examples of simple usage follow::

  import vo.conesearch as vcone
  tab = vcone.conesearch('USNO-A2', ra=10.2, dec=65.4, sr=0.1)

will return a `astropy.io.vo.tree.VOTableFile` object containing the
contents of the returned list of stars found in the USNO-A2 catalog.
