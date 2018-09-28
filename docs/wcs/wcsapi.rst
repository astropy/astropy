****************************************************
Shared Python interface for World Coordinate Systems
****************************************************

Background
==========

The :class:`~astropy.wcs.WCS` class implements what is considered the
most common 'standard' for representing world coordinate systems in
FITS files, but it cannot represent arbitrarily complex transformations
and there is no agreement on how to use the standard beyond FITS files.
Therefore, other world coordinate system transformation approaches exist,
such as the `gwcs <http://gwcs.readthedocs.io/>`_ package being developed
for the James Webb Space Telescope (which is also applicable to other data).

Since one of the goals of the Astropy Project is to improve interoperability
between packages, we have collaboratively defined a standardized application
programming interface (API) for world coordinate system objects to be used
in Python. This API is described in the Astropy Proposal for Enhancements (APE) 14:
`A shared Python interface for World Coordinate Systems
<https://doi.org/10.5281/zenodo.1188874>`_.

The core astropy package provides base classes that define the low- and high-
level APIs described in APE 14 in the :mod:`astropy.wcs.wcsapi` module, and
these are listed below.

Reference/API
=============

.. automodapi:: astropy.wcs.wcsapi
