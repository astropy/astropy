.. include:: references.rst

API documentation
=================

:mod:`astropy.wcs`
------------------

.. automodule:: astropy.wcs.wcs

Classes
-------

.. toctree::
   :maxdepth: 2

   api_wcs.rst
   api_wcsprm.rst
   api_distortion.rst
   api_sip.rst
   api_units.rst
   relax.rst

Testing astropy.wcs
===================

The unit tests are written for use with `nose
<http://code.google.com/p/python-nose/>`.  To run them, install astropy
and then at the commandline::

   nosetests astropy.wcs.tests.test

.. TODO: Move into new testing framework and document accordingly
