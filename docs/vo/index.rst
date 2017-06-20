.. _astropy_vo:

*******************************************
Virtual Observatory Access (``astropy.vo``)
*******************************************

.. module:: astropy.vo

.. warning::

    The ``astropy.vo`` module has been deprecated and will be removed
    in a future version. The SAMP functionality has been moved
    to ``astropy.samp``, and the conesearch functionality has been
    moved to `Astroquery <http://astroquery.readthedocs.io>`_.

The ``astropy.vo`` subpackage handles simple access for Virtual Observatory
(VO) services.

.. toctree::
   :maxdepth: 1

   conesearch/index

Other third-party Python packages and tools related to ``astropy.vo``:

* `PyVO <http://pyvo.readthedocs.io/en/latest/>`__
  provides further functionality to discover
  and query VO services. Its user guide contains a
  `good introduction <https://pyvo.readthedocs.io/en/latest/pyvo/vo.html>`__
  to how the VO works.

* `Astroquery <http://www.astropy.org/astroquery/>`__
  is an Astropy affiliated package that provides simply access to specific astronomical
  web services, many of which do not support the VO protocol.

* `Simple-Cone-Search-Creator <https://github.com/tboch/Simple-Cone-Search-Creator/>`_
  shows how to ingest a catalog into a cone search service and serve it in VO
  standard format using Python
  (using CSV files and `healpy <https://github.com/healpy/healpy/>`__).
