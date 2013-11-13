.. _astropy_vo:

*******************************************
Virtual Observatory Access (``astropy.vo``)
*******************************************

.. module:: astropy.vo

Introduction
============

The ``astropy.vo`` subpackage handles simple access for Virtual Observatory
(VO) services.

Current services include:

.. toctree::
   :maxdepth: 1

   conesearch

There are two third-party Python packages related to ``astropy.vo``:

* `PyVO <http://pyvo.readthedocs.org/en/latest/>`_
  provides further functionality to discover
  and query VO services. Its user guide contains a
  `good introduction <https://pyvo.readthedocs.org/en/latest/pyvo/vo.html>`_
  to how the VO works.

* `Astroquery <http://www.astropy.org/astroquery/>`_
  is an Astropy affiliated package that provides simply access to specific astronomical
  web services, many of which do not support the VO protocol.
