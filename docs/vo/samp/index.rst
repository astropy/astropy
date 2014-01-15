.. include:: references.txt

.. doctest-skip-all

.. _vo-samp:

***************************************************************
SAMP (Simple Application Messaging Protocol (`astropy.vo.samp`)
***************************************************************

Introduction
============

`astropy.vo.samp` is an IVOA SAMP (Simple Application Messaging Protocol)
messaging system implementation in Python. It provides classes to easily:

1. instantiate one or multiple Hubs;
2. interface an application or script to a running Hub;
3. create and manage a SAMP client.

`astropy.vo.samp` provides also a stand-alone program ``samp_hub`` capable to
instantiate a persistent hub.

SAMP is a protocol that is used by a number of other tools such as
`TOPCAT <http://www.star.bris.ac.uk/~mbt/topcat/>`_,
`SAO Ds9 <http://hea-www.harvard.edu/RD/ds9>`_,
and `Aladin <http://aladin.u-strasbg.fr>`_, which means that it is possible to
send and receive data to and from these tools. The `astropy.vo.samp` package
also supports the 'web profile' for SAMP, which means that it can be used to
communicate with web SAMP clients. See the `sampjs
<http://astrojs.github.io/sampjs/>`_ library examples for more details.

The following classes are available in `astropy.vo.samp`:

* |SAMPHubServer|, which is used to instantiate a hub server that clients can
  then connect to.
* |SAMPHubProxy|, which is used to connect to an existing hub (including hubs
  started from other applications such as
  `TOPCAT <http://www.star.bris.ac.uk/~mbt/topcat/>`_).
* |SAMPClient|, which is used to create a SAMP client
* |SAMPIntegratedClient|, which is the same as |SAMPClient| except that it has
  a self-contained |SAMPHubProxy| to provide a simpler user interface.

.. _IVOA Simple Application Messaging Protocol: http://www.ivoa.net/Documents/latest/SAMP.html

Using `astropy.vo.samp`
=======================

.. toctree::
   :maxdepth: 2

   example_hub
   example_table_image
   example_clients
   advanced_embed_samp_hub

Reference/API
=============

.. automodapi:: astropy.vo.samp

Acknowledgments
===============

This code is adapted from the `SAMPy <https://pypi.python.org/pypi/sampy>`__
package written by Luigi Paioro, who has granted the Astropy project permission
to use the code under a BSD license.
