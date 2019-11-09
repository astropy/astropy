.. include:: references.txt

.. doctest-skip-all

.. _vo-samp:

*************************************************************
SAMP (Simple Application Messaging Protocol) (`astropy.samp`)
*************************************************************

`astropy.samp` is a Python implementation of the SAMP messaging system.

Simple Application Messaging Protocol (SAMP) is an inter-process communication
system that allows different client programs, usually running on the same
computer, to communicate with each other by exchanging short messages that may
reference external data files. The protocol has been developed within the
International Virtual Observatory Alliance (IVOA) and is understood by many
desktop astronomy tools, including `TOPCAT
<http://www.star.bris.ac.uk/~mbt/topcat/>`_, `SAO DS9 <http://ds9.si.edu/>`_,
and `Aladin <http://aladin.u-strasbg.fr>`_.

So by using the classes in `astropy.samp`, Python code can interact with
other running desktop clients, for instance displaying a named FITS file in DS9,
prompting Aladin to recenter on a given sky position, or receiving a message
identifying the row when a user highlights a plotted point in TOPCAT.

The way the protocol works is that a SAMP "Hub" process must be running on the
local host, and then various client programs can connect to it. Once connected,
these clients can send messages to each other via the hub. The details are
described in the `SAMP standard <http://www.ivoa.net/documents/SAMP/>`_.

`astropy.samp` provides classes both to set up such a hub process, and to
help implement a client that can send and receive messages. It also provides a
stand-alone program ``samp_hub`` which can run a persistent hub in its own
process. Note that setting up the hub from Python is not always necessary, since
various other SAMP-aware applications may start up a hub independently; in most
cases, only one running hub is used during a SAMP session.

The following classes are available in `astropy.samp`:

* |SAMPHubServer|, which is used to instantiate a hub server that clients can
  then connect to.
* |SAMPHubProxy|, which is used to connect to an existing hub (including hubs
  started from other applications such as
  `TOPCAT <http://www.star.bris.ac.uk/~mbt/topcat/>`_).
* |SAMPClient|, which is used to create a SAMP client.
* |SAMPIntegratedClient|, which is the same as |SAMPClient| except that it has
  a self-contained |SAMPHubProxy| to provide a simpler user interface.

`astropy.samp` is a full implementation of `SAMP V1.3
<http://www.ivoa.net/documents/SAMP/20120411/>`_. As well as the Standard
Profile, it supports the Web Profile, which means that it can be used to also
communicate with web SAMP clients; see the `sampjs
<http://astrojs.github.io/sampjs/>`_ library examples for more details.

.. _IVOA Simple Application Messaging Protocol: http://www.ivoa.net/documents/latest/SAMP.html

Using `astropy.samp`
====================

.. toctree::
   :maxdepth: 2

   example_hub
   example_table_image
   example_clients
   advanced_embed_samp_hub

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. automodapi:: astropy.samp

Acknowledgments
===============

This code is adapted from the `SAMPy <https://pypi.org/project/sampy>`__
package written by Luigi Paioro, who has granted the Astropy Project permission
to use the code under a BSD license.
