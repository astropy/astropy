.. include:: references.txt

.. _vo-samp:

***************************************************************
SAMP (Simple Application Messaging Protocol (`astropy.vo.samp`)
***************************************************************

Introduction
============

``astropy.vo.samp`` is an IVOA SAMP (Simple Application Messaging Protocol) messaging system implementation in Python.

Messaging Protocol version 1.3.

It provides classes to easily:

1. instantiate one or multiple Hubs;
2. interface an application or script to a running Hub;
3. create and manage a SAMP client.

``astropy.vo.samp`` provides also a stand-alone program ``sampy`` capable to
instantiate a persistent hub.

.. _IVOA Simple Application Messaging Protocol: http://www.ivoa.net/Documents/latest/SAMP.html


Getting Started
===============

Using `astropy.vo.samp`
=======================

To learn how to use ``astropy.vo.samp`` have a look at these examples: 

.. toctree::
   :maxdepth: 2

   example_hub
   example_clients
   example_table
   example_image

Reference/API
=============

.. automodapi:: astropy.vo.samp
