.. include:: references.txt

.. _vo-samp-example_hub:

*****************************************
Starting and stopping a ``SAMPHubServer``
*****************************************

This module contains classes to create a SAMP Hub and/or interface 
an application to a running SAMP Hub (Standard Profile).

SAMPy is very simple to use for the creation of a running SAMP Hub. The only thing
you have to do is to create a |SAMPHubServer| instance and start it:

>>> from astropy.vo import samp
>>> hub = samp.SAMPHubServer()
>>> hub.start()

To start a hub from the command line use the ``sampy`` command line utility::
 
   $ sampy

To get more help on available options for ``sampy``::

   $ sampy -h

TODO: give info how to stop the server and explain that it is a separate process from Python.