.. include:: references.txt

.. doctest-skip-all

.. _vo-samp-example_hub:

*****************************************
Starting and stopping a ``SAMPHubServer``
*****************************************

This module contains classes to create a SAMP Hub and/or interface 
an application to a running SAMP Hub (Standard Profile).

Creating a SAMP hub is very simple. The only thing you have to do is to create
a |SAMPHubServer| instance and start it::

    >>> from astropy.vo import samp
    >>> hub = samp.SAMPHubServer()
    >>> hub.start()

You can also start a hub from the command line using the ``samp_hub`` command
line utility::
 
   $ samp_hub

To get more help on available options for ``samp_hub``::

   $ samp_hub -h

To stop the server, you can simply press control-C.

.. TODO: give info how to stop the server and explain that it is a separate process from Python.