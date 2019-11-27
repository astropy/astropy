.. include:: references.txt

.. doctest-skip-all

.. _vo-samp-example_hub:

Starting and Stopping a SAMP Hub Server
***************************************

There are several ways you can start up a SAMP hub:

Using an Existing Hub
=====================

You can start up another application that includes a hub, such as
`TOPCAT <http://www.star.bris.ac.uk/~mbt/topcat/>`_,
`SAO Ds9 <http://ds9.si.edu/>`_, or
`Aladin Desktop <http://aladin.u-strasbg.fr>`_.

Using the Command-Line Hub Utility
==================================

You can make use of the ``samp_hub`` command-line utility, which is included in
``astropy``::

    $ samp_hub

To get more help on available options for ``samp_hub``::

    $ samp_hub -h

To stop the server, press control-C.

Starting a Hub Programmatically (Advanced)
==========================================

You can start up a hub by creating a |SAMPHubServer| instance and starting it,
either from the interactive Python prompt, or from a Python script::

    >>> from astropy.samp import SAMPHubServer
    >>> hub = SAMPHubServer()
    >>> hub.start()

You can then stop the hub by calling::

    >>> hub.stop()

However, this method is generally not recommended for average users because it
does not work correctly when web SAMP clients try to connect. Instead, this
should be reserved for developers who want to embed a SAMP hub in a GUI, for
example. For more information, see :doc:`advanced_embed_samp_hub`.
