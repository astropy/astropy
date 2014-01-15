.. include:: references.txt

.. doctest-skip-all

.. _vo-samp-example_hub:

Starting and stopping a SAMP hub server
---------------------------------------

There are several ways you can start up a SAMP hub:

Using an existing hub
^^^^^^^^^^^^^^^^^^^^^

You can start up another application that includes a hub, such as
`TOPCAT <http://www.star.bris.ac.uk/~mbt/topcat/>`_,
`SAO Ds9 <http://hea-www.harvard.edu/RD/ds9>`_, or
`Aladin Desktop <http://aladin.u-strasbg.fr>`_.

Using the command-line hub utility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can make use of the ``samp_hub`` command-line utility, which is included in
Astropy::

    $ samp_hub

To get more help on available options for ``samp_hub``::

    $ samp_hub -h

To stop the server, you can simply press control-C.

Starting a hub programmatically (advanced)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can start up a hub by creating a |SAMPHubServer| instance and starting it,
either from the interactive Python prompt, or from a Python script::

    >>> from astropy.vo.samp import SAMPHubServer
    >>> hub = SAMPHubServer()
    >>> hub.start()

You can then stop the hub by calling::

    >>> hub.stop()

However, this method is generally not recommended for average users because it
does not work correctly when web SAMP clients try and connect. Instead, this
should be reserved for developers who want to embed a SAMP hub in a GUI for
example. For more information, see :doc:`advanced_embed_samp_hub`.
