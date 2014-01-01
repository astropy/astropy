.. include:: references.txt

.. doctest-skip-all

.. _vo-samp-example_hub:

*****************************************
Starting and stopping a ``SAMPHubServer``
*****************************************

There are several ways you can start up a SAMP hub:

* By starting up another application that includes a hub, such as
  `TOPCAT <http://www.star.bris.ac.uk/~mbt/topcat/>`_, then connecting
  to it from `astropy.vo.samp`.

* By creating a |SAMPHubServer| instance and starting it, either from the
  interactive Python prompt, or from a Python script::

        >>> from astropy.vo.samp import SAMPHubServer
        >>> hub = SAMPHubServer()
        >>> hub.start()

  You can then stop the hub by calling::

        >>> hub.stop()

* By using the ``samp_hub`` command-line utility, which is included in Astropy::

       $ samp_hub

  To get more help on available options for ``samp_hub``::

       $ samp_hub -h

  To stop the server, you can simply press control-C.
