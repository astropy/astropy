Logging system Documentation
============================

The Astropy logging system is designed to give users the option to filter out
messages depending on their importance.

.. todo:: Move this portion of the documentation to a more appropriate place     
          when a doc reorg happens

For Users
---------

All messages printed by Astropy routines use the built-in logging facility.
Messages can have one of several levels:

* ``DEBUG`` (10): detailed information, typically of interest only when
  diagnosing problems.

* ``INFO`` (20): confirmation that things are working as expected

* ``WARNING`` (30): An indication that something unexpected happened, or
  indicative of some problem in the near future. The program is still working
  as expected.

* ``ERROR`` (40): due to a more serious problem, the program has not been able
  to perform some function (but no exception is being raised).

* ``CRITICAL`` (50): a serious error, indicating that the program itself may
  be unable to continue running (but no exception is being raised).

Note that the ``CRITICAL`` level is unlikely to be used, since critical errors
should raise Exceptions in practice.

It is possible to specify a threshold for logging messages in the astropy
configuration file. Messages which are less severe than this level will be
ignored. The default threshold in Astropy is 20 (``INFO``), indicating that
all the above messages will be shown except ``DEBUG``. Using 40 for example
would cause only ``ERROR`` and ``CRITICAL`` messages to be shown.

.. todo:: Add a link to the page describing the astropy configuration file
          once it exists.

For Developers
--------------

The logging system using the built-in Python `logging
<http://docs.python.org/library/logging.html>`_ module. The logger can be
imported using::

    from astropy.config import logger

Messages are then printed using::

    logger.debug("This is a debug message")
    logger.info("This is an information message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    logger.critical("This is a critical error message")

Note that ``print`` statements should only ever be used for routines where the
user specifically requests printed output, e.g.
``some_object.print_values()``. See `When to use logging
<http://docs.python.org/howto/logging.html#when-to-use-logging>`_ to
understand when to use ``warnings.warn``, the logging facility, and
exceptions.

