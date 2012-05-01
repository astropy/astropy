Logging for developers
======================

The logging system using the built-in Python `logging
<http://docs.python.org/library/logging.html>`_ module. The logger can be
imported using::

    from astropy.config import log

Messages are then printed using::

    log.debug("This is a debug message")
    log.info("This is an information message")
    log.warning("This is a warning message")

While it is possible to use ``log.error`` and ``log.critical``, such messages
should use Exceptions instead. Note that ``print`` statements should only ever
be used for routines where the user specifically requests printed output, e.g.
``some_object.print_values()``.
