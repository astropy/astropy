**************
Logging system
**************

Overview
========

The Astropy logging system is designed to give users flexibility in deciding
which log messages to show, to capture them, and to send them to a file.

All messages printed by Astropy routines should use the built-in logging
facility (normal ``print()`` calls should only be done by routines that are
explicitly requested to print output). Messages can have one of several
levels:

* DEBUG: Detailed information, typically of interest only when diagnosing
  problems.

* INFO: An message conveying information about the current task, and
  confirming that things are working as expected

* WARNING: An indication that something unexpected happened, and that user
  action may be required.

* ERROR: indicates a more serious issue, including exceptions

By default, only WARNING and ERROR messages are displayed, and are sent to a
log file located at ``~/.astropy/astropy.log`` (if the file is writeable).

Configuring the logging system
==============================

First, import the logger::

    from astropy import log

The threshold level (defined above) for messages can be set with e.g.::

    log.setLevel('INFO')

Color (enabled by default) can be disabled with::

    log.setColor(False)

Warnings from ``warnings.warn`` can be logged with::

    log.enable_warnings_logging()

which can be disabled with::

    log.disable_warnings_logging()

and exceptions can be included in the log with::

    log.set_exception_logging()

which can be disabled with::

    log.disable_exception_logging()

It is also possible to set these settings from the Astropy configuration file,
which also allows an overall log file to be specified. See
`Using the configuration file`_ for more information.

Context managers
================

In some cases, you may want to capture the log messages, for example to check
whether a specific message was output, or to log the messages from a specific
section of code to a file. Both of these are possible using context managers.

To add the log messages to a list, first import the logger if you have not
already done so::

    from astropy import log

then enclose the code in which you want to log the messages to a list in a
``with`` statement::

    with log.log_to_list() as log_list:
        # your code here

In the above example, once the block of code has executed, ``log_list`` will
be a Python list containing all the Astropy logging messages that were raised.
Note that messages continue to be output as normal.

Similarly, you can output the log messages of a specific section of code to a
file using::

    with log.log_to_file('myfile.log'):
        # your code here

which will add all the messages to ``myfile.log`` (this is in addition to the
overall log file mentioned in `Using the configuration file`_).

While these context managers will include all the messages emitted by the
logger (using the global level set by ``log.setLevel``), it is possible to
filter a subset of these using ``filter_level=``, and specifying one of
``'DEBUG'``, ``'INFO'``, ``'WARN'``, ``'ERROR'``. Note that if
``filter_level`` is a lower level than that set via ``setLevel``, only
messages with the level set by ``setLevel`` or higher will be included (i.e.
``filter_level`` is only filtering a subset of the messages normally emitted
by the logger).

Similarly, it is possible to filter a subset of the messages by origin by
specifying ``filter_origin=`` followed by a string. If the origin of a message
starts with that string, the message will be included in the context manager.
For example, ``filter_origin='astropy.wcs'`` will include only messages
emitted in the ``astropy.wcs`` sub-package.

Using the configuration file
============================

Options for the logger can be set in the ``[config.logging_helper]`` section
of the Astropy configuration file::

    [config.logging_helper]

    # Threshold for the logging messages. Logging messages that are less severe
    # than this level will be ignored. The levels are 'DEBUG', 'INFO', 'WARNING',
    # 'ERROR'
    log_level = 'INFO'

    # Whether to use color for the level names
    use_color = True

    # Whether to log warnings.warn calls
    log_warnings = False

    # Whether to log exceptions before raising them
    log_exceptions = False

    # Whether to always log messages to a log file
    log_to_file = True

    # The file to log messages to
    log_file_path = '~/.astropy/astropy.log'

    # Threshold for logging messages to log_file_path
    log_file_level = 'INFO'

    # Format for log file entries
    log_file_format = '%(asctime)s, %(origin)s, %(levelname)s, %(message)s'


Reference/API
=============

.. automodapi:: astropy.logger
    :no-inheritance-diagram:

