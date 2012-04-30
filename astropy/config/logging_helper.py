# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module defines a logging class based on the built-in logging module"""

from __future__ import print_function

import os
import sys
import inspect
import logging
import warnings
from contextlib import contextmanager

from . import ConfigurationItem
from ..utils.console import color_print
from ..utils.misc import find_current_module


__all__ = ['log', 'AstropyLogger', 'LoggingError']


class LoggingError(Exception):
    pass

# Read in configuration

LOG_LEVEL = ConfigurationItem('log_level', 'WARN',
                              "Threshold for the logging messages. Logging "
                              "messages that are less severe than this level "
                              "will be ignored. The levels are 'DEBUG', "
                              "'INFO', 'WARNING', 'ERROR'")

USE_COLOR = ConfigurationItem('use_color', True,
                              "Whether to use color for the level names")

LOG_WARNINGS = ConfigurationItem('log_warnings', False,
                                 "Whether to log warnings.warn calls")

LOG_EXCEPTIONS = ConfigurationItem('log_exceptions', False,
                                   "Whether to log exceptions before raising them")

LOG_TO_FILE = ConfigurationItem('log_to_file', True,
                                "Whether to always log messages to a log "
                                "file")

LOG_FILE_PATH = ConfigurationItem('log_file_path', '~/.astropy/astropy.log',
                                  "The file to log messages to")

LOG_FILE_LEVEL = ConfigurationItem('log_file_level', 'WARN',
                                   "Threshold for logging messages to "
                                   "log_file_path")

LOG_FILE_FORMAT = ConfigurationItem('log_file_format', "%(asctime)r, "
                                    "%(origin)r, %(levelname)r, %(message)r",
                                    "Format for log file entries")


class FilterOrigin(object):
    '''A filter for the record origin'''
    def __init__(self, origin):
        self.origin = origin

    def filter(self, record):
        return record.origin.startswith(self.origin)


class ListHandler(logging.Handler):
    '''A handler that can be used to capture the records in a list'''

    def __init__(self, filter_level=None, filter_origin=None):
        logging.Handler.__init__(self)
        self.log_list = []

    def emit(self, record):
        self.log_list.append(record)

Logger = logging.getLoggerClass()


class AstropyLogger(Logger):
    '''
    This class is used to set up the Astropy logging.

    The main functionality added by this class over the built-in
    logging.Logger class is the ability to keep track of the origin of the
    messages, the ability to enable logging of warnings.warn calls and
    exceptions, and the addition of colorized output and context managers to
    easily capture messages to a file or list.
    '''

    def makeRecord(self, name, level, pathname, lineno, msg, args, exc_info, func=None, extra=None):
        if extra is None:
            extra = {}
        if 'origin' not in extra:
            current_module = find_current_module(1, finddiff=[True, 'logging'])
            if current_module is not None:
                extra['origin'] = current_module.__name__
            else:
                extra['origin'] = 'unknown'
        return Logger.makeRecord(self, name, level, pathname, lineno, msg, args, exc_info, func, extra)

    _showwarning_orig = None

    def _showwarning(self, *args, **kwargs):
        self.warn(args[0].message)

    def enable_warnings_logging(self):
        '''
        Enable logging of warnings.warn() calls

        Once called, any subsequent calls to ``warnings.warn()`` are
        redirected to this logger and emitted with level ``WARN``. Note that
        this replaces the output from ``warnings.warn``.

        This can be disabled with ``disable_warnings_logging``.
        '''
        if self._showwarning_orig is not None:
            raise LoggingError("Warnings logging has already been enabled")
        self._showwarning_orig = warnings.showwarning
        warnings.showwarning = self._showwarning

    def disable_warnings_logging(self):
        '''
        Disable logging of warnings.warn() calls

        Once called, any subsequent calls to ``warnings.warn()`` are no longer
        redirected to this logger.

        This can be re-enabled with ``enable_warnings_logging``.
        '''
        if self._showwarning_orig is None:
            raise LoggingError("Warnings logging has not been enabled")
        if warnings.showwarning != self._showwarning:
            raise LoggingError("Cannot disable warnings logging: warnings.showwarning was not set by this logger, or has been overridden")
        warnings.showwarning = self._showwarning_orig
        self._showwarning_orig = None

    _excepthook_orig = None

    def _excepthook(self, type, value, traceback):
        try:
            origin = inspect.getmodule(traceback.tb_next).__name__
        except:
            origin = inspect.getmodule(traceback).__name__
        self.error(value.message, extra={'origin': origin})
        self._excepthook_orig(type, value, traceback)

    def enable_exception_logging(self):
        '''
        Enable logging of exceptions

        Once called, any uncaught exceptions will be emitted with level
        ``ERROR`` by this logger, before being raised.

        This can be disabled with ``disable_exception_logging``.
        '''
        if self._excepthook_orig is not None:
            raise LoggingError("Exception logging has already been enabled")
        self._excepthook_orig = sys.excepthook
        sys.excepthook = self._excepthook

    def disable_exception_logging(self):
        '''
        Disable logging of exceptions

        Once called, any uncaught exceptions will no longer be emitted by this
        logger.

        This can be re-enabled with ``enable_exception_logging``.
        '''
        if self._excepthook_orig is None:
            raise LoggingError("Exception logging has not been enabled")
        if sys.excepthook != self._excepthook:
            raise LoggingError("Cannot disable exception logging: sys.excepthook was not set by this logger, or has been overridden")
        sys.excepthook = self._excepthook_orig
        self._excepthook_orig = None

    def enable_color(self):
        '''
        Enable colorized output
        '''
        self._use_color = True

    def disable_color(self):
        '''
        Disable colorized output
        '''
        self._use_color = False

    def _stream_formatter(self, record):
        '''
        The formatter for standard output
        '''
        if record.levelno < logging.DEBUG or not self._use_color:
            print(record.levelname, end='')
        elif(record.levelno < logging.INFO):
            color_print(record.levelname, 'magenta', end='')
        elif(record.levelno < logging.WARN):
            color_print(record.levelname, 'green', end='')
        elif(record.levelno < logging.ERROR):
            color_print(record.levelname, 'brown', end='')
        else:
            color_print(record.levelname, 'red', end='')
        print(": " + record.msg + " [{:s}]".format(record.origin))

    @contextmanager
    def log_to_file(self, filename, filter_level=None, filter_origin=None):
        '''
        Context manager to temporarily log messages to a file

        Parameters
        ----------
        filename : str
            The file to log messages to.
        filter_level : str
            If set, any log messages less important than ``filter_level`` will
            not be output to the file. Note that this is in addition to the
            top-level filtering for the logger, so if the logger has level
            'INFO', then setting ``filter_level`` to ``INFO`` or ``DEBUG``
            will have no effect, since these messages are already filtered
            out.
        filter_origin : str
            If set, only log messages with an origin starting with
            ``filter_origin`` will be output to the file.

        Notes
        -----

        By default, the logger already outputs log messages to a file set in
        the Astropy configuration file. Using this context manager does not
        stop log messages from being output to that file, nor does it stop log
        messages from being printed to standard output.

        Examples
        --------

        The context manager is used as::

            with logger.log_to_file('myfile.log'):
                # your code here
        '''

        fh = logging.FileHandler(filename)
        if filter_level is not None:
            fh.setLevel(filter_level)
        if filter_origin is not None:
            fh.addFilter(FilterOrigin(filter_origin))
        f = logging.Formatter(LOG_FILE_FORMAT())
        fh.setFormatter(f)
        self.addHandler(fh)
        yield
        self.removeHandler(fh)

    @contextmanager
    def log_to_list(self, filter_level=None, filter_origin=None):
        '''
        Context manager to temporarily log messages to a list

        Parameters
        ----------
        filename : str
            The file to log messages to.
        filter_level : str
            If set, any log messages less important than ``filter_level`` will
            not be output to the file. Note that this is in addition to the
            top-level filtering for the logger, so if the logger has level
            'INFO', then setting ``filter_level`` to ``INFO`` or ``DEBUG``
            will have no effect, since these messages are already filtered
            out.
        filter_origin : str
            If set, only log messages with an origin starting with
            ``filter_origin`` will be output to the file.

        Notes
        -----

        Using this context manager does not stop log messages from being
        output to standard output.

        Examples
        --------

        The context manager is used as::

            with logger.log_to_list() as log_list:
                # your code here
        '''
        lh = ListHandler()
        if filter_level is not None:
            lh.setLevel(filter_level)
        if filter_origin is not None:
            lh.addFilter(FilterOrigin(filter_origin))
        self.addHandler(lh)
        yield lh.log_list
        self.removeHandler(lh)

    def _set_defaults(self):
        '''
        Reset logger to its initial state
        '''

        # Reset any previously installed hooks
        if self._showwarning_orig is not None:
            self.disable_warnings_logging()
        if self._excepthook_orig is not None:
            self.disable_exception_logging()

        # Remove all previous handlers
        for handler in self.handlers[:]:
            self.removeHandler(handler)

        # Set levels
        self.setLevel(LOG_LEVEL())
        if USE_COLOR():
            self.enable_color()
        else:
            self.disable_color()

        # Set up the stdout handler
        sh = logging.StreamHandler()
        sh.emit = self._stream_formatter
        self.addHandler(sh)

        # Set up the main log file handler if requested (but this might fail if
        # configuration directory or log file is not writeable).
        if LOG_TO_FILE():
            try:
                fh = logging.FileHandler(os.path.expanduser(LOG_FILE_PATH()))
            except IOError:
                pass
            else:
                formatter = logging.Formatter(LOG_FILE_FORMAT())
                fh.setFormatter(formatter)
                fh.setLevel(LOG_FILE_LEVEL())
                self.addHandler(fh)

        if LOG_WARNINGS():
            self.enable_warnings_logging()

        if LOG_EXCEPTIONS():
            self.enable_exception_logging()

# Set up the class and initialize logger
_orig_logger_cls = logging.getLoggerClass()
logging.setLoggerClass(AstropyLogger)
try:
    log = logging.getLogger('astropy')
    log._set_defaults()
finally:
    logging.setLoggerClass(_orig_logger_cls)
    del _orig_logger_cls
