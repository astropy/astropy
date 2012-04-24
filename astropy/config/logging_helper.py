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


__all__ = ['log', 'LoggingError']


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

LOG_WARNINGS = ConfigurationItem('catch_warnings', False,
                                 "Whether to log warnings.warn calls")

LOG_EXCEPTIONS = ConfigurationItem('catch_exceptions', False,
                                   "Whether to log exceptions before raising them")

LOG_TO_FILE = ConfigurationItem('log_to_file', True,
                                "Whether to always log messages to a log "
                                "file")

LOG_FILE_PATH = ConfigurationItem('log_file_path', '~/.astropy/astropy.log',
                                  "The file to log messages to")

LOG_FILE_LEVEL = ConfigurationItem('log_file_level', 'WARN',
                                   "Threshold for logging messages to "
                                   "log_file_path")

LOG_FILE_FORMAT = ConfigurationItem('log_file_format', "%(asctime)s, "
                                    "%(origin)s, %(levelname)s, %(message)s",
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

    def makeRecord(self, name, level, pathname, lineno, msg, args, exc_info, func=None, extra=None):
        if extra is None:
            extra = {}
        if 'origin' not in extra:
            extra['origin'] = find_current_module(1, finddiff=[True, 'logging']).__name__
        return Logger.makeRecord(self, name, level, pathname, lineno, msg, args, exc_info, func, extra)

    _showwarning_orig = None

    def _showwarning(self, *args, **kwargs):
        self.warn(args[0].message)

    def enable_warnings_logging(self):
        if self._showwarning_orig is not None:
            raise LoggingError("Warnings logging has already been enabled")
        self._showwarning_orig = warnings.showwarning
        warnings.showwarning = self._showwarning

    def disable_warnings_logging(self):
        if self._showwarning_orig is None:
            raise LoggingError("Warnings logging has not been enabled")
        if warnings.showwarning != self._showwarning:
            raise LoggingError("Cannot disable warnings logging: warnings.showwarning was not set by this logger, or has been overridden")
        warnings.showwarning = self._showwarning_orig
        self._showwarning_orig = None

    _excepthook_orig = None

    def _excepthook(self, type, value, exception):
        try:
            origin = inspect.getmodule(exception.tb_next).__name__
        except:
            origin = inspect.getmodule(exception).__name__
        self.error(value.message, extra={'origin': origin})
        self._excepthook_orig(type, value, exception)

    def enable_exception_logging(self):
        if self._excepthook_orig is not None:
            raise LoggingError("Exception logging has already been enabled")
        self._excepthook_orig = sys.excepthook
        sys.excepthook = self._excepthook

    def disable_exception_logging(self):
        if self._excepthook_orig is None:
            raise LoggingError("Exception logging has not been enabled")
        if sys.excepthook != self._excepthook:
            raise LoggingError("Cannot disable exception logging: sys.excepthook was not set by this logger, or has been overridden")
        sys.excepthook = self._excepthook_orig
        self._excepthook_orig = None

    def setColor(self, use_color):
        self._use_color = use_color

    def stream_formatter(self, record):
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
        lh = ListHandler()
        if filter_level is not None:
            lh.setLevel(filter_level)
        if filter_origin is not None:
            lh.addFilter(FilterOrigin(filter_origin))
        self.addHandler(lh)
        yield lh.log_list
        self.removeHandler(lh)

    def set_defaults(self):

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
        self.setColor(USE_COLOR())

        # Set up the stdout handler
        sh = logging.StreamHandler()
        sh.emit = self.stream_formatter
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

# Set up the class
logging.setLoggerClass(AstropyLogger)

# Initialize logger
log = logging.getLogger('astropy')
log.set_defaults()
