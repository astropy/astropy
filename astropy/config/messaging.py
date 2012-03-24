# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module defines a messaging class similar but more advanced than the
built-in logging class."""

from __future__ import print_function

import sys

from ..utils.misc import find_current_module

# Set fixed levels
DEBUG = 10
INFO = 20
WARN = 30
ERROR = 40

# Set default prefixes
DEBUG_PREFIX = 'DEBUG: '
INFO_PREFIX = 'INFO: '
WARN_PREFIX = 'WARN: '
ERROR_PREFIX = 'ERROR: '
DEBUG_PREFIX_COLOR = '\033[95mDEBUG: \x1b[0m'
INFO_PREFIX_COLOR = '\033[92mINFO: \x1b[0m'
WARN_PREFIX_COLOR = '\033[93mWARN: \x1b[0m'
ERROR_PREFIX_COLOR = '\033[91mERROR: \x1b[0m'

# Set default level
_threshold_level = INFO

# Set whether to use color
_use_color = True


class Message(object):

    def __init__(self, content, level, origin, color=True):

        self.content = content
        self.level = level
        self.origin = origin
        self.color = color

    def __str__(self):

        if self.color:
            if self.level >= ERROR:
                prefix = ERROR_PREFIX_COLOR
            elif self.level >= WARN:
                prefix = WARN_PREFIX_COLOR
            elif self.level >= INFO:
                prefix = INFO_PREFIX_COLOR
            elif self.level >= DEBUG:
                prefix = DEBUG_PREFIX_COLOR
            else:
                prefix = ""
        else:
            if self.level >= ERROR:
                prefix = ERROR_PREFIX
            elif self.level >= WARN:
                prefix = WARN_PREFIX
            elif self.level >= INFO:
                prefix = INFO_PREFIX
            elif self.level >= DEBUG:
                prefix = DEBUG_PREFIX
            else:
                prefix = ""

        return "{:s} {:s}".format(prefix, self.content)


def _show_message(message, level, origin):
    if level >= _threshold_level:
        m = Message(message, level, origin, color=_use_color)
        print(m)


def debug(message):
    return _show_message(message, DEBUG, find_current_module(2))


def info(message):
    return _show_message(message, INFO, find_current_module(2))


def warn(message):
    return _show_message(message, WARN, find_current_module(2))


def error(message):
    return _show_message(message, ERROR, find_current_module(2))


def set_level(level):
    global _threshold_level
    _threshold_level = level


def set_color(use_color):
    global _use_color
    _use_color = use_color


class catch_messages(object):

    def __init__(self, record=True, filter_level=None, filter_origin=None):
        self.record = record
        self.filter_level = filter_level
        self.filter_origin = filter_origin
        self.current_module = sys.modules[__name__]

    def __enter__(self):
        if self.record:
            self._show_message = self.current_module._show_message
            log = []
            def _log_message(message, level, origin):
                global _threshold_level
                if (self.filter_level is None and level >= _threshold_level) or self.filter_level == level:
                    if level >= _threshold_level:
                        if self.filter_origin is None or self.filter_origin == origin:
                            m = Message(message, level, origin, color=_use_color)
                            log.append(m)
                else:
                    if level >= _threshold_level:
                        m = Message(message, level, origin, color=_use_color)
                        print(m)

            self.current_module._show_message = _log_message
            return log
        else:
            return None

    def __exit__(self, *exc_info):
        self.current_module._show_message = self._show_message
