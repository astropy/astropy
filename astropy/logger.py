# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module defines a logging class based on the built-in logging module"""

from __future__ import print_function

import inspect
import os
import sys
import logging
import warnings
from contextlib import contextmanager

from . import config as _config
from . import conf as _conf
from .extern.six import PY3, text_type
from .utils import find_current_module
from .utils.exceptions import AstropyWarning, AstropyUserWarning

__all__ = ['Conf', 'conf', 'log', 'AstropyLogger', 'LoggingError']

# import the logging levels from logging so that one can do:
# log.setLevel(log.DEBUG), for example
logging_levels = ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL',
                  'FATAL', ]
for level in logging_levels:
    globals()[level] = getattr(logging, level)
__all__ += logging_levels


# Initialize by calling _init_log()
log = None


class LoggingError(Exception):
    """
    This exception is for various errors that occur in the astropy logger,
    typically when activating or deactivating logger-related features.
    """


class _AstLogIPYExc(Exception):
    """
    An exception that is used only as a placeholder to indicate to the
    IPython exception-catching mechanism that the astropy
    exception-capturing is activated. It should not actually be used as
    an exception anywhere.
    """


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.logger`.
    """
    log_level = _config.ConfigItem(
        'INFO',
        "Threshold for the logging messages. Logging "
        "messages that are less severe than this level "
        "will be ignored. The levels are ``'DEBUG'``, "
        "``'INFO'``, ``'WARNING'``, ``'ERROR'``.")
    log_warnings = _config.ConfigItem(
        True,
        "Whether to log `warnings.warn` calls.")
    log_exceptions = _config.ConfigItem(
        False,
        "Whether to log exceptions before raising "
        "them.")
    log_to_file = _config.ConfigItem(
        False,
        "Whether to always log messages to a log "
        "file.")
    log_file_path = _config.ConfigItem(
        '',
        "The file to log messages to. When ``''``, "
        "it defaults to a file ``'astropy.log'`` in "
        "the astropy config directory.")
    log_file_level = _config.ConfigItem(
        'INFO',
        "Threshold for logging messages to "
        "`log_file_path`.")
    log_file_format = _config.ConfigItem(
        "%(asctime)r, "
        "%(origin)r, %(levelname)r, %(message)r",
        "Format for log file entries.")


conf = Conf()


def _init_log():
    """Initializes the Astropy log--in most circumstances this is called
    automatically when importing astropy.
    """

    global log

    orig_logger_cls = logging.getLoggerClass()
    logging.setLoggerClass(AstropyLogger)
    try:
        log = logging.getLogger('astropy')
        log._set_defaults()
    finally:
        logging.setLoggerClass(orig_logger_cls)

    return log


def _teardown_log():
    """Shut down exception and warning logging (if enabled) and clear all
    Astropy loggers from the logging module's cache.

    This involves poking some logging module internals, so much if it is 'at
    your own risk' and is allowed to pass silently if any exceptions occur.
    """

    global log

    if log.exception_logging_enabled():
        log.disable_exception_logging()

    if log.warnings_logging_enabled():
        log.disable_warnings_logging()

    del log

    # Now for the fun stuff...
    try:
        logging._acquireLock()
        try:
            loggerDict = logging.Logger.manager.loggerDict
            for key in loggerDict.keys():
                if key == 'astropy' or key.startswith('astropy.'):
                    del loggerDict[key]
        finally:
            logging._releaseLock()
    except Exception:
        pass


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

    def makeRecord(self, name, level, pathname, lineno, msg, args, exc_info,
                   func=None, extra=None, sinfo=None):
        if extra is None:
            extra = {}
        if 'origin' not in extra:
            current_module = find_current_module(1, finddiff=[True, 'logging'])
            if current_module is not None:
                extra['origin'] = current_module.__name__
            else:
                extra['origin'] = 'unknown'
        if PY3:
            return Logger.makeRecord(self, name, level, pathname, lineno, msg,
                                     args, exc_info, func=func, extra=extra,
                                     sinfo=sinfo)
        else:
            return Logger.makeRecord(self, name, level, pathname, lineno, msg,
                                     args, exc_info, func=func, extra=extra)

    _showwarning_orig = None

    def _showwarning(self, *args, **kwargs):

        # Bail out if we are not catching a warning from Astropy
        if not isinstance(args[0], AstropyWarning):
            return self._showwarning_orig(*args, **kwargs)

        warning = args[0]
        # Deliberately not using isinstance here: We want to display
        # the class name only when it's not the default class,
        # AstropyWarning.  The name of subclasses of AstropyWarning should
        # be displayed.
        if type(warning) not in (AstropyWarning, AstropyUserWarning):
            message = '{0}: {1}'.format(warning.__class__.__name__, args[0])
        else:
            message = str(args[0])

        mod_path = args[2]
        # Now that we have the module's path, we look through
        # sys.modules to find the module object and thus the
        # fully-package-specified module name.  On Python 2, the
        # module.__file__ is the compiled file name, not the .py, so
        # we have to ignore the extension.  On Python 3,
        # module.__file__ is the original source file name, so things
        # are more direct.
        mod_name = None
        if PY3:
            mod_path, ext = os.path.splitext(mod_path)
            for name, mod in list(sys.modules.items()):
                try:
                    # Believe it or not this can fail in some cases:
                    # https://github.com/astropy/astropy/issues/2671
                    path = os.path.splitext(getattr(mod, '__file__', ''))[0]
                except Exception:
                    continue
                if path == mod_path:
                    mod_name = mod.__name__
                    break
        else:  # pragma: py2
            for name, mod in list(sys.modules.items()):
                try:
                    if getattr(mod, '__file__', '') == mod_path:
                        mod_name = mod.__name__
                        break
                except Exception:
                    continue

        if mod_name is not None:
            self.warning(message, extra={'origin': mod_name})
        else:
            self.warning(message)

    def warnings_logging_enabled(self):
        return self._showwarning_orig is not None

    def enable_warnings_logging(self):
        '''
        Enable logging of warnings.warn() calls

        Once called, any subsequent calls to ``warnings.warn()`` are
        redirected to this logger and emitted with level ``WARN``. Note that
        this replaces the output from ``warnings.warn``.

        This can be disabled with ``disable_warnings_logging``.
        '''
        if self.warnings_logging_enabled():
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
        if not self.warnings_logging_enabled():
            raise LoggingError("Warnings logging has not been enabled")
        if warnings.showwarning != self._showwarning:
            raise LoggingError("Cannot disable warnings logging: "
                               "warnings.showwarning was not set by this "
                               "logger, or has been overridden")
        warnings.showwarning = self._showwarning_orig
        self._showwarning_orig = None

    _excepthook_orig = None

    def _excepthook(self, etype, value, traceback):

        if traceback is None:
            mod = None
        else:
            tb = traceback
            while tb.tb_next is not None:
                tb = tb.tb_next
            mod = inspect.getmodule(tb)

        # include the the error type in the message.
        if len(value.args) > 0:
            message = '{0}: {1}'.format(etype.__name__, str(value))
        else:
            message = text_type(etype.__name__)

        if mod is not None:
            self.error(message, extra={'origin': mod.__name__})
        else:
            self.error(message)
        self._excepthook_orig(etype, value, traceback)

    def exception_logging_enabled(self):
        '''
        Determine if the exception-logging mechanism is enabled.

        Returns
        -------
        exclog : bool
            True if exception logging is on, False if not.
        '''
        try:
            ip = get_ipython()
        except NameError:
            ip = None

        if ip is None:
            return self._excepthook_orig is not None
        else:
            return _AstLogIPYExc in ip.custom_exceptions

    def enable_exception_logging(self):
        '''
        Enable logging of exceptions

        Once called, any uncaught exceptions will be emitted with level
        ``ERROR`` by this logger, before being raised.

        This can be disabled with ``disable_exception_logging``.
        '''
        try:
            ip = get_ipython()
        except NameError:
            ip = None

        if self.exception_logging_enabled():
            raise LoggingError("Exception logging has already been enabled")

        if ip is None:
            # standard python interpreter
            self._excepthook_orig = sys.excepthook
            sys.excepthook = self._excepthook
        else:
            # IPython has its own way of dealing with excepthook

            # We need to locally define the function here, because IPython
            # actually makes this a member function of their own class
            def ipy_exc_handler(ipyshell, etype, evalue, tb, tb_offset=None):
                # First use our excepthook
                self._excepthook(etype, evalue, tb)

                # Now also do IPython's traceback
                ipyshell.showtraceback((etype, evalue, tb), tb_offset=tb_offset)

            # now register the function with IPython
            # note that we include _AstLogIPYExc so `disable_exception_logging`
            # knows that it's disabling the right thing
            ip.set_custom_exc((BaseException, _AstLogIPYExc), ipy_exc_handler)

            # and set self._excepthook_orig to a no-op
            self._excepthook_orig = lambda etype, evalue, tb: None

    def disable_exception_logging(self):
        '''
        Disable logging of exceptions

        Once called, any uncaught exceptions will no longer be emitted by this
        logger.

        This can be re-enabled with ``enable_exception_logging``.
        '''
        try:
            ip = get_ipython()
        except NameError:
            ip = None

        if not self.exception_logging_enabled():
            raise LoggingError("Exception logging has not been enabled")

        if ip is None:
            # standard python interpreter
            if sys.excepthook != self._excepthook:
                raise LoggingError("Cannot disable exception logging: "
                                   "sys.excepthook was not set by this logger, "
                                   "or has been overridden")
            sys.excepthook = self._excepthook_orig
            self._excepthook_orig = None
        else:
            # IPython has its own way of dealing with exceptions
            ip.set_custom_exc(tuple(), None)

    def enable_color(self):
        '''
        Enable colorized output
        '''
        _conf.use_color = True

    def disable_color(self):
        '''
        Disable colorized output
        '''
        _conf.use_color = False

    @contextmanager
    def log_to_file(self, filename, filter_level=None, filter_origin=None):
        '''
        Context manager to temporarily log messages to a file.

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
        f = logging.Formatter(conf.log_file_format)
        fh.setFormatter(f)
        self.addHandler(fh)
        yield
        fh.close()
        self.removeHandler(fh)

    @contextmanager
    def log_to_list(self, filter_level=None, filter_origin=None):
        '''
        Context manager to temporarily log messages to a list.

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
        if self.warnings_logging_enabled():
            self.disable_warnings_logging()
        if self.exception_logging_enabled():
            self.disable_exception_logging()

        # Remove all previous handlers
        for handler in self.handlers[:]:
            self.removeHandler(handler)

        # Set levels
        self.setLevel(conf.log_level)

        # Set up the stdout handler
        sh = StreamHandler()
        self.addHandler(sh)

        # Set up the main log file handler if requested (but this might fail if
        # configuration directory or log file is not writeable).
        if conf.log_to_file:
            log_file_path = conf.log_file_path

            # "None" as a string because it comes from config
            try:
                _ASTROPY_TEST_
                testing_mode = True
            except NameError:
                testing_mode = False

            try:
                if log_file_path == '' or testing_mode:
                    log_file_path = os.path.join(
                        _config.get_config_dir(), "astropy.log")
                else:
                    log_file_path = os.path.expanduser(log_file_path)

                fh = logging.FileHandler(log_file_path)
            except (IOError, OSError) as e:
                warnings.warn(
                    'log file {0!r} could not be opened for writing: '
                    '{1}'.format(log_file_path, text_type(e)), RuntimeWarning)
            else:
                formatter = logging.Formatter(conf.log_file_format)
                fh.setFormatter(formatter)
                fh.setLevel(conf.log_file_level)
                self.addHandler(fh)

        if conf.log_warnings:
            self.enable_warnings_logging()

        if conf.log_exceptions:
            self.enable_exception_logging()


class StreamHandler(logging.StreamHandler):
    """
    A specialized StreamHandler that logs INFO and DEBUG messages to
    stdout, and all other messages to stderr.  Also provides coloring
    of the output, if enabled in the parent logger.
    """

    def emit(self, record):
        '''
        The formatter for stderr
        '''
        if record.levelno <= logging.INFO:
            stream = sys.stdout
        else:
            stream = sys.stderr

        if record.levelno < logging.DEBUG or not _conf.use_color:
            print(record.levelname, end='', file=stream)
        else:
            # Import utils.console only if necessary and at the latest because
            # the import takes a significant time [#4649]
            from .utils.console import color_print
            if record.levelno < logging.INFO:
                color_print(record.levelname, 'magenta', end='', file=stream)
            elif record.levelno < logging.WARN:
                color_print(record.levelname, 'green', end='', file=stream)
            elif record.levelno < logging.ERROR:
                color_print(record.levelname, 'brown', end='', file=stream)
            else:
                color_print(record.levelname, 'red', end='', file=stream)
        record.message = "{0} [{1:s}]".format(record.msg, record.origin)
        print(": " + record.message, file=stream)


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
