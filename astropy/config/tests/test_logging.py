import sys
import warnings
from tempfile import NamedTemporaryFile

import pytest

from .. import log, LoggingError


# Save original values of hooks. These are not the system values, but the
# already overwritten values since the logger already gets imported before
# this file gets executed.
_excepthook = sys.__excepthook__
_showwarning = warnings.showwarning


def setup_function(function):

    # Reset hooks to original values
    sys.excepthook = _excepthook
    warnings.showwarning = _showwarning

    # Reset internal original hooks
    log._showwarning_orig = None
    log._excepthook_orig = None

    # Set up the logger
    log._set_defaults()

    # Reset hooks
    if log.warnings_logging_enabled():
        log.disable_warnings_logging()
    if log.exception_logging_enabled():
        log.disable_exception_logging()


def teardown_module(function):

    # Ensure that hooks are restored to original values
    sys.excepthook = _excepthook
    warnings.showwarning = _showwarning


def test_warnings_logging_disable_no_enable():
    with pytest.raises(LoggingError) as e:
        log.disable_warnings_logging()
    assert e.value.args[0] == 'Warnings logging has not been enabled'


def test_warnings_logging_enable_twice():
    log.enable_warnings_logging()
    with pytest.raises(LoggingError) as e:
        log.enable_warnings_logging()
    assert e.value.args[0] == 'Warnings logging has already been enabled'


def test_warnings_logging_overridden():
    log.enable_warnings_logging()
    warnings.showwarning = lambda: None
    with pytest.raises(LoggingError) as e:
        log.disable_warnings_logging()
    assert e.value.args[0] == 'Cannot disable warnings logging: warnings.showwarning was not set by this logger, or has been overridden'


def test_warnings_logging():

    # Without warnings logging
    with warnings.catch_warnings(record=True) as warn_list:
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning")
    assert len(log_list) == 0
    assert len(warn_list) == 1
    assert warn_list[0].message.args[0] == "This is a warning"

    # With warnings logging
    with warnings.catch_warnings(record=True) as warn_list:
        log.enable_warnings_logging()
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning")
        log.disable_warnings_logging()
    assert len(log_list) == 1
    assert len(warn_list) == 0
    assert log_list[0].levelname == 'WARNING'
    assert log_list[0].message == 'This is a warning'
    assert log_list[0].origin == 'astropy.config.tests.test_logging'

    # Without warnings logging
    with warnings.catch_warnings(record=True) as warn_list:
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning")
    assert len(log_list) == 0
    assert len(warn_list) == 1
    assert warn_list[0].message.args[0] == "This is a warning"


def test_exception_logging_disable_no_enable():
    with pytest.raises(LoggingError) as e:
        log.disable_exception_logging()
    assert e.value.args[0] == 'Exception logging has not been enabled'


def test_exception_logging_enable_twice():
    log.enable_exception_logging()
    with pytest.raises(LoggingError) as e:
        log.enable_exception_logging()
    assert e.value.args[0] == 'Exception logging has already been enabled'


def test_exception_logging_overridden():
    log.enable_exception_logging()
    sys.excepthook = lambda: None
    with pytest.raises(LoggingError) as e:
        log.disable_exception_logging()
    assert e.value.args[0] == 'Cannot disable exception logging: sys.excepthook was not set by this logger, or has been overridden'


def test_exception_logging():

    # Without exception logging
    try:
        with log.log_to_list() as log_list:
            raise Exception("This is an Exception")
    except Exception as exc:
        sys.excepthook(*sys.exc_info())
        assert exc.args[0] == "This is an Exception"
    else:
        assert False  # exception should have been raised
    assert len(log_list) == 0

    # With exception logging
    try:
        log.enable_exception_logging()
        with log.log_to_list() as log_list:
            raise Exception("This is an Exception")
    except Exception as exc:
        sys.excepthook(*sys.exc_info())
        assert exc.args[0] == "This is an Exception"
    else:
        assert False  # exception should have been raised
    assert len(log_list) == 1
    assert log_list[0].levelname == 'ERROR'
    assert log_list[0].message == 'This is an Exception'
    assert log_list[0].origin == 'astropy.config.tests.test_logging'

    # Without exception logging
    log.disable_exception_logging()
    try:
        with log.log_to_list() as log_list:
            raise Exception("This is an Exception")
    except Exception as exc:
        sys.excepthook(*sys.exc_info())
        assert exc.args[0] == "This is an Exception"
    else:
        assert False  # exception should have been raised
    assert len(log_list) == 0


@pytest.mark.parametrize(('level'), [None, 'DEBUG', 'INFO', 'WARN', 'ERROR'])
def test_log_to_list(level):

    if level is not None:
        log.setLevel(level)

    with log.log_to_list() as log_list:
        log.error("Error message")
        log.warn("Warning message")
        log.info("Information message")
        log.debug("Debug message")

    # Check list length
    if level == 'DEBUG':
        assert len(log_list) == 4
    elif level is None or level == 'INFO':
        assert len(log_list) == 3
    elif level == 'WARN':
        assert len(log_list) == 2
    elif level == 'ERROR':
        assert len(log_list) == 1

    # Check list content

    assert log_list[0].levelname == 'ERROR'
    assert log_list[0].message == 'Error message'
    assert log_list[0].origin == 'astropy.config.tests.test_logging'

    if len(log_list) >= 2:
        assert log_list[1].levelname == 'WARNING'
        assert log_list[1].message == 'Warning message'
        assert log_list[1].origin == 'astropy.config.tests.test_logging'

    if len(log_list) >= 3:
        assert log_list[2].levelname == 'INFO'
        assert log_list[2].msg == 'Information message'
        assert log_list[2].origin == 'astropy.config.tests.test_logging'

    if len(log_list) >= 4:
        assert log_list[3].levelname == 'DEBUG'
        assert log_list[3].msg == 'Debug message'
        assert log_list[3].origin == 'astropy.config.tests.test_logging'


def test_log_to_list_level():

    with log.log_to_list(filter_level='ERROR') as log_list:
        log.error("Error message")
        log.warn("Warning message")

    assert len(log_list) == 1 and log_list[0].levelname == 'ERROR'


def test_log_to_list_origin1():

    with log.log_to_list(filter_origin='astropy.config.tests') as log_list:
        log.error("Error message")
        log.warn("Warning message")

    assert len(log_list) == 2


def test_log_to_list_origin2():

    with log.log_to_list(filter_origin='astropy.wcs') as log_list:
        log.error("Error message")
        log.warn("Warning message")

    assert len(log_list) == 0


@pytest.mark.parametrize(('level'), [None, 'DEBUG', 'INFO', 'WARN', 'ERROR'])
def test_log_to_file(level):

    log_file = NamedTemporaryFile()

    if level is not None:
        log.setLevel(level)

    with log.log_to_file(log_file.name):
        log.error("Error message")
        log.warn("Warning message")
        log.info("Information message")
        log.debug("Debug message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    # Check list length
    if level == 'DEBUG':
        assert len(log_entries) == 4
    elif level is None or level == 'INFO':
        assert len(log_entries) == 3
    elif level == 'WARN':
        assert len(log_entries) == 2
    elif level == 'ERROR':
        assert len(log_entries) == 1

    # Check list content

    assert log_entries[0].endswith(b"'astropy.config.tests.test_logging', 'ERROR', 'Error message'\n")

    if len(log_entries) >= 2:
        assert log_entries[1].endswith(b"'astropy.config.tests.test_logging', 'WARNING', 'Warning message'\n")

    if len(log_entries) >= 3:
        assert log_entries[2].endswith(b"'astropy.config.tests.test_logging', 'INFO', 'Information message'\n")

    if len(log_entries) >= 4:
        assert log_entries[3].endswith(b"'astropy.config.tests.test_logging', 'DEBUG', 'Debug message'\n")


def test_log_to_file_level():

    log_file = NamedTemporaryFile()

    with log.log_to_file(log_file.name, filter_level='ERROR'):
        log.error("Error message")
        log.warn("Warning message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 1 and log_entries[0].endswith(b"'ERROR', 'Error message'\n")


def test_log_to_file_origin1():

    log_file = NamedTemporaryFile()

    with log.log_to_file(log_file.name, filter_origin='astropy.config.tests'):
        log.error("Error message")
        log.warn("Warning message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 2


def test_log_to_file_origin2():

    log_file = NamedTemporaryFile()

    with log.log_to_file(log_file.name, filter_origin='astropy.wcs'):
        log.error("Error message")
        log.warn("Warning message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 0
