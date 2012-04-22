import sys
import warnings
from tempfile import NamedTemporaryFile

import pytest

from .. import logger, LoggingError


# Save original values of hooks
_excepthook = sys.excepthook
_showwarning = warnings.showwarning


def setup_function(function):

    # Reset hooks to original values
    sys.excepthook = _excepthook
    warnings.showwarning = _showwarning

    # Reset internal original hooks
    logger._showwarning_orig = None
    logger._excepthook_orig = None

    # Set up the logger
    logger.set_defaults()


def teardown_module(function):

    # Ensure that hooks are restored to original values
    sys.excepthook = _excepthook
    warnings.showwarning = _showwarning


def test_warnings_logging_disable_no_enable():
    with pytest.raises(LoggingError) as e:
        logger.disable_warnings_logging()
    assert e.value.args[0] == 'Warnings logging has not been enabled'


def test_warnings_logging_enable_twice():
    logger.enable_warnings_logging()
    with pytest.raises(LoggingError) as e:
        logger.enable_warnings_logging()
    assert e.value.args[0] == 'Warnings logging has already been enabled'


def test_warnings_logging_overridden():
    logger.enable_warnings_logging()
    warnings.showwarning = lambda: None
    with pytest.raises(LoggingError) as e:
        logger.disable_warnings_logging()
    assert e.value.args[0] == 'Cannot disable warnings logging: warnings.showwarning was not set by this logger, or has been overridden'


@pytest.mark.xfail  # origin is not set correctly
def test_warnings_logging():

    # Without warnings logging
    with warnings.catch_warnings(record=True) as warn_list:
        with logger.log_to_list() as log_list:
            warnings.warn("This is a warning")
    assert len(log_list) == 0
    assert len(warn_list) == 1
    assert warn_list[0].message.args[0] == "This is a warning"

    # With warnings logging
    with warnings.catch_warnings(record=True) as warn_list:
        logger.enable_warnings_logging()
        with logger.log_to_list() as log_list:
            warnings.warn("This is a warning")
    assert len(log_list) == 1
    assert len(warn_list) == 0
    assert log_list[0].levelname == 'WARNING'
    assert log_list[0].message == 'This is a warning'
    assert log_list[0].origin == 'astropy.config.tests.test_logging'

    # Without warnings logging
    logger.disable_warnings_logging()
    with warnings.catch_warnings(record=True) as warn_list:
        with logger.log_to_list() as log_list:
            raise Exception("This is a warning")
    assert len(log_list) == 0
    assert len(warn_list) == 1
    assert warn_list[0].message.args[0] == "This is a warning"


def test_exception_logging_disable_no_enable():
    with pytest.raises(LoggingError) as e:
        logger.disable_exception_logging()
    assert e.value.args[0] == 'Exception logging has not been enabled'


def test_exception_logging_enable_twice():
    logger.enable_exception_logging()
    with pytest.raises(LoggingError) as e:
        logger.enable_exception_logging()
    assert e.value.args[0] == 'Exception logging has already been enabled'


def test_exception_logging_overridden():
    logger.enable_exception_logging()
    sys.excepthook = lambda: None
    with pytest.raises(LoggingError) as e:
        logger.disable_exception_logging()
    assert e.value.args[0] == 'Cannot disable exception logging: sys.excepthook was not set by this logger, or has been overridden'


@pytest.mark.xfail  # pytest.raises interferes with excepthook
def test_exception_logging():

    # Without exception logging
    with pytest.raises(Exception) as e:
        with logger.log_to_list() as log_list:
            raise Exception("This is an Exception")
    assert len(log_list) == 0
    assert e.value.args[0] == "This is an Exception"

    # With exception logging
    with pytest.raises(LoggingError) as e:
        logger.enable_exception_logging()
        with logger.log_to_list() as log_list:
            raise Exception("This is an Exception")
    assert len(log_list) == 1
    assert log_list[0].levelname == 'ERROR'
    assert log_list[0].message == 'This is an Exception'
    assert log_list[0].origin == 'astropy.config.tests.test_logging'
    assert e.value.args[0] == "This is an Exception"

    # Without exception logging
    logger.disable_exception_logging()
    with pytest.raises(LoggingError) as e:
        with logger.log_to_list() as log_list:
            raise Exception("This is an Exception")
    assert len(log_list) == 0
    assert e.value.args[0] == "This is an Exception"


@pytest.mark.parametrize(('level'), [None, 'DEBUG', 'INFO', 'WARN', 'ERROR'])
def test_log_to_list(level):

    if level is not None:
        logger.setLevel(level)

    with logger.log_to_list() as log_list:
        logger.error("Error message")
        logger.warn("Warning message")
        logger.info("Information message")
        logger.debug("Debug message")

    # Check list length
    if level == 'DEBUG':
        assert len(log_list) == 4
    elif level == 'INFO':
        assert len(log_list) == 3
    elif level is None or level == 'WARN':
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

    with logger.log_to_list(filter_level='ERROR') as log_list:
        logger.error("Error message")
        logger.warn("Warning message")

    assert len(log_list) == 1 and log_list[0].levelname == 'ERROR'


def test_log_to_list_origin1():

    with logger.log_to_list(filter_origin='astropy.config.tests') as log_list:
        logger.error("Error message")
        logger.warn("Warning message")

    assert len(log_list) == 2


def test_log_to_list_origin2():

    with logger.log_to_list(filter_origin='astropy.wcs') as log_list:
        logger.error("Error message")
        logger.warn("Warning message")

    assert len(log_list) == 0


@pytest.mark.parametrize(('level'), [None, 'DEBUG', 'INFO', 'WARN', 'ERROR'])
def test_log_to_file(level):

    log_file = NamedTemporaryFile()

    if level is not None:
        logger.setLevel(level)

    with logger.log_to_file(log_file.name):
        logger.error("Error message")
        logger.warn("Warning message")
        logger.info("Information message")
        logger.debug("Debug message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    # Check list length
    if level == 'DEBUG':
        assert len(log_entries) == 4
    elif level == 'INFO':
        assert len(log_entries) == 3
    elif level is None or level == 'WARN':
        assert len(log_entries) == 2
    elif level == 'ERROR':
        assert len(log_entries) == 1

    # Check list content

    assert log_entries[0].endswith('astropy.config.tests.test_logging, ERROR, Error message\n')

    if len(log_entries) >= 2:
        assert log_entries[1].endswith('astropy.config.tests.test_logging, WARNING, Warning message\n')

    if len(log_entries) >= 3:
        assert log_entries[2].endswith('astropy.config.tests.test_logging, INFO, Information message\n')

    if len(log_entries) >= 4:
        assert log_entries[3].endswith('astropy.config.tests.test_logging, DEBUG, Debug message\n')


def test_log_to_file_level():

    log_file = NamedTemporaryFile()

    with logger.log_to_file(log_file.name, filter_level='ERROR'):
        logger.error("Error message")
        logger.warn("Warning message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 1 and log_entries[0].endswith('ERROR, Error message\n')


def test_log_to_file_origin1():

    log_file = NamedTemporaryFile()

    with logger.log_to_file(log_file.name, filter_origin='astropy.config.tests'):
        logger.error("Error message")
        logger.warn("Warning message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 2


def test_log_to_file_origin2():

    log_file = NamedTemporaryFile()

    with logger.log_to_file(log_file.name, filter_origin='astropy.wcs'):
        logger.error("Error message")
        logger.warn("Warning message")

    log_file.seek(0)
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 0
