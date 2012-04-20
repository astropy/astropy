import sys
import warnings

import pytest

from .. import logger

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
    with pytest.raises(Exception) as e:
        logger.disable_warnings_logging()
    assert e.value.args[0] == 'Warnings logging has not been enabled'


def test_warnings_logging_enable_twice():
    logger.enable_warnings_logging()
    with pytest.raises(Exception) as e:
        logger.enable_warnings_logging()
    assert e.value.args[0] == 'Warnings logging has already been enabled'


def test_warnings_logging_overridden():
    logger.enable_warnings_logging()
    warnings.showwarning = lambda: None
    with pytest.raises(Exception) as e:
        logger.disable_warnings_logging()
    assert e.value.args[0] == 'Cannot disable warnings logging: warnings.showwarning was not set by this logger, or has been overridden'


def test_exception_logging_disable_no_enable():
    with pytest.raises(Exception) as e:
        logger.disable_exception_logging()
    assert e.value.args[0] == 'Exception logging has not been enabled'


def test_exception_logging_enable_twice():
    logger.enable_exception_logging()
    with pytest.raises(Exception) as e:
        logger.enable_exception_logging()
    assert e.value.args[0] == 'Exception logging has already been enabled'


def test_exception_logging_overridden():
    logger.enable_exception_logging()
    sys.excepthook = lambda: None
    with pytest.raises(Exception) as e:
        logger.disable_exception_logging()
    assert e.value.args[0] == 'Cannot disable exception logging: sys.excepthook was not set by this logger, or has been overridden'


@pytest.mark.parametrize(('level'), [None, 'DEBUG', 'INFO', 'WARN', 'ERROR'])
def test_catch_to_list_default(level):

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
