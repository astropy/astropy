# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import sys
import warnings
import logging
import locale

import pytest

from astropy import log
from astropy.logger import LoggingError, conf
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning

# Save original values of hooks. These are not the system values, but the
# already overwritten values since the logger already gets imported before
# this file gets executed.
_excepthook = sys.__excepthook__
_showwarning = warnings.showwarning

try:
    ip = get_ipython()
except NameError:
    ip = None


def setup_function(function):

    # Reset modules to default
    importlib.reload(warnings)
    importlib.reload(sys)

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


teardown_module = setup_function


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
    with pytest.raises(LoggingError, match=r'Cannot disable warnings logging: '
                       r'warnings\.showwarning was not set by this logger, or has been overridden'):
        log.disable_warnings_logging()


def test_warnings_logging():

    # Without warnings logging
    with pytest.warns(AstropyUserWarning, match="This is a warning") as warn_list:
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning", AstropyUserWarning)
    assert len(log_list) == 0
    assert len(warn_list) == 1

    # With warnings logging
    with warnings.catch_warnings(record=True) as warn_list:
        log.enable_warnings_logging()
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning", AstropyUserWarning)
        log.disable_warnings_logging()
    assert len(log_list) == 1
    assert len(warn_list) == 0
    assert log_list[0].levelname == 'WARNING'
    assert log_list[0].message.startswith('This is a warning')
    assert log_list[0].origin == 'astropy.tests.test_logger'

    # With warnings logging (differentiate between Astropy and non-Astropy)
    with pytest.warns(UserWarning, match="This is another warning, not "
                      "from Astropy") as warn_list:
        log.enable_warnings_logging()
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning", AstropyUserWarning)
            warnings.warn("This is another warning, not from Astropy")
        log.disable_warnings_logging()
    assert len(log_list) == 1
    assert len(warn_list) == 1
    assert log_list[0].levelname == 'WARNING'
    assert log_list[0].message.startswith('This is a warning')
    assert log_list[0].origin == 'astropy.tests.test_logger'

    # Without warnings logging
    with pytest.warns(AstropyUserWarning, match="This is a warning") as warn_list:
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning", AstropyUserWarning)
    assert len(log_list) == 0
    assert len(warn_list) == 1


def test_warnings_logging_with_custom_class():
    class CustomAstropyWarningClass(AstropyWarning):
        pass

    # With warnings logging
    with warnings.catch_warnings(record=True) as warn_list:
        log.enable_warnings_logging()
        with log.log_to_list() as log_list:
            warnings.warn("This is a warning", CustomAstropyWarningClass)
        log.disable_warnings_logging()
    assert len(log_list) == 1
    assert len(warn_list) == 0
    assert log_list[0].levelname == 'WARNING'
    assert log_list[0].message.startswith('CustomAstropyWarningClass: This is a warning')
    assert log_list[0].origin == 'astropy.tests.test_logger'


def test_warning_logging_with_io_votable_warning():
    from astropy.io.votable.exceptions import W02, vo_warn

    with warnings.catch_warnings(record=True) as warn_list:
        log.enable_warnings_logging()
        with log.log_to_list() as log_list:
            vo_warn(W02, ('a', 'b'))
        log.disable_warnings_logging()
    assert len(log_list) == 1
    assert len(warn_list) == 0
    assert log_list[0].levelname == 'WARNING'
    x = log_list[0].message.startswith("W02: ?:?:?: W02: a attribute 'b' is "
                                       "invalid.  Must be a standard XML id")
    assert x
    assert log_list[0].origin == 'astropy.tests.test_logger'


def test_import_error_in_warning_logging():
    """
    Regression test for https://github.com/astropy/astropy/issues/2671

    This test actually puts a goofy fake module into ``sys.modules`` to test
    this problem.
    """

    class FakeModule:
        def __getattr__(self, attr):
            raise ImportError('_showwarning should ignore any exceptions '
                              'here')

    log.enable_warnings_logging()

    sys.modules['<test fake module>'] = FakeModule()
    try:
        warnings.showwarning(AstropyWarning('Regression test for #2671'),
                             AstropyWarning, '<this is only a test>', 1)
    finally:
        del sys.modules['<test fake module>']


def test_exception_logging_disable_no_enable():
    with pytest.raises(LoggingError) as e:
        log.disable_exception_logging()
    assert e.value.args[0] == 'Exception logging has not been enabled'


def test_exception_logging_enable_twice():
    log.enable_exception_logging()
    with pytest.raises(LoggingError) as e:
        log.enable_exception_logging()
    assert e.value.args[0] == 'Exception logging has already been enabled'


# You can't really override the exception handler in IPython this way, so
# this test doesn't really make sense in the IPython context.
@pytest.mark.skipif("ip is not None")
def test_exception_logging_overridden():
    log.enable_exception_logging()
    sys.excepthook = lambda etype, evalue, tb: None
    with pytest.raises(LoggingError, match='Cannot disable exception logging: '
                       'sys.excepthook was not set by this logger, or has been overridden'):
        log.disable_exception_logging()


@pytest.mark.xfail("ip is not None")
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
    assert log_list[0].message.startswith('Exception: This is an Exception')
    assert log_list[0].origin == 'astropy.tests.test_logger'

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


@pytest.mark.xfail("ip is not None")
def test_exception_logging_origin():
    # The point here is to get an exception raised from another location
    # and make sure the error's origin is reported correctly

    from astropy.utils.collections import HomogeneousList

    l = HomogeneousList(int)  # noqa
    try:
        log.enable_exception_logging()
        with log.log_to_list() as log_list:
            l.append('foo')
    except TypeError as exc:
        sys.excepthook(*sys.exc_info())
        assert exc.args[0].startswith(
            "homogeneous list must contain only objects of type ")
    else:
        assert False
    assert len(log_list) == 1
    assert log_list[0].levelname == 'ERROR'
    assert log_list[0].message.startswith(
        "TypeError: homogeneous list must contain only objects of type ")
    assert log_list[0].origin == 'astropy.utils.collections'


@pytest.mark.skip(reason="Infinite recursion on Python 3.5+, probably a real issue")
# @pytest.mark.xfail("ip is not None")
def test_exception_logging_argless_exception():
    """
    Regression test for a crash that occurred on Python 3 when logging an
    exception that was instantiated with no arguments (no message, etc.)

    Regression test for https://github.com/astropy/astropy/pull/4056
    """

    try:
        log.enable_exception_logging()
        with log.log_to_list() as log_list:
            raise Exception()
    except Exception:
        sys.excepthook(*sys.exc_info())
    else:
        assert False  # exception should have been raised
    assert len(log_list) == 1
    assert log_list[0].levelname == 'ERROR'
    assert log_list[0].message == 'Exception [astropy.tests.test_logger]'
    assert log_list[0].origin == 'astropy.tests.test_logger'


@pytest.mark.parametrize(('level'), [None, 'DEBUG', 'INFO', 'WARN', 'ERROR'])
def test_log_to_list(level):

    orig_level = log.level

    try:
        if level is not None:
            log.setLevel(level)

        with log.log_to_list() as log_list:
            log.error("Error message")
            log.warning("Warning message")
            log.info("Information message")
            log.debug("Debug message")
    finally:
        log.setLevel(orig_level)

    if level is None:
        # The log level *should* be set to whatever it was in the config
        level = conf.log_level

    # Check list length
    if level == 'DEBUG':
        assert len(log_list) == 4
    elif level == 'INFO':
        assert len(log_list) == 3
    elif level == 'WARN':
        assert len(log_list) == 2
    elif level == 'ERROR':
        assert len(log_list) == 1

    # Check list content

    assert log_list[0].levelname == 'ERROR'
    assert log_list[0].message.startswith('Error message')
    assert log_list[0].origin == 'astropy.tests.test_logger'

    if len(log_list) >= 2:
        assert log_list[1].levelname == 'WARNING'
        assert log_list[1].message.startswith('Warning message')
        assert log_list[1].origin == 'astropy.tests.test_logger'

    if len(log_list) >= 3:
        assert log_list[2].levelname == 'INFO'
        assert log_list[2].message.startswith('Information message')
        assert log_list[2].origin == 'astropy.tests.test_logger'

    if len(log_list) >= 4:
        assert log_list[3].levelname == 'DEBUG'
        assert log_list[3].message.startswith('Debug message')
        assert log_list[3].origin == 'astropy.tests.test_logger'


def test_log_to_list_level():

    with log.log_to_list(filter_level='ERROR') as log_list:
        log.error("Error message")
        log.warning("Warning message")

    assert len(log_list) == 1 and log_list[0].levelname == 'ERROR'


def test_log_to_list_origin1():

    with log.log_to_list(filter_origin='astropy.tests') as log_list:
        log.error("Error message")
        log.warning("Warning message")

    assert len(log_list) == 2


def test_log_to_list_origin2():

    with log.log_to_list(filter_origin='astropy.wcs') as log_list:
        log.error("Error message")
        log.warning("Warning message")

    assert len(log_list) == 0


@pytest.mark.parametrize(('level'), [None, 'DEBUG', 'INFO', 'WARN', 'ERROR'])
def test_log_to_file(tmpdir, level):

    local_path = tmpdir.join('test.log')
    log_file = local_path.open('wb')
    log_path = str(local_path.realpath())
    orig_level = log.level

    try:
        if level is not None:
            log.setLevel(level)

        with log.log_to_file(log_path):
            log.error("Error message")
            log.warning("Warning message")
            log.info("Information message")
            log.debug("Debug message")

        log_file.close()
    finally:
        log.setLevel(orig_level)

    log_file = local_path.open('rb')
    log_entries = log_file.readlines()
    log_file.close()

    if level is None:
        # The log level *should* be set to whatever it was in the config
        level = conf.log_level

    # Check list length
    if level == 'DEBUG':
        assert len(log_entries) == 4
    elif level == 'INFO':
        assert len(log_entries) == 3
    elif level == 'WARN':
        assert len(log_entries) == 2
    elif level == 'ERROR':
        assert len(log_entries) == 1

    # Check list content

    assert eval(log_entries[0].strip())[-3:] == (
        'astropy.tests.test_logger', 'ERROR', 'Error message')

    if len(log_entries) >= 2:
        assert eval(log_entries[1].strip())[-3:] == (
            'astropy.tests.test_logger', 'WARNING', 'Warning message')

    if len(log_entries) >= 3:
        assert eval(log_entries[2].strip())[-3:] == (
            'astropy.tests.test_logger', 'INFO', 'Information message')

    if len(log_entries) >= 4:
        assert eval(log_entries[3].strip())[-3:] == (
            'astropy.tests.test_logger', 'DEBUG', 'Debug message')


def test_log_to_file_level(tmpdir):

    local_path = tmpdir.join('test.log')
    log_file = local_path.open('wb')
    log_path = str(local_path.realpath())

    with log.log_to_file(log_path, filter_level='ERROR'):
        log.error("Error message")
        log.warning("Warning message")

    log_file.close()

    log_file = local_path.open('rb')
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 1
    assert eval(log_entries[0].strip())[-2:] == (
        'ERROR', 'Error message')


def test_log_to_file_origin1(tmpdir):

    local_path = tmpdir.join('test.log')
    log_file = local_path.open('wb')
    log_path = str(local_path.realpath())

    with log.log_to_file(log_path, filter_origin='astropy.tests'):
        log.error("Error message")
        log.warning("Warning message")

    log_file.close()

    log_file = local_path.open('rb')
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 2


def test_log_to_file_origin2(tmpdir):

    local_path = tmpdir.join('test.log')
    log_file = local_path.open('wb')
    log_path = str(local_path.realpath())

    with log.log_to_file(log_path, filter_origin='astropy.wcs'):
        log.error("Error message")
        log.warning("Warning message")

    log_file.close()

    log_file = local_path.open('rb')
    log_entries = log_file.readlines()
    log_file.close()

    assert len(log_entries) == 0


@pytest.mark.parametrize(('encoding'), ['', 'utf-8', 'cp1252'])
def test_log_to_file_encoding(tmpdir, encoding):

    local_path = tmpdir.join('test.log')
    log_path = str(local_path.realpath())

    orig_encoding = conf.log_file_encoding

    conf.log_file_encoding = encoding

    with log.log_to_file(log_path):
        for handler in log.handlers:
            if isinstance(handler, logging.FileHandler):
                if encoding:
                    assert handler.stream.encoding == encoding
                else:
                    assert handler.stream.encoding == locale.getpreferredencoding()

    conf.log_file_encoding = orig_encoding
