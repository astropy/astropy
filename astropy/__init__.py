# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them.
"""

# this indicates whether or not we are in astropy's setup.py
try:
    _ASTROPY_SETUP_
except NameError:
    from sys import version_info
    if version_info[0] >= 3:
        import builtins
    else:
        import __builtin__ as builtins
    builtins._ASTROPY_SETUP_ = False
    del version_info
    del builtins

try:
    from .version import version as __version__
except ImportError:
    # TODO: Issue a warning using the logging framework
    __version__ = ''
try:
    from .version import githash as __githash__
except ImportError:
    # TODO: Issue a warning using the logging framework
    __githash__ = ''


__minimum_numpy_version__ = '1.5.1'


# The location of the online documentation for astropy
# This location will normally point to the current released version of astropy
if 'dev' in __version__:
    online_docs_root = 'http://docs.astropy.org/en/latest/'
else:
    online_docs_root = 'http://docs.astropy.org/en/{0}/'.format(__version__)


def _check_numpy():
    """
    Check that Numpy is installed and it is of the minimum version we
    require.
    """
    # Note: We could have used distutils.version for this comparison,
    # but it seems like overkill to import distutils at runtime.
    requirement_met = False

    try:
        import numpy
    except ImportError:
        pass
    else:
        major, minor, rest = numpy.__version__.split(".", 2)
        rmajor, rminor, rest = __minimum_numpy_version__.split(".", 2)
        requirement_met = ((int(major), int(minor)) >=
                           (int(rmajor), int(rminor)))

    if not requirement_met:
        msg = ("numpy version {0} or later must be installed to use "
               "astropy".format(
                   __minimum_numpy_version__))
        raise ImportError(msg)

    return numpy


if not _ASTROPY_SETUP_:
    _check_numpy()


from . import config as _config
import sys


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy`.
    """

    unicode_output = _config.ConfigItem(
        False,
        'When True, use Unicode characters when outputting values, and '
        'displaying widgets at the console.')
    use_color = _config.ConfigItem(
        sys.platform != 'win32',
        'When True, use ANSI color escape sequences when writing to the console.',
        aliases=['astropy.utils.console.USE_COLOR', 'astropy.logger.USE_COLOR'])
    max_lines = _config.ConfigItem(
        None,
        description='Maximum number of lines in the display of pretty-printed '
        'objects. If not provided, try to determine automatically from the '
        'terminal size.  Negative numbers mean no limit.',
        cfgtype='integer(default=None)',
        aliases=['astropy.table.pprint.max_lines'])
    max_width = _config.ConfigItem(
        None,
        description='Maximum number of characters per line in the display of '
        'pretty-printed objects.  If not provided, try to determine '
        'automatically from the terminal size. Negative numbers mean no '
        'limit.',
        cfgtype='integer(default=None)',
        aliases=['astropy.table.pprint.max_width'])

conf = Conf()


UNICODE_OUTPUT = _config.ConfigAlias(
    '0.4', 'UNICODE_OUTPUT', 'unicode_output')


del sys


# set up the test command
def _get_test_runner():
    from .tests.helper import TestRunner
    return TestRunner(__path__[0])


def test(package=None, test_path=None, args=None, plugins=None,
         verbose=False, pastebin=None, remote_data=False, pep8=False,
         pdb=False, open_files=False, parallel=0, docs_path=None,
         skip_docs=False, repeat=None):
    """
    Run Astropy tests using py.test. A proper set of arguments is
    constructed and passed to `pytest.main`.

    Parameters
    ----------
    package : str, optional
        The name of a specific package to test, e.g. 'io.fits' or 'utils'.
        If nothing is specified all default Astropy tests are run.

    test_path : str, optional
        Specify location to test by path. May be a single file or
        directory. Must be specified absolutely or relative to the
        calling directory.

    args : str, optional
        Additional arguments to be passed to `pytest.main` in the `args`
        keyword argument.

    plugins : list, optional
        Plugins to be passed to `pytest.main` in the `plugins` keyword
        argument.

    verbose : bool, optional
        Convenience option to turn on verbose output from py.test. Passing
        True is the same as specifying `-v` in `args`.

    pastebin : {'failed','all',None}, optional
        Convenience option for turning on py.test pastebin output. Set to
        'failed' to upload info for failed tests, or 'all' to upload info
        for all tests.

    remote_data : bool, optional
        Controls whether to run tests marked with @remote_data. These
        tests use online data and are not run by default. Set to True to
        run these tests.

    pep8 : bool, optional
        Turn on PEP8 checking via the pytest-pep8 plugin and disable normal
        tests. Same as specifying `--pep8 -k pep8` in `args`.

    pdb : bool, optional
        Turn on PDB post-mortem analysis for failing tests. Same as
        specifying `--pdb` in `args`.

    open_files : bool, optional
        Fail when any tests leave files open.  Off by default, because
        this adds extra run time to the test suite.  Works only on
        platforms with a working `lsof` command.

    parallel : int, optional
        When provided, run the tests in parallel on the specified
        number of CPUs.  If parallel is negative, it will use the all
        the cores on the machine.  Requires the `pytest-xdist` plugin
        is installed.

    docs_path : str, optional
        The path to the documentation .rst files.

    skip_docs : bool, optional
        When `True`, skips running the doctests in the .rst files.

    repeat : int, optional
        If set, specifies how many times each test should be run. This is
        useful for diagnosing sporadic failures.

    See Also
    --------
    pytest.main : py.test function wrapped by `run_tests`.

    """
    test_runner = _get_test_runner()
    return test_runner.run_tests(
        package=package, test_path=test_path, args=args,
        plugins=plugins, verbose=verbose, pastebin=pastebin,
        remote_data=remote_data, pep8=pep8, pdb=pdb,
        open_files=open_files, parallel=parallel, docs_path=docs_path,
        skip_docs=skip_docs, repeat=repeat)


# if we are *not* in setup mode, import the logger and possibly populate the
# configuration file with the defaults
def _initialize_astropy():
    from . import config

    import os
    import sys
    from warnings import warn

    # If this __init__.py file is in ./astropy/ then import is within a source dir
    is_astropy_source_dir = (os.path.abspath(os.path.dirname(__file__)) ==
                             os.path.abspath('astropy') and os.path.exists('setup.py'))

    def _rollback_import(message):
        log.error(message)
        # Now disable exception logging to avoid an annoying error in the
        # exception logger before we raise the import error:
        _teardown_log()

        # Roll back any astropy sub-modules that have been imported thus
        # far

        for key in list(sys.modules):
            if key.startswith('astropy.'):
                del sys.modules[key]
        raise ImportError('astropy')

    try:
        from .utils import _compiler
    except ImportError:
        if is_astropy_source_dir:
            _rollback_import(
                'You appear to be trying to import astropy from within a '
                'source checkout; please run `./setup.py develop` or '
                '`./setup.py build_ext --inplace` first so that extension '
                'modules can be compiled and made importable.')
        else:
            # Outright broken installation; don't be nice.
            raise

    # add these here so we only need to cleanup the namespace at the end
    config_dir = os.path.dirname(__file__)

    try:
        config.configuration.update_default_config(__package__, config_dir)
    except config.configuration.ConfigurationDefaultMissingError as e:
        wmsg = (e.args[0] + " Cannot install default profile. If you are "
                "importing from source, this is expected.")
        warn(config.configuration.ConfigurationDefaultMissingWarning(wmsg))


import logging

# Use the root logger as a dummy log before initilizing Astropy's logger
log = logging.getLogger()


if not _ASTROPY_SETUP_:
    from .logger import _init_log, _teardown_log

    log = _init_log()

    _initialize_astropy()

    from .utils.misc import find_api_page
