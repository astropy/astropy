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


import logging


# The location of the online documentation for astropy
# This location will normally point to the current released version of astropy
if 'dev' in __version__:
    online_docs_root = 'http://docs.astropy.org/en/latest/'
else:
    online_docs_root = 'http://docs.astropy.org/en/{0}/'.format(__version__)


# set up the test command
def _get_test_runner():
    from .tests.helper import TestRunner
    return TestRunner(__path__[0])


def test(package=None, test_path=None, args=None, plugins=None,
         verbose=False, pastebin=None, remote_data=False, pep8=False,
         pdb=False, coverage=False, open_files=False):
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

    coverage : bool, optional
        Generate a test coverage report.  The result will be placed in
        the directory htmlcov.

    open_files : bool, optional
        Fail when any tests leave files open.  Off by default, because
        this adds extra run time to the test suite.  Works only on
        platforms with a working `lsof` command.

    See Also
    --------
    pytest.main : py.test function wrapped by `run_tests`.

    """
    test_runner = _get_test_runner()
    return test_runner.run_tests(
        package=package, test_path=test_path, args=args,
        plugins=plugins, verbose=verbose, pastebin=pastebin,
        remote_data=remote_data, pep8=pep8, pdb=pdb,
        coverage=coverage, open_files=open_files)


# Use the root logger as a dummy log before initilizing Astropy's logger
log = logging.getLogger()

# if we are *not* in setup mode, import the logger and possibly populate the
# configuration file with the defaults
if not _ASTROPY_SETUP_:
    from .logger import _init_log
    from . import config

    import os
    import sys
    from warnings import warn

    try:
        from .utils import _compiler
    except ImportError:
        if os.path.exists('setup.py'):
            log.error('You appear to be trying to import astropy from within '
                      'a source checkout; please run `./setup.py develop` or '
                      '`./setup.py build_ext --inplace` first so that '
                      'extension modules can be compiled and made importable.')
            sys.exit(1)
        else:
            # Outright broken installation; don't be nice.
            raise



    # add these here so we only need to cleanup the namespace at the end
    config_dir = None

    log = _init_log()

    config_dir = os.path.dirname(__file__)
    try:
        config.configuration.update_default_config(__package__, config_dir)
    except config.configuration.ConfigurationDefaultMissingError as e:
        wmsg = (e.args[0] + " Cannot install default profile. If you are "
                "importing from source, this is expected.")
        warn(config.configuration.ConfigurationDefaultMissingWarning(wmsg))
        del e

    del _init_log, os, warn, config_dir  # clean up namespace
