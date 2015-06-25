# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['__version__', '__githash__', 'test']

# this indicates whether or not we are in the package's setup.py
try:
    _ASTROPY_SETUP_
except NameError:
    from sys import version_info
    if version_info[0] >= 3:
        import builtins
    else:
        import __builtin__ as builtins
    builtins._ASTROPY_SETUP_ = False

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''
try:
    from .version import githash as __githash__
except ImportError:
    __githash__ = ''

# set up the test command
def _get_test_runner():
    import os
    from astropy.tests.helper import TestRunner
    return TestRunner(os.path.dirname(__file__))

def test(package=None, test_path=None, args=None, plugins=None, verbose=False,
         pastebin=None, remote_data=False, pep8=False, pdb=False,
         coverage=False, open_files=False, **kwargs):
    """
    Run the tests using `py.test <http://pytest.org/latest>`__. A proper set
    of arguments is constructed and passed to `pytest.main`_.

    .. _py.test: http://pytest.org/latest/
    .. _pytest.main: http://pytest.org/latest/builtin.html#pytest.main

    Parameters
    ----------
    package : str, optional
        The name of a specific package to test, e.g. 'io.fits' or 'utils'.
        If nothing is specified all default tests are run.

    test_path : str, optional
        Specify location to test by path. May be a single file or
        directory. Must be specified absolutely or relative to the
        calling directory.

    args : str, optional
        Additional arguments to be passed to pytest.main_ in the ``args``
        keyword argument.

    plugins : list, optional
        Plugins to be passed to pytest.main_ in the ``plugins`` keyword
        argument.

    verbose : bool, optional
        Convenience option to turn on verbose output from py.test_. Passing
        True is the same as specifying ``'-v'`` in ``args``.

    pastebin : {'failed','all',None}, optional
        Convenience option for turning on py.test_ pastebin output. Set to
        ``'failed'`` to upload info for failed tests, or ``'all'`` to upload
        info for all tests.

    remote_data : bool, optional
        Controls whether to run tests marked with @remote_data. These
        tests use online data and are not run by default. Set to True to
        run these tests.

    pep8 : bool, optional
        Turn on PEP8 checking via the `pytest-pep8 plugin
        <http://pypi.python.org/pypi/pytest-pep8>`_ and disable normal
        tests. Same as specifying ``'--pep8 -k pep8'`` in ``args``.

    pdb : bool, optional
        Turn on PDB post-mortem analysis for failing tests. Same as
        specifying ``'--pdb'`` in ``args``.

    coverage : bool, optional
        Generate a test coverage report.  The result will be placed in
        the directory htmlcov.

    open_files : bool, optional
        Fail when any tests leave files open.  Off by default, because
        this adds extra run time to the test suite.  Requires the
        `psutil <https://pypi.python.org/pypi/psutil>`_ package.

    parallel : int, optional
        When provided, run the tests in parallel on the specified
        number of CPUs.  If parallel is negative, it will use the all
        the cores on the machine.  Requires the
        `pytest-xdist <https://pypi.python.org/pypi/pytest-xdist>`_ plugin
        installed. Only available when using Astropy 0.3 or later.

    kwargs
        Any additional keywords passed into this function will be passed
        on to the astropy test runner.  This allows use of test-related
        functionality implemented in later versions of astropy without
        explicitly updating the package template.

    """
    test_runner = _get_test_runner()
    return test_runner.run_tests(
        package=package, test_path=test_path, args=args,
        plugins=plugins, verbose=verbose, pastebin=pastebin,
        remote_data=remote_data, pep8=pep8, pdb=pdb,
        coverage=coverage, open_files=open_files, **kwargs)

if not _ASTROPY_SETUP_:
    import os
    from warnings import warn
    from astropy import config

    # add these here so we only need to cleanup the namespace at the end
    config_dir = None

    if not os.environ.get('ASTROPY_SKIP_CONFIG_UPDATE', False):
        config_dir = os.path.dirname(__file__)
        config_template = os.path.join(config_dir, __package__ + ".cfg")
        if os.path.isfile(config_template):
            try:
                config.configuration.update_default_config(
                    __package__, config_dir, version=__version__)
            except TypeError as orig_error:
                try:
                    config.configuration.update_default_config(
                        __package__, config_dir)
                except config.configuration.ConfigurationDefaultMissingError as e:
                    wmsg = (e.args[0] + " Cannot install default profile. If you are "
                            "importing from source, this is expected.")
                    warn(config.configuration.ConfigurationDefaultMissingWarning(wmsg))
                    del e
                except:
                    raise orig_error
