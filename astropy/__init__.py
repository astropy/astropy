# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them.
"""

from __future__ import absolute_import

import sys
import os
from warnings import warn

if sys.version_info[:2] < (2, 7):
    warn("Astropy does not support Python 2.6 (in v1.2 and later)")


def _is_astropy_source(path=None):
    """
    Returns whether the source for this module is directly in an astropy
    source distribution or checkout.
    """

    # If this __init__.py file is in ./astropy/ then import is within a source
    # dir .astropy-root is a file distributed with the source, but that should
    # not installed
    if path is None:
        path = os.path.join(os.path.dirname(__file__), os.pardir)
    elif os.path.isfile(path):
        path = os.path.dirname(path)

    source_dir = os.path.abspath(path)
    return os.path.exists(os.path.join(source_dir, '.astropy-root'))


def _is_astropy_setup():
    """
    Returns whether we are currently being imported in the context of running
    Astropy's setup.py.
    """

    main_mod = sys.modules.get('__main__')
    if not main_mod:
        return False

    return (getattr(main_mod, '__file__', False) and
            os.path.basename(main_mod.__file__).rstrip('co') == 'setup.py' and
            _is_astropy_source(main_mod.__file__))


# this indicates whether or not we are in astropy's setup.py
try:
    _ASTROPY_SETUP_
except NameError:
    from sys import version_info
    if version_info[0] >= 3:
        import builtins
    else:
        import __builtin__ as builtins

    # This will set the _ASTROPY_SETUP_ to True by default if
    # we are running Astropy's setup.py
    builtins._ASTROPY_SETUP_ = _is_astropy_setup()


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


__minimum_numpy_version__ = '1.9.0'


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
        from .utils import minversion
        requirement_met = minversion(numpy, __minimum_numpy_version__)

    if not requirement_met:
        msg = ("Numpy version {0} or later must be installed to use "
               "Astropy".format(__minimum_numpy_version__))
        raise ImportError(msg)

    return numpy


if not _ASTROPY_SETUP_:
    _check_numpy()


from . import config as _config


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

# Create the test() function
from .tests.runner import TestRunner
test = TestRunner.make_test_runner_in(__path__[0])


# if we are *not* in setup mode, import the logger and possibly populate the
# configuration file with the defaults
def _initialize_astropy():
    from . import config

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
        if _is_astropy_source():
            log.warning('You appear to be trying to import astropy from '
                        'within a source checkout without building the '
                        'extension modules first.  Attempting to (re)build '
                        'extension modules:')

            try:
                _rebuild_extensions()
            except BaseException as exc:
                _rollback_import(
                    'An error occurred while attempting to rebuild the '
                    'extension modules.  Please try manually running '
                    '`./setup.py develop` or `./setup.py build_ext '
                    '--inplace` to see what the issue was.  Extension '
                    'modules must be successfully compiled and importable '
                    'in order to import astropy.')
                # Reraise the Exception only in case it wasn't an Exception,
                # for example if a "SystemExit" or "KeyboardInterrupt" was
                # invoked.
                if not isinstance(exc, Exception):
                    raise

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


def _rebuild_extensions():
    global __version__
    global __githash__

    import subprocess
    import time

    from .utils.console import Spinner
    from .extern.six import next

    devnull = open(os.devnull, 'w')
    old_cwd = os.getcwd()
    os.chdir(os.path.join(os.path.dirname(__file__), os.pardir))
    try:
        sp = subprocess.Popen([sys.executable, 'setup.py', 'build_ext',
                               '--inplace'], stdout=devnull,
                               stderr=devnull)
        with Spinner('Rebuilding extension modules') as spinner:
            while sp.poll() is None:
                next(spinner)
                time.sleep(0.05)
    finally:
        os.chdir(old_cwd)
        devnull.close()

    if sp.returncode != 0:
        raise OSError('Running setup.py build_ext --inplace failed '
                      'with error code {0}: try rerunning this command '
                      'manually to check what the error was.'.format(
                          sp.returncode))

    # Try re-loading module-level globals from the astropy.version module,
    # which may not have existed before this function ran
    try:
        from .version import version as __version__
    except ImportError:
        pass

    try:
        from .version import githash as __githash__
    except ImportError:
        pass


# Set the bibtex entry to the article referenced in CITATION
def _get_bibtex():
    import re

    citation_file = os.path.join(os.path.dirname(__file__), 'CITATION')

    with open(citation_file, 'r') as citation:
        refs = re.findall(r'\{[^()]*\}', citation.read())
        if len(refs) == 0: return ''
        bibtexreference = "@ARTICLE{0}".format(refs[0])
    return bibtexreference


__citation__ = __bibtex__ = _get_bibtex()

import logging

# Use the root logger as a dummy log before initilizing Astropy's logger
log = logging.getLogger()


if not _ASTROPY_SETUP_:
    from .logger import _init_log, _teardown_log

    log = _init_log()

    _initialize_astropy()

    from .utils.misc import find_api_page


def online_help(query):
    """
    Search the online Astropy documentation for the given query.
    Opens the results in the default web browser.  Requires an active
    Internet connection.

    Parameters
    ----------
    query : str
        The search query.
    """
    from .extern.six.moves.urllib.parse import urlencode
    import webbrowser

    version = __version__
    if 'dev' in version:
        version = 'latest'
    else:
        version = 'v' + version

    url = 'http://docs.astropy.org/en/{0}/search.html?{1}'.format(
        version, urlencode({'q': query}))

    webbrowser.open(url)


__dir_inc__ = ['__version__', '__githash__', '__minimum_numpy_version__',
               '__bibtex__', 'test', 'log', 'find_api_page', 'online_help',
               'online_docs_root', 'conf']


from types import ModuleType as __module_type__
# Clean up top-level namespace--delete everything that isn't in __dir_inc__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not ((varname.startswith('__') and varname.endswith('__')) or
            varname in __dir_inc__ or
            (varname[0] != '_' and
                isinstance(locals()[varname], __module_type__) and
                locals()[varname].__name__.startswith(__name__ + '.'))):
        # The last clause in the the above disjunction deserves explanation:
        # When using relative imports like ``from .. import config``, the
        # ``config`` variable is automatically created in the namespace of
        # whatever module ``..`` resolves to (in this case astropy).  This
        # happens a few times just in the module setup above.  This allows
        # the cleanup to keep any public submodules of the astropy package
        del locals()[varname]

del varname, __module_type__
