# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them.
"""

# Prior to Astropy 3.2, astropy was imported during setup.py commands. If we are
# in setup mode, then astropy-helpers defines an _ASTROPY_SETUP_ variable, which
# we used to use to conditionally import C extensions for example. However, the
# behavior of importing the package during the setup process is not good
# practice and we therefore now explicitly prevent the package from being
# imported in that case to prevent any regressions. We use _ASTROPY_CORE_SETUP_
# (defined in setup.py) rather than _ASTROPY_SETUP_ since the latter is also
# set up for affiliated packages, and those need to be able to import the
# (installed) core package during e.g. python setup.py test.
try:
    _ASTROPY_CORE_SETUP_
except NameError:
    pass
else:
    raise RuntimeError("The astropy package cannot be imported during setup")

import sys
import os
from warnings import warn

__minimum_python_version__ = '3.6'
__minimum_numpy_version__ = '1.16.0'
__minimum_scipy_version__ = '0.18'
# ASDF is an optional dependency, but this is the minimum version that is
# compatible with Astropy when it is installed.
__minimum_asdf_version__ = '2.5.2'
# PyYAML is an optional dependency, but this is the minimum version that is
# advertised to be supported.
__minimum_yaml_version__ = '3.12'


class UnsupportedPythonError(Exception):
    pass


# This is the same check as the one at the top of setup.py
if sys.version_info < tuple(int(val) for val in __minimum_python_version__.split('.')):
    raise UnsupportedPythonError(f"Astropy does not support Python < {__minimum_python_version__}")


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


# The location of the online documentation for astropy
# This location will normally point to the current released version of astropy
if 'dev' in __version__:
    online_docs_root = 'https://docs.astropy.org/en/latest/'
else:
    online_docs_root = f'https://docs.astropy.org/en/{__version__}/'


def _check_numpy():
    """
    Check that Numpy is installed and it is of the minimum version we
    require.
    """
    # Note: We could have used distutils.version for this comparison,
    # but it seems like overkill to import distutils at runtime.
    requirement_met = False
    import_fail = ''
    try:
        import numpy
    except ImportError:
        import_fail = 'Numpy is not installed.'
    else:
        from .utils import minversion
        requirement_met = minversion(numpy, __minimum_numpy_version__)

    if not requirement_met:
        msg = (f"Numpy version {__minimum_numpy_version__} or later must "
               f"be installed to use Astropy. {import_fail}")
        raise ImportError(msg)

    return numpy


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


# Define a base ScienceState for configuring constants and units
from .utils.state import ScienceState


class base_constants_version(ScienceState):
    """
    Base class for the real version-setters below
    """
    _value = 'test'

    _versions = dict(test='test')

    @classmethod
    def validate(cls, value):
        if value not in cls._versions:
            raise ValueError('Must be one of {}'
                             .format(list(cls._versions.keys())))
        return cls._versions[value]

    @classmethod
    def set(cls, value):
        """
        Set the current constants value.
        """
        import sys
        if 'astropy.units' in sys.modules:
            raise RuntimeError('astropy.units is already imported')
        if 'astropy.constants' in sys.modules:
            raise RuntimeError('astropy.constants is already imported')

        class _Context:
            def __init__(self, parent, value):
                self._value = value
                self._parent = parent

            def __enter__(self):
                pass

            def __exit__(self, type, value, tb):
                self._parent._value = self._value

            def __repr__(self):
                return ('<ScienceState {}: {!r}>'
                        .format(self._parent.__name__, self._parent._value))

        ctx = _Context(cls, cls._value)
        value = cls.validate(value)
        cls._value = value
        return ctx


class physical_constants(base_constants_version):
    """
    The version of physical constants to use
    """
    # Maintainers: update when new constants are added
    _value = 'codata2018'

    _versions = dict(codata2018='codata2018', codata2014='codata2014',
                     codata2010='codata2010', astropyconst40='codata2018',
                     astropyconst20='codata2014', astropyconst13='codata2010')


class astronomical_constants(base_constants_version):
    """
    The version of astronomical constants to use
    """
    # Maintainers: update when new constants are added
    _value = 'iau2015'

    _versions = dict(iau2015='iau2015', iau2012='iau2012',
                     astropyconst40='iau2015', astropyconst20='iau2015',
                     astropyconst13='iau2012')


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
                      'with error code {}: try rerunning this command '
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


# Set the bibtex entry to the article referenced in CITATION.
def _get_bibtex():
    citation_file = os.path.join(os.path.dirname(__file__), 'CITATION')

    with open(citation_file, 'r') as citation:
        refs = citation.read().split('@ARTICLE')[1:]
        if len(refs) == 0: return ''
        bibtexreference = "@ARTICLE{}".format(refs[0])
    return bibtexreference


__citation__ = __bibtex__ = _get_bibtex()

import logging

# Use the root logger as a dummy log before initilizing Astropy's logger
log = logging.getLogger()


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
    from urllib.parse import urlencode
    import webbrowser

    version = __version__
    if 'dev' in version:
        version = 'latest'
    else:
        version = 'v' + version

    url = 'https://docs.astropy.org/en/{}/search.html?{}'.format(
        version, urlencode({'q': query}))

    webbrowser.open(url)


__dir_inc__ = ['__version__', '__githash__', '__minimum_numpy_version__',
               '__bibtex__', 'test', 'log', 'find_api_page', 'online_help',
               'online_docs_root', 'conf', 'physical_constants',
               'astronomical_constants']


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
