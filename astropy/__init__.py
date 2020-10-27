# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them.
"""

import sys
import os
from warnings import warn

from .version import version as __version__

__minimum_python_version__ = '3.7'
__minimum_numpy_version__ = '1.17.0'
__minimum_erfa_version__ = '1.7'
__minimum_scipy_version__ = '1.1'
# ASDF is an optional dependency, but this is the minimum version that is
# compatible with Astropy when it is installed.
__minimum_asdf_version__ = '2.6.0'
# PyYAML is an optional dependency, but this is the minimum version that is
# advertised to be supported.
__minimum_yaml_version__ = '3.13'
# Matplotlib is an optional dependency, but this is the minimum version that is
# advertised to be supported.
__minimum_matplotlib_version__ = '3.0'


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


# The location of the online documentation for astropy
# This location will normally point to the current released version of astropy
if 'dev' in __version__:
    online_docs_root = 'https://docs.astropy.org/en/latest/'
else:
    online_docs_root = f'https://docs.astropy.org/en/{__version__}/'


def _check_requirement(name, minimum_version):
    """
    Check that ``name`` is installed and it is of the ``minimum_version`` we
    require.
    """
    # Note: We could have used distutils.version for this comparison,
    # but it seems like overkill to import distutils at runtime, and
    # our own utils.introspection.minversion indirectly needs requirements.
    requirement_met = False
    import_fail = ''
    try:
        module = __import__(name)
    except ImportError:
        import_fail = f'{name} is not installed.'
    else:
        version = getattr(module, '__version__')
        requirement_met = version.split('.') >= minimum_version.split('.')

    if not requirement_met:
        msg = (f"{name} version {minimum_version} or later must "
               f"be installed to use Astropy. {import_fail}")
        raise ImportError(msg)

    return module


_check_requirement('numpy', __minimum_numpy_version__)
_check_requirement('erfa', __minimum_erfa_version__)


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

        return super().set(cls, value)


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

    try:
        from .utils import _compiler
    except ImportError:
        if _is_astropy_source():
            raise ImportError('You appear to be trying to import astropy from '
                              'within a source checkout or from an editable '
                              'installation without building the extension '
                              'modules first. Either run:\n\n'
                              '  pip install -e .\n\nor\n\n'
                              '  python setup.py build_ext --inplace\n\n'
                              'to make sure the extension modules are built ')
        else:
            # Outright broken installation, just raise standard error
            raise

    # add these here so we only need to cleanup the namespace at the end
    config_dir = os.path.dirname(__file__)

    try:
        config.configuration.update_default_config(__package__, config_dir)
    except config.configuration.ConfigurationDefaultMissingError as e:
        wmsg = (e.args[0] + " Cannot install default profile. If you are "
                "importing from source, this is expected.")
        warn(config.configuration.ConfigurationDefaultMissingWarning(wmsg))


# Set the bibtex entry to the article referenced in CITATION.
def _get_bibtex():
    citation_file = os.path.join(os.path.dirname(__file__), 'CITATION')

    with open(citation_file, 'r') as citation:
        refs = citation.read().split('@ARTICLE')[1:]
        if len(refs) == 0: return ''
        bibtexreference = f'@ARTICLE{refs[0]}'
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

    url = f"https://docs.astropy.org/en/{version}/search.html?{urlencode({'q': query})}"
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
