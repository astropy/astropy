# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handle loading six package from system or from the bundled copy
"""

try:
    import importlib
except ImportError:
    importlib = None
    import imp
import io
import sys

from distutils.version import StrictVersion


_SIX_MIN_VERSION = StrictVersion('1.10.0')

# Update this to prevent Astropy from using its bundled copy of six
# (but only if some other version of at least _SIX_MIN_VERSION can
# be provided)
_SIX_SEARCH_PATH = ['astropy.extern.bundled.six', 'six']


def _find_module(name, path=None):
    """
    Alternative to `imp.find_module` that can also search in subpackages.
    """

    parts = name.split('.')

    for part in parts:
        if path is not None:
            path = [path]

        fh, path, descr = imp.find_module(part, path)

    return fh, path, descr


def _import_six(search_path=_SIX_SEARCH_PATH):
    for mod_name in search_path:
        if importlib is not None:
            try:
                six_mod = importlib.import_module(mod_name)
            except ImportError:
                continue
        else:
            try:
                mod_info = _find_module(mod_name)
            except ImportError:
                continue
            else:
                try:
                    # Using __name__ causes the import to effectively overwrite
                    # this shim.
                    six_mod = imp.load_module(__name__, *mod_info)
                finally:
                    if mod_info[0] is not None:
                        mod_info[0].close()

        try:
            if StrictVersion(six_mod.__version__) >= _SIX_MIN_VERSION:
                break
        except (AttributeError, ValueError):
            # Attribute error if the six module isn't what it should be and
            # doesn't have a .__version__; ValueError if the version string
            # exists but is somehow bogus/unparseable
            continue
    else:
        raise ImportError(
            "Astropy requires the 'six' module of minimum version {0}; "
            "normally this is bundled with the astropy package so if you get "
            "this warning consult the packager of your Astropy "
            "distribution.".format(_SIX_MIN_VERSION))

    # Using importlib does not overwrite this shim, so do it ourselves.
    this_module = sys.modules[__name__]
    if not hasattr(this_module, '_importer'):
        # Copy all main six attributes.
        for name, value in six_mod.__dict__.items():
            if name.startswith('__'):
                continue
            this_module.__dict__[name] = value

        # Tell six's importer to accept this shim's name as its own.
        importer = six_mod._importer
        known_modules = list(importer.known_modules.items())
        for name, mod in known_modules:
            this_name = __name__ + name[len(mod_name):]
            importer.known_modules[this_name] = mod

        # Turn this shim into a package just like six does.
        this_module.__path__ = []  # required for PEP 302 and PEP 451
        this_module.__package__ = __name__  # see PEP 366
        if this_module.__dict__.get('__spec__') is not None:
            this_module.__spec__.submodule_search_locations = []  # PEP 451

_import_six()
