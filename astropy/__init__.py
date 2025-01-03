# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Astropy is a package intended to contain core functionality and some
common tools needed for performing astronomy and astrophysics research with
Python. It also provides an index for other astronomy packages and tools for
managing them.
"""

import sys
from pathlib import Path

from astropy.utils.system_info import system_info

__all__ = [  # noqa: RUF100, RUF022
    "__version__",
    "__bibtex__",
    # Subpackages (mostly lazy-loaded)
    "config",
    "constants",
    "convolution",
    "coordinates",
    "cosmology",
    "io",
    "modeling",
    "nddata",
    "samp",
    "stats",
    "table",
    "tests",
    "time",
    "timeseries",
    "uncertainty",
    "units",
    "utils",
    "visualization",
    "wcs",
    # Functions
    "test",
    "log",
    "find_api_page",
    "online_help",
    "online_docs_root",
    "conf",
    "physical_constants",
    "astronomical_constants",
    "system_info",
]


def __getattr__(attr):
    if attr in __all__:
        from importlib import import_module

        return import_module("astropy." + attr)

    raise AttributeError(f"module 'astropy' has no attribute {attr!r}")


def __dir__():
    return sorted(set(globals()).union(__all__))


from .version import version as __version__

# The location of the online documentation for astropy
# This location will normally point to the current released version of astropy
online_docs_root = "https://docs.astropy.org/en/{}/".format(
    "latest" if "dev" in __version__ else f"v{__version__}"
)


from . import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy`.
    """

    unicode_output = _config.ConfigItem(
        False,
        "When True, use Unicode characters when outputting values, and "
        "displaying widgets at the console.",
    )
    use_color = _config.ConfigItem(
        sys.platform != "win32",
        "When True, use ANSI color escape sequences when writing to the console.",
        aliases=["astropy.utils.console.USE_COLOR", "astropy.logger.USE_COLOR"],
    )
    max_lines = _config.ConfigItem(
        None,
        description=(
            "Maximum number of lines in the display of pretty-printed "
            "objects. If not provided, try to determine automatically from the "
            "terminal size.  Negative numbers mean no limit."
        ),
        cfgtype="integer(default=None)",
        aliases=["astropy.table.pprint.max_lines"],
    )
    max_width = _config.ConfigItem(
        None,
        description=(
            "Maximum number of characters per line in the display of "
            "pretty-printed objects.  If not provided, try to determine "
            "automatically from the terminal size. Negative numbers mean no "
            "limit."
        ),
        cfgtype="integer(default=None)",
        aliases=["astropy.table.pprint.max_width"],
    )


conf = Conf()


# Define a base ScienceState for configuring constants and units
from .utils.state import ScienceState


class base_constants_version(ScienceState):
    """
    Base class for the real version-setters below.
    """

    _value = "test"

    _versions = dict(test="test")

    @classmethod
    def validate(cls, value):
        if value not in cls._versions:
            raise ValueError(f"Must be one of {list(cls._versions.keys())}")
        return cls._versions[value]

    @classmethod
    def set(cls, value):
        """
        Set the current constants value.
        """
        import sys

        if "astropy.units" in sys.modules:
            raise RuntimeError("astropy.units is already imported")
        if "astropy.constants" in sys.modules:
            raise RuntimeError("astropy.constants is already imported")

        return super().set(value)


class physical_constants(base_constants_version):
    """
    The version of physical constants to use.
    """

    # Maintainers: update when new constants are added
    _value = "codata2018"

    _versions = dict(
        codata2018="codata2018",
        codata2014="codata2014",
        codata2010="codata2010",
        astropyconst40="codata2018",
        astropyconst20="codata2014",
        astropyconst13="codata2010",
    )


class astronomical_constants(base_constants_version):
    """
    The version of astronomical constants to use.
    """

    # Maintainers: update when new constants are added
    _value = "iau2015"

    _versions = dict(
        iau2015="iau2015",
        iau2012="iau2012",
        astropyconst40="iau2015",
        astropyconst20="iau2015",
        astropyconst13="iau2012",
    )


# Create the test() function
from .tests.runner import TestRunner

test = TestRunner.make_test_runner_in(__path__[0])


# if we are *not* in setup mode, import the logger and possibly populate the
# configuration file with the defaults
def _initialize_astropy():
    try:
        from .utils import _compiler
    except ImportError:
        # If this __init__.py file is in ./astropy/ then import is within a source
        # dir .astropy-root is a file distributed with the source, but that should
        # not installed
        if (Path(__file__).parent.parent / ".astropy-root").exists():
            raise ImportError(
                "You appear to be trying to import astropy from "
                "within a source checkout or from an editable "
                "installation without building the extension "
                "modules first. Either run:\n\n"
                "  pip install -e .\n\nor\n\n"
                "  python setup.py build_ext --inplace\n\n"
                "to make sure the extension modules are built "
            ) from None

        # Outright broken installation, just raise standard error
        raise


# Set the bibtex entry to the article referenced in CITATION.
def _get_bibtex():
    refs = (Path(__file__).parent / "CITATION").read_text().split("@ARTICLE")[1:]
    return f"@ARTICLE{refs[0]}" if refs else ""


__citation__ = __bibtex__ = _get_bibtex()

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
    import webbrowser
    from urllib.parse import urlencode

    url = online_docs_root + f"search.html?{urlencode({'q': query})}"
    webbrowser.open(url)


from types import ModuleType as __module_type__

# Clean up top-level namespace--delete everything that isn't in __all__
# or is a magic attribute, and that isn't a submodule of this package
for varname in dir():
    if not (
        (varname.startswith("__") and varname.endswith("__"))
        or varname in __all__
        or (
            varname[0] != "_"
            and isinstance(locals()[varname], __module_type__)
            and locals()[varname].__name__.startswith(__name__ + ".")
        )
    ):
        # The last clause in the above disjunction deserves explanation:
        # When using relative imports like ``from .. import config``, the
        # ``config`` variable is automatically created in the namespace of
        # whatever module ``..`` resolves to (in this case astropy).  This
        # happens a few times just in the module setup above.  This allows
        # the cleanup to keep any public submodules of the astropy package
        del locals()[varname]

del varname, __module_type__
