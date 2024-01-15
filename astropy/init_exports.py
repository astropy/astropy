"""Steps to initialize on astropy on import. Copied from __init__.py to
support lazy_loader"""

from pathlib import Path

from .version import version as __version__

# The location of the online documentation for astropy
# This location will normally point to the current released version of astropy
online_docs_root = "https://docs.astropy.org/en/{}/".format(
    "latest" if "dev" in __version__ else f"v{__version__}"
)


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
from . import __path__

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


__all__ = [
    "__bibtex__",
    "test",
    "log",
    "find_api_page",
    "online_help",
    "online_docs_root",
    "physical_constants",
    "astronomical_constants",
]
