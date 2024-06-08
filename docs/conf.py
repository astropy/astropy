# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Astropy documentation build configuration file.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this file.
#
# All configuration values have a default. Some values are defined in
# the global Astropy configuration which is loaded here before anything else.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath('..'))
# IMPORTANT: the above commented section was generated by sphinx-quickstart, but
# is *NOT* appropriate for astropy or Astropy affiliated packages. It is left
# commented out with this explanation to make it clear why this should not be
# done. If the sys.path entry above is added, when the astropy.sphinx.conf
# import occurs, it will import the *source* version of astropy instead of the
# version installed (if invoked as "make html" or directly with sphinx), or the
# version in the build directory.
# Thus, any C-extensions that are needed to build the documentation will *not*
# be accessible, and the documentation will not build correctly.
# See sphinx_astropy.conf for which values are set there.

import doctest
import os
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from sphinx.util import logging

# from docs import global_substitutions


if sys.version_info < (3, 11):
    import tomli as tomllib
else:
    import tomllib

logger = logging.getLogger(__name__)

# -- Check for missing dependencies -------------------------------------------
missing_requirements = {}
for line in metadata.requires("astropy"):
    if 'extra == "docs"' in line:
        req = Requirement(line.split(";")[0])
        req_package = req.name.lower()
        req_specifier = str(req.specifier)

        try:
            version = metadata.version(req_package)
        except metadata.PackageNotFoundError:
            missing_requirements[req_package] = req_specifier

        if version not in SpecifierSet(req_specifier, prereleases=True):
            missing_requirements[req_package] = req_specifier

if missing_requirements:
    msg = (
        "The following packages could not be found and are required to "
        "build the documentation:\n"
        "%s"
        '\nPlease install the "docs" requirements.',
        "\n".join([f"    * {key} {val}" for key, val in missing_requirements.items()]),
    )
    logger.error(msg)
    sys.exit(1)

from sphinx_astropy.conf.v2 import *  # noqa: E402, F403
from sphinx_astropy.conf.v2 import (  # noqa: E402
    exclude_patterns,
    extensions,
    html_theme_options,
    intersphinx_mapping,
    numpydoc_xref_aliases,
    numpydoc_xref_astropy_aliases,
    numpydoc_xref_ignore,
)

# -- Plot configuration -------------------------------------------------------
plot_rcparams = {
    "axes.labelsize": "large",
    "figure.figsize": (6, 6),
    "figure.subplot.hspace": 0.5,
    "savefig.bbox": "tight",
    "savefig.facecolor": "none",
}
plot_apply_rcparams = True
plot_html_show_source_link = False
plot_formats = ["png", "svg", "pdf"]
# Don't use the default - which includes a numpy and matplotlib import
plot_pre_code = ""

# -- General configuration ----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "3.0"

# The intersphinx_mapping in sphinx_astropy.sphinx refers to astropy for
# the benefit of other packages who want to refer to objects in the
# astropy core.  However, we don't want to cyclically reference astropy in its
# own build so we remove it here.
del intersphinx_mapping["astropy"]

# add any custom intersphinx for astropy
intersphinx_mapping.update(
    {
        "pyerfa": ("https://pyerfa.readthedocs.io/en/stable/", None),
        "pytest": ("https://docs.pytest.org/en/stable/", None),
        "ipython": ("https://ipython.readthedocs.io/en/stable/", None),
        "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
        "sphinx_automodapi": (
            "https://sphinx-automodapi.readthedocs.io/en/stable/",
            None,
        ),
        "asdf-astropy": ("https://asdf-astropy.readthedocs.io/en/latest/", None),
        "fsspec": ("https://filesystem-spec.readthedocs.io/en/latest/", None),
        "cycler": ("https://matplotlib.org/cycler/", None),
    }
)

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# .inc.rst mean *include* files, don't have sphinx process them
exclude_patterns += ["_templates", "changes", "_pkgtemplate.rst", "**/*.inc.rst"]

# Add any paths that contain templates here, relative to this directory.
if "templates_path" not in locals():  # in case parent conf.py defines it
    templates_path = []
templates_path.append("_templates")

extensions += ["sphinx_changelog", "sphinx_design", "sphinxcontrib.globalsubs"]

# Grab minversion from pyproject.toml
with (Path(__file__).parents[1] / "pyproject.toml").open("rb") as f:
    pyproject = tomllib.load(f)

# Manually register doctest options since matplotlib 3.5 messed up allowing them
# from pytest-doctestplus
IGNORE_OUTPUT = doctest.register_optionflag("IGNORE_OUTPUT")
REMOTE_DATA = doctest.register_optionflag("REMOTE_DATA")
FLOAT_CMP = doctest.register_optionflag("FLOAT_CMP")

# Whether to create cross-references for the parameter types in the
# Parameters, Other Parameters, Returns and Yields sections of the docstring.
numpydoc_xref_param_type = True

# Words not to cross-reference. Most likely, these are common words used in
# parameter type descriptions that may be confused for classes of the same
# name. The base set comes from sphinx-astropy. We add more here.
numpydoc_xref_ignore.update(
    {
        "mixin",
        "Any",  # aka something that would be annotated with `typing.Any`
        # needed in subclassing numpy  # TODO! revisit
        "Arguments",
        "Path",
        # TODO! not need to ignore.
        "flag",
        "bits",
    }
)

# Mappings to fully qualified paths (or correct ReST references) for the
# aliases/shortcuts used when specifying the types of parameters.
# Numpy provides some defaults
# https://github.com/numpy/numpydoc/blob/b352cd7635f2ea7748722f410a31f937d92545cc/numpydoc/xref.py#L62-L94
# and a base set comes from sphinx-astropy.
# so here we mostly need to define Astropy-specific x-refs
numpydoc_xref_aliases.update(
    {
        # python & adjacent
        "Any": "`~typing.Any`",
        "file-like": ":term:`python:file-like object`",
        "file": ":term:`python:file object`",
        "path-like": ":term:`python:path-like object`",
        "module": ":term:`python:module`",
        "buffer-like": ":term:buffer-like",
        "hashable": ":term:`python:hashable`",
        # for matplotlib
        "color": ":term:`color`",
        # for numpy
        "ints": ":class:`python:int`",
        # for astropy
        "number": ":term:`number`",
        "Representation": ":class:`~astropy.coordinates.BaseRepresentation`",
        "writable": ":term:`writable file-like object`",
        "readable": ":term:`readable file-like object`",
        "BaseHDU": ":doc:`HDU </io/fits/api/hdus>`",
    }
)
# Add from sphinx-astropy 1) glossary aliases 2) physical types.
numpydoc_xref_aliases.update(numpydoc_xref_astropy_aliases)

# Turn off table of contents entries for functions and classes
toc_object_entries = False

# -- Project information ------------------------------------------------------

project = "Astropy"
author = "The Astropy Developers"
copyright = f"2011–{datetime.now(tz=timezone.utc).year}, " + author

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

# The full version, including alpha/beta/rc tags.
release = metadata.version(project)
# The short X.Y version.
version = ".".join(release.split(".")[:2])

# Only include dev docs in dev version.
dev = "dev" in release
if not dev:
    exclude_patterns += ["development/*"]

# -- Options for the module index ---------------------------------------------

modindex_common_prefix = ["astropy."]


# -- Options for HTML output ---------------------------------------------------

html_theme_options.update(
    {
        "github_url": "https://github.com/astropy/astropy",
        "external_links": [
            {"name": "Tutorials", "url": "https://learn.astropy.org/"},
        ],
        "use_edit_page_button": True,
        "logo": {
            "image_light": "_static/astropy_banner_96.png",
            "image_dark": "_static/astropy_banner_96_dark.png",
        },
        # https://github.com/pydata/pydata-sphinx-theme/issues/1492
        "navigation_with_keys": False,
    }
)

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = f"{project} v{release}"

html_favicon = "_static/astropy_logo.ico"

html_css_files = ["astropy.css"]
html_copy_source = False

# Output file base name for HTML help builder.
htmlhelp_basename = project + "doc"

# A dictionary of values to pass into the template engine's context for all pages.
html_context = {
    "default_mode": "light",
    "to_be_indexed": ["stable", "latest"],
    "is_development": dev,
    "github_user": "astropy",
    "github_repo": "astropy",
    "github_version": "main",
    "doc_path": "docs",
}

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
html_extra_path = ["robots.txt"]

# -- Options for LaTeX output --------------------------------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
    ("index", project + ".tex", project + " Documentation", author, "manual")
]

latex_logo = "_static/astropy_logo.pdf"


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [("index", project.lower(), project + " Documentation", [author], 1)]

# Setting this URL is required by sphinx-astropy
github_issues_url = "https://github.com/astropy/astropy/issues/"
edit_on_github_branch = "main"

# Enable nitpicky mode - which ensures that all references in the docs
# resolve.

nitpicky = True
show_warning_types = True
# See docs/nitpick-exceptions file for the actual listing.
nitpick_ignore = []
for line in open("nitpick-exceptions"):
    if line.strip() == "" or line.startswith("#"):
        continue
    dtype, target = line.split(None, 1)
    nitpick_ignore.append((dtype, target.strip()))

suppress_warnings = [
    "config.cache",  # our rebuild is okay
]

# -- Options for the Sphinx gallery -------------------------------------------

try:
    import warnings

    import sphinx_gallery

    extensions += ["sphinx_gallery.gen_gallery"]

    sphinx_gallery_conf = {
        "backreferences_dir": "generated/modules",  # path to store the module using example template
        "filename_pattern": "^((?!skip_).)*$",  # execute all examples except those that start with "skip_"
        "examples_dirs": f"..{os.sep}examples",  # path to the examples scripts
        "gallery_dirs": "generated/examples",  # path to save gallery generated examples
        "reference_url": {
            "astropy": None,
            "matplotlib": "https://matplotlib.org/stable/",
            "numpy": "https://numpy.org/doc/stable/",
        },
        "abort_on_example_error": True,
    }

    # Filter out backend-related warnings as described in
    # https://github.com/sphinx-gallery/sphinx-gallery/pull/564
    warnings.filterwarnings(
        "ignore",
        category=UserWarning,
        message=(
            "Matplotlib is currently using agg, which is a"
            " non-GUI backend, so cannot show the figure."
        ),
    )

except ImportError:
    sphinx_gallery = None


# -- Options for linkcheck output -------------------------------------------
linkcheck_retry = 5
linkcheck_ignore = [
    "https://journals.aas.org/manuscript-preparation/",
    "https://maia.usno.navy.mil/",
    "https://www.usno.navy.mil/USNO/time/gps/usno-gps-time-transfer",
    "https://aa.usno.navy.mil/publications/docs/Circular_179.php",
    "http://data.astropy.org",
    "https://doi.org/",  # CI blocked by service provider
    "https://ui.adsabs.harvard.edu",  # CI blocked by service provider
    "https://www.tandfonline.com/",  # 403 Client Error: Forbidden
    "https://stackoverflow.com/",  # 403 Client Error: Forbidden
    "https://ieeexplore.ieee.org/",  # 418 Client Error: I'm a teapot
    "https://pyfits.readthedocs.io/en/v3.2.1/",  # defunct page in CHANGES.rst
    "https://pkgs.dev.azure.com/astropy-project",  # defunct page in CHANGES.rst
    r"https://github\.com/astropy/astropy/(?:issues|pull)/\d+",
]
linkcheck_timeout = 180
linkcheck_anchors = False
linkcheck_report_timeouts_as_broken = True
linkcheck_allow_unauthorized = False


def rstjinja(app, docname, source):
    """Render pages as a jinja template to hide/show dev docs."""
    # Make sure we're outputting HTML
    if app.builder.format != "html":
        return
    files_to_render = ["index_dev", "install"]
    if docname in files_to_render:
        logger.info("Jinja rendering %s", docname)
        rendered = app.builder.templates.render_string(
            source[0], app.config.html_context
        )
        source[0] = rendered


__minimum_python_version__ = pyproject["project"]["requires-python"].replace(">=", "")

min_versions = {}
for line in metadata.requires("astropy"):
    req = Requirement(line.split(";")[0])
    min_versions[req.name.lower()] = str(req.specifier)

# The following global_substitutions can be used throughout the
# documentation via sphinxcontrib-globalsubs. The key to the dictionary
# is the name of the case-sensitive substitution. For example, if the
# key is `"SkyCoord"`, then it can be used as `|SkyCoord|` throughout
# the documentation.

global_substitutions: dict[str, str] = {
    # NumPy
    "ndarray": ":class:`numpy.ndarray`",
    # Coordinates
    "EarthLocation": ":class:`~astropy.coordinates.EarthLocation`",
    "Angle": "`~astropy.coordinates.Angle`",
    "Latitude": "`~astropy.coordinates.Latitude`",
    "Longitude": ":class:`~astropy.coordinates.Longitude`",
    "BaseFrame": "`~astropy.coordinates.BaseCoordinateFrame`",
    "SkyCoord": ":class:`~astropy.coordinates.SkyCoord`",
    "SpectralCoord": "`~astropy.coordinates.SpectralCoord`",
    # Cosmology
    "Cosmology": ":class:`~astropy.cosmology.Cosmology`",
    "Cosmology.read": ":meth:`~astropy.cosmology.Cosmology.read`",
    "Cosmology.write": ":meth:`~astropy.cosmology.Cosmology.write`",
    "Cosmology.from_format": ":meth:`~astropy.cosmology.Cosmology.from_format`",
    "Cosmology.to_format": ":meth:`~astropy.cosmology.Cosmology.to_format`",
    "FLRW": ":class:`~astropy.cosmology.FLRW`",
    "LambdaCDM": ":class:`~astropy.cosmology.LambdaCDM`",
    "FlatLambdaCDM": ":class:`~astropy.cosmology.FlatLambdaCDM`",
    "WMAP1": ":ref:`astropy_cosmology_realizations_WMAP1`",
    "WMAP3": ":ref:`astropy_cosmology_realizations_WMAP3`",
    "WMAP5": ":ref:`astropy_cosmology_realizations_WMAP5`",
    "WMAP7": ":ref:`astropy_cosmology_realizations_WMAP7`",
    "WMAP9": ":ref:`astropy_cosmology_realizations_WMAP9`",
    "Planck13": ":ref:`astropy_cosmology_realizations_Planck13`",
    "Planck15": ":ref:`astropy_cosmology_realizations_Planck15`",
    "Planck18": ":ref:`astropy_cosmology_realizations_Planck18`",
    "FlatCosmologyMixin": ":class:`~astropy.cosmology.FlatCosmologyMixin`",
    "FlatFLRWMixin": ":class:`~astropy.cosmology.FlatFLRWMixin`",
    "default_cosmology": ":class:`~astropy.cosmology.default_cosmology`",
    # SAMP
    "SAMPClient": ":class:`~astropy.samp.SAMPClient`",
    "SAMPIntegratedClient": ":class:`~astropy.samp.SAMPIntegratedClient`",
    "SAMPHubServer": ":class:`~astropy.samp.SAMPHubServer`",
    "SAMPHubProxy": ":class:`~astropy.samp.SAMPHubProxy`",
    # Table
    "Column": ":class:`~astropy.table.Column`",
    "MaskedColumn": ":class:`~astropy.table.MaskedColumn`",
    "TableColumns": ":class:`~astropy.table.TableColumns`",
    "Row": ":class:`~astropy.table.Row`",
    "Table": ":class:`~astropy.table.Table`",
    "QTable": ":class:`~astropy.table.QTable`",
    # Time
    "Time": ":class:`~astropy.time.Time`",
    "TimeDelta": ":class:`~astropy.time.TimeDelta`",
    # Timeseries
    "TimeSeries": ":class:`~astropy.timeseries.TimeSeries`",
    "BinnedTimeSeries": ":class:`~astropy.timeseries.BinnedTimeSeries`",
    # Distribution
    "Distribution": ":class:`~astropy.uncertainty.Distribution`",
    # Units
    "PhysicalType": ":class:`~astropy.units.PhysicalType`",
    "Quantity": ":class:`~astropy.units.Quantity`",
    "Unit": ":class:`~astropy.units.UnitBase`",
    "StructuredUnit": ":class:`~astropy.units.StructuredUnit`",
    # Utils
    "Masked": ":class:`~astropy.utils.masked.Masked`",
    # Minimum versions
    "minimum_python_version": f"{__minimum_python_version__}",
    "minimum_numpy_version": f"{min_versions['numpy']}",
    "minimum_pyerfa_version": f"{min_versions['pyerfa']}",
    "minimum_matplotlib_version": f"{min_versions['matplotlib']}",
    "minimum_scipy_version": f"{min_versions['scipy']}",
    "minimum_asdf_astropy_version": f"{min_versions['asdf-astropy']}",
    "minimum_packaging_version": f"{min_versions['packaging']}",
    "minimum_pyyaml_version": f"{min_versions['pyyaml']}",
    "minimum_ipython_version": f"{min_versions['ipython']}",
    "minimum_pyarrow_version": f"{min_versions['pyarrow']}",
    "minimum_fsspec_version": f"{min_versions['fsspec']}",
    "minimum_s3fs_version": f"{min_versions['s3fs']}",
}
# Because sphinxcontrib-globalsubs does not work for regular reStructuredText
# links, we first define the links and then process them into the form
# of a reStructuredText external link.

links_to_become_substitutions: dict[str, str] = {
    # Python
    "Python": "https://www.python.org",
    "PEP8": "https://www.python.org/dev/peps/pep-0008",
    # Astropy
    "Astropy mailing list": "https://mail.python.org/mailman/listinfo/astropy",
    "astropy-dev mailing list": "http://groups.google.com/group/astropy-dev",
    # NumPy
    "NumPy": "https://numpy.org",
    "numpydoc": "https://pypi.org/project/numpydoc",
    # erfa
    "ERFA": "https://github.com/liberfa/erfa",
    "PyERFA": "http://pyerfa.readthedocs.org",
    # matplotlib
    "Matplotlib": "https://matplotlib.org",
    # sofa
    "SOFA": "http://www.iausofa.org/index.html",
    # scipy
    "SciPy": "https://www.scipy.org",
    # packaging
    "packaging": "https://packaging.pypa.io",
    # IPython
    "IPython": "https://ipython.org",
    # pip
    "pip": "https://pip.pypa.io",
    # pipenv
    "pipenv": "https://pipenv.pypa.io/en/latest",
    # virtualenv
    "virtualenv": "https://pypi.org/project/virtualenv",
    "virtualenvwrapper": "https://pypi.org/project/virtualenvwrapper",
    # conda
    "conda": "https://conda.io/docs",
    "miniconda": "https://docs.conda.io/en/latest/miniconda.html",
    # pytest
    "pytest": "https://pytest.org/en/latest/index.html",
    "pytest-astropy": "https://github.com/astropy/pytest-astropy",
    "pytest-doctestplus": "https://github.com/astropy/pytest-doctestplus",
    "pytest-remotedata": "https://github.com/astropy/pytest-remotedata",
    # fsspec
    "fsspec": "https://filesystem-spec.readthedocs.io",
    # s3fs
    "s3fs": "https://s3fs.readthedocs.io",
    # TOPCAT
    "STIL": "http://www.starlink.ac.uk/stil",
    "STILTS": "http://www.starlink.ac.uk/stilts",
    "TOPCAT": "http://www.starlink.ac.uk/topcat",
    # OpenAstronomy
    "OpenAstronomy Packaging Guide": "https://packaging-guide.openastronomy.org/en/latest",
}

processed_links = {
    key: f"`{key} <{value}>`_" for key, value in links_to_become_substitutions.items()
}

global_substitutions |= processed_links


def setup(app):
    if sphinx_gallery is None:
        logger.warning(
            "The sphinx_gallery extension is not installed, so the "
            "gallery will not be built.  You will probably see "
            "additional warnings about undefined references due "
            "to this."
        )

    # Generate the page from Jinja template
    app.connect("source-read", rstjinja)
