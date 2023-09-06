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

import configparser
import doctest
import os
import sys
from datetime import datetime
from importlib import metadata

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet

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
    print(
        "The following packages could not be found and are required to "
        "build the documentation:"
    )
    for key, val in missing_requirements.items():
        print(f"    * {key} {val}")
    print('Please install the "docs" requirements.')
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
    rst_epilog,
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
        "astropy-dev": ("https://docs.astropy.org/en/latest/", None),
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

extensions += ["sphinx_changelog", "sphinx_design"]

# Grab minversion from setup.cfg
setup_cfg = configparser.ConfigParser()
setup_cfg.read(os.path.join(os.path.pardir, "setup.cfg"))
__minimum_python_version__ = setup_cfg["options"]["python_requires"].replace(">=", "")

min_versions = {}
for line in metadata.requires("astropy"):
    req = Requirement(line.split(";")[0])
    min_versions[req.name.lower()] = str(req.specifier)


# This is added to the end of RST files - a good place to put substitutions to
# be used globally.
with open("common_links.txt") as cl:
    rst_epilog += cl.read().format(
        minimum_python=__minimum_python_version__, **min_versions
    )

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
copyright = f"2011–{datetime.utcnow().year}, " + author

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
    exclude_patterns += ["development/*", "testhelpers.rst"]

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

# Setting this URL is requited by sphinx-astropy
github_issues_url = "https://github.com/astropy/astropy/issues/"
edit_on_github_branch = "main"

# Enable nitpicky mode - which ensures that all references in the docs
# resolve.

nitpicky = True
# See docs/nitpick-exceptions file for the actual listing.
nitpick_ignore = []
for line in open("nitpick-exceptions"):
    if line.strip() == "" or line.startswith("#"):
        continue
    dtype, target = line.split(None, 1)
    nitpick_ignore.append((dtype, target.strip()))

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
    "https://physics.nist.gov/",  # SSL: CERTIFICATE_VERIFY_FAILED
    "https://ieeexplore.ieee.org/",  # 418 Client Error: I'm a teapot
    "https://pyfits.readthedocs.io/en/v3.2.1/",  # defunct page in CHANGES.rst
    r"https://github\.com/astropy/astropy/(?:issues|pull)/\d+",
]
linkcheck_timeout = 180
linkcheck_anchors = False


def rstjinja(app, docname, source):
    """Render pages as a jinja template to hide/show dev docs."""
    # Make sure we're outputting HTML
    if app.builder.format != "html":
        return
    files_to_render = ["index", "install", "development/index"]
    if docname in files_to_render:
        print(f"Jinja rendering {docname}")
        rendered = app.builder.templates.render_string(
            source[0], app.config.html_context
        )
        source[0] = rendered


def resolve_astropy_and_dev_reference(app, env, node, contnode):
    """
    Reference targets for ``astropy:`` and ``astropy-dev:`` are special cases.

    Documentation links in astropy can be set up as intersphinx links so that
    affiliate packages do not have to override the docstrings when building
    the docs.

    If we are building the development docs it is a local ref targeting the
    label ``astropy-dev:<label>``, but for stable docs it should be an
    intersphinx resolution to the development docs.

    See https://github.com/astropy/astropy/issues/11366
    """
    # should the node be processed?
    reftarget = node.get("reftarget")  # str or None
    if str(reftarget).startswith("astropy:"):
        # This allows Astropy to use intersphinx links to itself and have
        # them resolve to local links. Downstream packages will see intersphinx.
        # TODO! deprecate this if sphinx-doc/sphinx/issues/9169 is implemented.
        process, replace = True, "astropy:"
    elif dev and str(reftarget).startswith("astropy-dev:"):
        process, replace = True, "astropy-dev:"
    else:
        process, replace = False, ""

    # make link local
    if process:
        reftype = node.get("reftype")
        refdoc = node.get("refdoc", app.env.docname)
        # convert astropy intersphinx targets to local links.
        # there are a few types of intersphinx link patterns, as described in
        # https://docs.readthedocs.io/en/stable/guides/intersphinx.html
        reftarget = reftarget.replace(replace, "")
        if reftype == "doc":  # also need to replace the doc link
            node.replace_attr("reftarget", reftarget)
        # Delegate to the ref node's original domain/target (typically :ref:)
        try:
            domain = app.env.domains[node["refdomain"]]
            return domain.resolve_xref(
                app.env, refdoc, app.builder, reftype, reftarget, node, contnode
            )
        except Exception:
            pass

        # Otherwise return None which should delegate to intersphinx


def setup(app):
    if sphinx_gallery is None:
        msg = (
            "The sphinx_gallery extension is not installed, so the "
            "gallery will not be built.  You will probably see "
            "additional warnings about undefined references due "
            "to this."
        )
        try:
            app.warn(msg)
        except AttributeError:
            # Sphinx 1.6+
            from sphinx.util import logging

            logger = logging.getLogger(__name__)
            logger.warning(msg)

    # Generate the page from Jinja template
    app.connect("source-read", rstjinja)
    # Set this to higher priority than intersphinx; this way when building
    # dev docs astropy-dev: targets will go to the local docs instead of the
    # intersphinx mapping
    app.connect("missing-reference", resolve_astropy_and_dev_reference, priority=400)
