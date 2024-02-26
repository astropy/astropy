.. _dev-build-astropy-subpkg:

************************************
Building Astropy and its Subpackages
************************************

The build process currently uses the `setuptools
<https://setuptools.readthedocs.io>`_ package to build and install the astropy
core (and any affiliated packages that use the template). As is typical, there
is a single ``setup.py`` file that is used for the whole ``astropy`` package. To
make it easier to set up C extensions for individual sub-packages, we use
`extension-helpers <https://extension-helpers.readthedocs.io/>`_, which allows
extensions to be defined inside each sub-package.

The way extension-helpers works is that it looks for ``setup_package.py`` files
anywhere in the package, and then looks for a function called ``get_extensions``
inside each of these files. This function should return a list of
``setuptools.Extension`` objects, and these are combined into an
overall list of extensions to build.

For certain string-parsing tasks, Astropy uses the
`PLY <http://www.dabeaz.com/ply/>`_ tool.  PLY generates tables that speed up
the parsing process, which are checked into source code so they don't have to
be regenerated.  These tables can be recognized by having either ``lextab`` or
``parsetab`` in their names.  To regenerate these files (e.g. if a new version
of PLY is bundled with Astropy or some of the parsing code changes), the tables
need to be deleted and the appropriate parts of astropy re-imported and run. For
exact details, see the comments in the headers of the ``parsetab`` and
``lextab`` files.

.. _dev-build-astropy-subpkg-win:

Building on Windows
*******************

The most convenient option is to use Python installation from Miniconda. If you like
Unix-like commands, Git Bash, which comes installed with Git, complements
Miniconda pretty well, as long as Miniconda is installed with the option for
it to be available system-wide (the option that is not recommended by the
installer).

Since ``astropy`` contains C extensions, you also need to install Microsoft
Visual Studio (the latest available should work) so Python can access the
system C compiler.

Once everything is set up as above, you can proceed to build ``astropy``
from source in the ``conda`` environment in an OS-agnostic way. For example:

* Create a new ``conda`` environment.
* Go to the ``astropy`` code checkout directory.
* If you have not already, fetch all of the tags from the main repository.
  If you do not have the latest tag, your developer version number will be
  wrong.
* Run ``python -m pip install --editable .`` to build ``astropy``.
