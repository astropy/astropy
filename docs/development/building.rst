.. _dev-build-astropy-subpkg:

************************************
Building Astropy and its Subpackages
************************************

The build process currently uses the `setuptools
<https://setuptools.readthedocs.io>`_ package to build and install the
astropy core (and any affiliated packages that use the template).


Astropy-helpers
===============

As of Astropy v0.4, Astropy also uses an external package called
`astropy-helpers <https://github.com/astropy/astropy-helpers>`_ to provide some
of its build and installation functionality.  A copy of astropy-helpers is
included with the Astropy source distribution, but also includes a mechanism to
automatically download bug fixes from PyPI.  The reason for providing these
helpers as a separate package is that it makes it easier for affiliated
packages to take advantage of these same utilities without requiring Astropy to
be installed *first*.  See `APE4
<https://github.com/astropy/astropy-APEs/blob/master/APE4.rst>`_ for the full
background on this.

Astropy-helpers is automatically bootstrapped to the Astropy build/installation
script (``setup.py``) via a script called ``ah_bootstrap.py`` that is imported
by ``setup.py``.  This script will do its best to ensure that the user has an
up-to-date copy of astropy-helpers before building the package.  The
auto-upgrade mechanism in particular allows pushing platform-specific fixes for
the build process without releasing a new version of Astropy (or any affiliated
package that uses astropy-helpers).

The behavior of the ``ah_bootstrap.py`` script can also be modified by options
in the project's ``setup.cfg`` file under a section called ``[ah_boostrap]``.
APE4 provides `more details
<https://github.com/astropy/astropy-APEs/blob/master/APE4.rst#astropy_helpers-bootstrap-script>`_.

The astropy-helpers distribution provides a Python package called
``astropy_helpers``.  Code that previously referenced the modules
``astropy.setup_helpers`` and ``astropy.version_helpers`` should now depend on
astropy-helpers and use ``astrop_helpers.setup_helpers`` and
``astropy_helpers.version_helpers`` respectively.  Likewise, astropy-helpers
includes tools for building Astropy's documentation.  The ``astropy.sphinx``
package is deprecated in favor of ``astropy_helpers.sphinx``.  As such,
astropy-helpers is a dependency of building Astropy's documentation.


Customizing setup/build for subpackages
=======================================

As is typical, there is a single ``setup.py`` file that is used for the whole
``astropy`` package.  To customize setup parameters for a given sub-package, a
``setup_package.py`` file can be defined inside a package, and if it is present,
the setup process will look for the following functions to customize the build
process:

* ``get_package_data``
    This function, if defined, should return a dictionary mapping the name of
    the subpackage(s) that need package data to a list of data file paths
    (possibly including wildcards) relative to the path of the package's source
    code.  e.g. if the source distribution has a needed data file
    ``astropy/wcs/tests/data/3d_cd.hdr``, this function should return
    ``{'astropy.wcs.tests':['data/3d_cd.hdr']}``. See the ``package_data``
    option of the  :func:`distutils.core.setup` function.

    It is recommended that all such data be in a directory named ``data`` inside
    the package within which it is supposed to be used.  This package data should
    be accessed via the `astropy.utils.data.get_pkg_data_filename` and
    `astropy.utils.data.get_pkg_data_fileobj` functions.

* ``get_extensions``
    This provides information for building C or Cython extensions. If defined,
    it should return a list of `distutils.core.Extension` objects controlling
    the Cython/C build process (see below for more detail).

* ``get_build_options``
    This function allows a package to add extra build options.  It
    should return a list of tuples, where each element has:

    - *name*: The name of the option as it would appear on the
      commandline or in the ``setup.cfg`` file.

    - *doc*: A short doc string for the option, displayed by
      ``setup.py build --help``.

    - *is_bool* (optional): When `True`, the option is a boolean
      option and doesn't have an associated value.

    Once an option has been added, its value can be looked up using
    ``astropy_helpers.setup_helpers.get_distutils_build_option``.

* ``get_external_libraries``
    This function declares that the package uses libraries that are
    included in the astropy distribution that may also be distributed
    elsewhere on the users system.  It should return a list of library
    names.  For each library, a new build option is created,
    ``'--use-system-X'`` which allows the user to request to use the
    system's copy of the library.  The package would typically call
    ``astropy_helpers.setup_helpers.use_system_library`` from its
    ``get_extensions`` function to determine if the package should use
    the system library or the included one.

The ``astropy_helpers.setup_helpers`` modules includes an
``update_package_files`` function which automatically searches the given source
path for ``setup_package.py`` modules and calls each of the above functions, if
they exist.  This makes it easy for affiliated packages to use this machinery
in their own ``setup.py``.


PLY Parsing/Lexing tables
=========================

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
* Run ``pip install -e .`` to build ``astropy``.
