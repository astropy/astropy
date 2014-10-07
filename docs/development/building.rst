====================================
Building Astropy and its Subpackages
====================================

The build process currently uses the `setuptools
<https://bitbucket.org/pypa/setuptools>`_ package to build and install the
astropy core (and any affiliated packages that use the template).  The user
doesn't necessarily need to have `setuptools`_ installed, as it will
automatically bootstrap itself using the ``ez_setup.py`` file in the source
distribution if it isn't installed for the user.


Astropy-helpers
---------------

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
---------------------------------------

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

* ``requires_2to3``
    This function declares whether the package requires processing
    through the `2to3`_ tool to run on Python 3.  If not included, it
    defaults to `True`.  The use of `2to3`_ is being phased out in
    astropy, in favor of using `six`_ instead.  See :ref:`dev-portable`
    for more information.

The ``astropy_helpers.setup_helpers`` modules includes an
``update_package_files`` function which automatically searches the given source
path for ``setup_package.py`` modules and calls each of the above functions, if
they exist.  This makes it easy for affiliated packages to use this machinery
in their own ``setup.py``.

.. _six: http://pythonhosted.org/six/
.. _2to3: https://docs.python.org/2/library/2to3.html

