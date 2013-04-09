====================================
Building Astropy and its Subpackages
====================================

The build process currently uses the
`Distribute <http://packages.python.org/distribute/>`_ package to build and
install the astropy core (and any affiliated packages that use the template).
The user doesn't necessarily need to have `distribute` installed, as it will
automatically bootstrap itself using the ``distribute_setup.py`` file in the
source distribution if it isn't installed for the user.

Customizing setup/build for subpackages
---------------------------------------

As is typical, there is a single ``setup.py`` file that is used for the whole
`astropy` package.  To customize setup parameters for a given sub-package, a
``setup_package.py`` file can be defined inside a package, and if it is present,
the setup process will look for the following functions to customize the build
process:

* :func:`get_package_data`
    This function, if defined, should return a dictionary mapping the name of
    the subpackage(s) that need package data to a list of data file paths
    (possibly including wildcards) relative to the path of the package's source
    code.  e.g. if the source distribution has a needed data file
    ``astropy/wcs/tests/data/3d_cd.hdr``, this function should return
    ``{'astropy.wcs.tests':['data/3d_cd.hdr']}``. See the ``package_data``
    option of the  :func:`distutils.core.setup` function.

    It is recommended that all such data be in a directory named ``data`` inside
    the package within which it is supposed to be used.  This package data should
    be accessed via the `astropy.utils.data.get_data_filename` and
    `astropy.utils.data.get_data_fileobj` functions.

* :func:`get_extensions`
    This provides information for building C or Cython extensions. If defined,
    it should return a list of `distutils.core.Extension` objects controlling
    the Cython/C build process (see below for more detail).

* :func:`get_legacy_alias`
    This function allows for the creation of `shims` that allow a
    subpackage to be imported under another name.  For example,
    `astropy.io.fits` used to be available under the namespace
    `pyfits`.  For backward compatibility, it is helpful to have it
    still importable under the old name.  Under most circumstances,
    this function should call `astropy.setup_helpers.add_legacy_alias`
    to generate a legacy module and then return what it returns.

* :func:`get_build_options`
    This function allows a package to add extra build options.  It
    should return a list of tuples, where each element has:

    - *name*: The name of the option as it would appear on the
      commandline or in the `setup.cfg` file.

    - *doc*: A short doc string for the option, displayed by
      `setup.py build --help`.

    - *is_bool* (optional): When `True`, the option is a boolean
      option and doesn't have an associated value.

    Once an option has been added, its value can be looked up using
    `astropy.setup_helpers.get_distutils_build_option`.

* :func:`get_external_libraries`
    This function declares that the package uses libraries that are
    included in the astropy distribution that may also be distributed
    elsewhere on the users system.  It should return a list of library
    names.  For each library, a new build option is created,
    `--use-system-X` which allows the user to request to use the
    system's copy of the library.  The package would typically call
    `astropy.setup_helpers.use_system_library` from its
    `get_extensions` function to determine if the package should use
    the system library or the included one.

The `astropy.setup_helpers` modules includes a :func:`update_package_files`
function which automatically searches the given source path for
``setup_package.py`` modules and calls each of the above functions, if they
exist.  This makes it easy for affiliated packages to use this machinery in
their own ``setup.py``.
