============================================
Building, Cython/C Extensions, and Releasing
============================================

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
    ``{'astropy.wcs.tests:'['data/3d_cd.hdr']}``. See the ``package_data``
    option of the  :func:`distutils.core.setup` function.

    It is recommended that all such data be in a directory named "data" inside
    the package within which it is supposed to be used, and package data should
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

.. _building-c-or-cython-extensions:

C or Cython Extensions
----------------------

Astropy supports using C extensions for wrapping C libraries and Cython for
speeding up computationally-intensive calculations. Both Cython and C extension
building can be customized using the :func:`get_extensions` function of the
``setup_package.py`` file. If defined, this function must return a list of
`distutils.core.Extension` objects. The creation process is left to the
subpackage designer, and can be customized however is relevant for the
extensions in the subpackage.

While C extensions must always be defined through the :func:`get_extensions`
mechanism, Cython files (ending in ``.pyx``) are automatically located and
loaded in separate extensions if they are not in :func:`get_extensions`. For
Cython extensions located in this way, headers for numpy C functions are
included in the build, but no other external headers are included. ``.pyx``
files present in the extensions returned by :func:`get_extensions` are not
included in the list of extensions automatically generated extensions. Note
that this allows disabling a Cython file by providing an extension that
includes the Cython file, but giving it the special `name` 'cython_skip'. Any
extension with this package name will not be built by ``setup.py``.

.. note::

    If an :class:`~distutils.core.Extension` object is provided for Cython
    source files using the :func:`get_extensions` mechanism, it is very
    important that the ``.pyx`` files be given as the `source`, rather than the
    ``.c`` files generated by Cython.

Installing C header files
^^^^^^^^^^^^^^^^^^^^^^^^^

If your C extension needs to be linked from other third-party C code,
you probably want to install its header files along side the Python module.

    1) Create an `include` directory inside of your package for
       all of the header files.

    2) Use the :func:`get_package_data` hook in `setup_package.py` to
       install those header files.  For example, the `astropy.wcs`
       package has this::

           def get_package_data():
               return {'astropy.wcs': ['include/*.h']}

Preventing importing at build time
----------------------------------

In rare cases, some packages may need to be imported at build time.
Unfortunately, anything that requires a C or Cython extension or
processing through 2to3 will fail to import until the build phase has
completed.  In those cases, the `_ASTROPY_SETUP_` variable can be used
to determine if the package is being imported as part of the build and
choose to not import problematic modules.  `_ASTROPY_SETUP_` is
inserted into the builtins, and is `True` when inside of astropy's
`setup.py` script, and `False` otherwise.

For example, suppose there is a subpackage ``foo`` that needs to
import a module called ``version.py`` at build time in order to set
some version information, and also has a C extension, ``process``,
that will not be available in the source tree.  In this case,
``astropy/foo/__init__.py`` would probably want to check the value of
`_ASTROPY_SETUP_` before importing the C extension::

    if not _ASTROPY_SETUP_:
        from . import process

    from . import version

Release
-------

The current release procedure for Astropy involves a combination of an
automated release script and some manual steps.  Future versions will automate
more of the process, if not all.

One of the main steps in performing a release is to create a tag in the git
repository representing the exact state of the repository that represents the
version being released.  For Astropy we will always use `signed tags`_: A
signed tag is annotated with the name and e-mail address of the signer, a date
and time, and a checksum of the code in the tag.  This information is then
signed with a GPG private key and stored in the repository.

Using a signed tag ensures the integrity of the contents of that tag for the
future.  On a distributed VCS like git, anyone can create a tag of Astropy
called "0.1" in their repository--and where it's easy to monkey around even
after the tag has been created.  But only one "0.1" will be signed by one of
the Astropy project coordinators and will be verifiable with their public key.

Creating a GPG Signing Key and a Signed Tag
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Git uses GPG to created signed tags, so in order to perform an Astropy release
you will need GPG installed and will have to generated a signing key pair.
Most \*NIX installations come with GPG installed by default (as it is used to
verify the integrity of system packages).  If you don't have the ``gpg``
command, consult the documentation for your system on how to install it.

For OSX, GPG can be installed from MacPorts using ``sudo port install gnupg``.

To create a new public/private key pair, simply run::

    $ gpg --gen-key

This will take you through a few interactive steps. For the encryption
and expiry settings, it should be safe to use the default settings (I use
a key size of 4096 just because what does a couple extra kilobytes
hurt?) Enter your full name, preferably including your middle name or
middle initial, and an e-mail address that you expect to be active for a
decent amount of time. Note that this name and e-mail address must match
the info you provide as your git configuration, so you should either
choose the same name/e-mail address when you create your key, or update
your git configuration to match the key info. Finally, choose a very good
pass phrase that won't be easily subject to brute force attacks.


If you expect to use the same key for some time, it's good to make a backup of
both your public and private key::

    $ gpg --export --armor > public.key
    $ gpg --export-secret-key --armor > private.key

Back up these files to a trusted location--preferably a write-once physical
medium that can be stored safely somewhere.  One may also back up their keys to
a trusted online encrypted storage, though some might not find that secure
enough--it's up to you and what you're comfortable with.

Add your public key to a keyserver
""""""""""""""""""""""""""""""""""
Now that you have a public key, you can publish this anywhere you like--in your
e-mail, in a public code repository, etc.  You can also upload it to a
dedicated public OpenPGP keyserver.  This will store the public key
indefinitely (until you manually revoke it), and will be automatically synced
with other keyservers around the world.  That makes it easy to retrieve your
public key using the gpg command-line tool.

To do this you will need your public key's keyname.  To find this enter::

    $ gpg --list-keys

This will output something like::

    /path/to/.gnupg/pubring.gpg
    ---------------------------------------------
    pub   4096D/1234ABCD 2012-01-01
    uid                  Your Name <your_email>
    sub   4096g/567890EF 2012-01-01

The 8 digit hex number on the line starting with "pub"--in this example the
"1234ABCD" unique keyname for your public key.  To push it to a keyserver
enter::

    $ gpg --send-keys 1234ABCD

But replace the 1234ABCD with the keyname for your public key.  Most systems
come configured with a sensible default keyserver, so you shouldn't have to
specify any more than that.

Create a tag
""""""""""""
Now test creating a signed tag in git.  It's safe to experiment with this--you
can always delete the tag before pushing it to a remote repository::

    $ git tag -s v0.1 -m "Astropy version 0.1"

This will ask for the password to unlock your private key in order to sign
the tag with it.  Confirm that the default signing key selected by git is the
correct one (it will be if you only have one key).

Once the tag has been created, you can verify it with::

    $ git tag -v v0.1

This should output something like::

    object e8e3e3edc82b02f2088f4e974dbd2fe820c0d934
    type commit
    tag v0.1
    tagger Your Name <your_email> 1339779534 -0400

    Astropy version 0.1
    gpg: Signature made Fri 15 Jun 2012 12:59:04 PM EDT using DSA key ID 0123ABCD
    gpg: Good signature from "Your Name <your_email>"

You can use this to verify signed tags from any repository as long as you have
the signer's public key in your keyring.  In this case you signed the tag
yourself, so you already have your public key.

Note that if you are planning to do a release following the steps below, you
will want to delete the tag you just created, because the release script does
that for you.  You can delete this tag by doing::

    $ git tag -d v0.1

Release Procedure
^^^^^^^^^^^^^^^^^

The automated portion of the Astropy release procedure uses `zest.releaser`_
to create the tag and update the version.  zest.releaser is extendable through
hook functions--Astropy already includes a couple hook functions to modify the
default behavior, but future releases may be further automated through the
implementation of additional hook functions.  In order to use the hooks,
Astropy itself must be *installed* alongside zest.releaser.  It is recommended
to create a `virtualenv`_ specifically for this purpose.

This may seem like a lot of steps, but most of them won't be necessary to
repeat for each release.  The advantage of using an automated or semi-automated
procedure is that ensures a consistent release process each time.

 1. Update the list of contributors in the ``creditsandlicense.rst`` file. The
    easiest way to check this is do::

        $ git shortlog -s

    And just add anyone from that list who isn't already credited.

 2. Install virtualenv if you don't already have it.  See the linked virtualenv
    documentation for details.  Also, make sure that you have `cython`_
    installed, as you will need it to generate the .c files needed for the
    release.

 3. Create and activate a virtualenv::

    $ virtualenv --system-site-packages --distribute astropy-release
    $ source astropy-release/bin/activate

 4. Obtain a *clean* version of the Astropy repository.  That is, one
    where you don't have any intermediate build files.  Either use a fresh
    ``git clone`` or do ``git clean -dfx``.

 5. Be sure you're the "master" branch, and install Astropy into the
    virtualenv::

    $ python setup.py install

    This is necessary for two reasons.  First, the entry points for the
    releaser scripts need to be availale, and these are in the Astropy
    package. Second, the build process will generate .c files from the
    Cython .pyx files, and the .c files are necessary for the source
    distribution.

 6. Install zest.releaser into the virtualenv; use ``--upgrade --force`` to
    ensure that the latest version is installed in the virtualenv (if you're
    running a csh variant make sure to run ``rehash`` afterwards too)::

    $ pip install zest.releaser --upgrade --force

 7. Ensure that all changes to the code have been committed, then start the
    release by running::

    $ fullrelease

 8. You will be asked to enter the version to be released.  Press enter to
    accept the default (which will normally be correct) or enter a specific
    version string.  A diff will then be shown of CHANGES.rst and setup.py
    showing that a release date has been added to the changelog, and that the
    version has been updated in setup.py.  Enter 'Y' when asked to commit
    these changes.

 9. You will then be shown the command that will be run to tag the release.
    Enter 'Y' to confirm and run the command.

 10. When asked "Check out the tag (for tweaks or pypi/distutils server
     upload)" enter 'N': zest.releaser does not offer enough control yet over
     how the register and upload are performed so we will do this manually
     until the release scripts have been improved.

 11. You will be asked to enter a new development version.  Normally the next
     logical version will be selected--press enter to accept the default, or
     enter a specific version string.  Do not add ".dev" to the version, as
     this will be appended automatically (ignore the message that says ".dev0
     will be appended"--it will actually be ".dev" without the 0).  For
     example, if the just-released version was "0.1" the default next version
     will be "0.2".  If we want the next version to be, say "1.0" then that
     must be entered manually.

 12. You will be shown a diff of CHANGES.rst showing that a new section has
     been added for the new development version, and showing that the version
     has been updated in setup.py.  Enter 'Y' to commit these changes.

 13. When asked to push the changes to a remote repository, enter 'Y'.  This
     should complete the portion of the process that's automated at this point.

 14. Check out the tag of the released version.  For example::

     $ git checkout v0.1

 15. Create the source distribution by doing::

     $ python setup.py sdist

     Copy the produced ``.tar.gz`` somewhere and verify that you can unpack it,
     build it, and get all the tests to pass.  It would be best to create a new
     virtualenv in which to do this.

 16. Register the release on PyPI with::

     $ python setup.py register

 17. Upload the source distribution to PyPI; this is preceded by re-running
     the sdist command, which is necessary for the upload command to know
     which distribution to upload::

     $ python setup.py sdist upload

 18. Update the website to reflect the fact there is now a stable release.

 19. Update Readthedocs so that it builds docs for the corresponding github tag,
     and set the default page to the new release.

 20. Create a bug fix branch.  If the version just was released was a "X.Y.0"
     version ("0.1" or "0.2" for example--the final ".0" is typically ommitted)
     it is good to create a bug fix branch as well.  Starting from the tagged
     changset, just checkout a new branch and push it to the remote server.
     For example, after releasing version 0.1, do::

     $ git checkout -b v0.1.x

     Then edit ``setup.py`` so that the version is ``'0.1.1.dev'``, and commit
     that change. Then, do::

     $ git push upstream v0.1.x

    .. note::
        You may need to replace ``upstream`` here with ``astropy`` or
        whatever remote name you use for the main astropy repository.

     The purpose of this branch is for creating bug fix releases like "0.1.1"
     and "0.1.2", while allowing development of new features to continue in
     the master branch.  Only changesets that fix bugs without making
     significant API changes should be merged to the bug fix branches.

 21. Create a bug fix label on GitHub; this should have the same name as the
     just created bug fix branch.  This label should be applied to all issues
     that should be backported to the bug fix branch.


.. _signed tags: http://git-scm.com/book/en/Git-Basics-Tagging#Signed-Tags
.. _zest.releaser: http://pypi.python.org/pypi/zest.releaser
.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _cython: http://www.cython.org/

Creating a MacOS X Installer on a DMG
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``bdist_dmg`` command can be used to create a ``.dmg`` disk image for
MacOS X with a ``.pkg`` installer. In order to do this, you will need to
ensure that you have the following dependencies installed:

* `Numpy <http://www.numpy.org>`_
* `Sphinx <http://sphinx.pocoo.org>`_
* `bdist_mpkg <http://pypi.python.org/pypi/bdist_mpkg/>`_

To create a ``.dmg`` file, run::

    python setup.py bdist_dmg

Note that for the actual release version, you should do this with the Python
distribution from `python.org <http://python.org>`_ (not e.g. MacPorts, EPD,
etc.). The best way to ensure maximum compatibility is to make sure that
Python and Numpy are installed into ``/Library/Frameworks/Python.framework``
using the latest stable ``.dmg`` installers available for those packages. In
addition, the ``.dmg`` should be build on a MacOS 10.6 system, to ensure
compatibility with 10.6, 10.7, and 10.8.

Before distributing, you should test out an installation of Python, Numpy, and
Astropy from scratch using the ``.dmg`` installers, preferably on a clean
virtual machine.

Future directions
-----------------

We plan to switch to a newer packaging scheme when it's more stable, the
upcoming standard library `packaging` module, derived from the
`distutils2 <http://packages.python.org/Distutils2/library/distutils2.html>`_
project.  Until it's working right, however, we will be using `distribute` and
`distutils`.

.. include:: links.inc
