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


Maintaining Bug Fix Releases
----------------------------

Astropy releases, as recommended for most Python projects, follows a
<major>.<minor>.<micro> version scheme, where the "micro" version is also
known as a "bug fix" release.  Bug fix releases should not change any user-
visible interfaces.  They should only fix bugs on the previous major/minor
release and may also refactor internal APIs or include omissions from previous
releases--that is, features that were documented to exist but were accidentally
left out of the previous release.

Bug fix releases are typically managed by maintaining one or more bug fix
branches separate from the master branch (the release procedure below discusses
creating these branches).  Typically, whenever an issue is fixed on the Astropy
master branch a decision must be made whether this is a fix that should be
included in the Astropy bug fix release.  Usually the answer to this question
is "yes", though there are some issues that may not apply to the bug fix
branch.  For example, it is not necessary to backport a fix to a new feature
that did not exist when the bug fix branch was first created.  New features
are never merged into the bug fix branch--only bug fixes; hence the name.

In rare cases a bug fix may be made directly into the bug fix branch without
going into the master branch first.  This may occur if a fix is made to a
feature that has been removed or rewritten in the development version and no
longer has the issue being fixed.  However, depending on how critical the bug
is it may be worth including in a bug fix release, as some users can be slow to
upgrade to new major/micro versions due to API changes.

Issues are assigned to an Astropy release by way of the Milestone feature in
the GitHub issue tracker.  At any given time there are at least two versions
under development: The next major/minor version, and the next bug fix release.
For example, at the time of writing there are two release milestones open:
v0.2.2 and v0.3.0.  In this case, v0.2.2 is the next bug fix release and all
issues that should include fixes in that release should be assigned that
milestone.  Any issues that implement new features would go into the v0.3.0
milestone--this is any work that goes in the master branch that should not
be backported.  For a more detailed set of guidelines on using milestones, see
:ref:`milestones-and-labels`.

Backporting fixes from master
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most fixes are backported using the ``git cherry-pick`` command, which applies
the diff from a single commit like a patch.  For the sake of example, say the
current bug fix branch is 'v0.2.x', and that a bug was fixed in master in a
commit ``abcd1234``.  In order to backport the fix, simply checkout the v0.2.x
branch (it's also good to make sure it's in sync with the main Astropy
repository) and cherry-pick the appropriate commit::

    $ git checkout v0.2.x
    $ git pull upstream v0.2.x
    $ git cherry-pick abcd1234

Sometimes a cherry-pick does not apply cleanly, since the bug fix branch
represents a different line of development.  This can be resolved like any
other merge conflict:  Edit the conflicted files by hand, and then run
``git commit`` and accept the default commit message.

What if the issue required more than one commit to fix?  There are a few
possibilities for this.  The easiest is if the fix came in the form of a
pull request that was merged into the master branch.  Whenever GitHub merges
a pull request it generates a merge commit in the master branch.  This merge
commit represents the *full* difference of all the commits in the pull request
combined.  What this means is that it is only necessary to cherry-pick the
merge commit (this requires adding the ``-m 1`` option to the cherry-pick
command).  For example, if ``5678abcd`` is a merge commit::

    $ git checkout v0.2.x
    $ git pull upstream v0.2.x
    $ git cherry-pick -m 1 5678abcd

In fact, because Astropy emphasizes a pull request-based workflow, this is the
*most* common scenario for backporting bug fixes, and the one requiring the
least thought.  However, if you're not dealing with backporting a fix that was
not brought in as a pull request, read on.

.. seealso::

    :ref:`merge-commits-and-cherry-picks` for further explanation of the
    cherry-pick command and how it works with merge commits.

If not cherry-picking a merge commit there are still other options for dealing
with multiple commits.  The simplest, though potentially tedious, is to simply
run the cherry-pick command once for each commit in the correct order.
However, as of Git 1.7.2 it is possible to merge a range of commits like so::

    $ git cherry-pick 1234abcd..56789def

This works fine so long as the commits you want to pick are actually congruous
with each other.  In most cases this will be the case, though some bug fixes
will involve followup commits that need to back backported as well.  Most bug
fixes will have an issues associated with it in the issue tracker, so make sure
to reference all commits related to that issue in the commit message.  That way
it's harder for commits that need to be backported from getting lost.

Making fixes directly to the bug fix branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned earlier in this section, in some cases a fix only applies to a bug
fix release, and is not applicable in the mainline development.  In this case
there are two choices:

1. An Astropy developer with commit access to the main Astropy repository may
   check out the bug fix branch and commit and push your fix directly.

2. **Preferable**: You may also make a pull request through GitHub against the
   bug fix branch rather than against master.  Normally when making a pull
   request from a branch on your fork to the main Astropy repository GitHub
   compares your branch to Astropy's master.  If you look on the left-hand
   side of the pull request page, under "base repo: astropy/astropy" there is
   a drop-down list labeled "base branch: master".  You can click on this
   drop-down and instead select the bug fix branch ("v0.2.x" for example). Then
   GitHub will instead compare your fix against that branch, and merge into
   that branch when the PR is accepted.

Preparing the bug fix branch for release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two primary steps that need to be taken before creating a bug fix
release (the rest of the procedure is the same as any other release as
described in the release procedure below).

1. Any existing fixes to the issues assigned to the current bug fix release
   milestone, or labeled with the relevant "backport-x.y.z" label must be
   merged into the bug fix branch.

2. The Astropy changelog must be updated to list all issues--especially
   user-visible issues--fixed for the current release.  The changelog should
   be updated in the master branch, and then merged into the bug fix branch.

To aid in this process there is a script called ``suggest_backports.py`` at
https://gist.github.com/iguananaut/4497178.  The script is not perfect and
still needs a little work, but it will get most of the work done.  For example,
if the current bug fix branch is called 'v0.2.x' run it like so::

    $ suggest_backports.py astropy astropy v0.2.x -f backport.sh

This will search GitHub for all issues assigned to the next bug fix release
milestone that's associated with the given bug fix branch ('v0.2.2' for
example), find the commits that fix those issues, and will generate a shell
script called ``backport.sh`` containing all the ``git cherry-pick`` commands
to backport all those fixes.

The ``suggest_backports.py`` script will typically take a couple minutes to
run, but once it's done simply execute the generated script from within your
local clone of the Astropy repository::

    $ ./backport.sh

This will checkout the appropriate bug fix branch ('v0.2.x' in this example),
do a ``git pull upstream v0.2.x`` to make sure it's up to date, and then start
doing cherry-picks into the bug fix branch.

.. note::

    As discussed earlier, cherry-pick may result in merge conflicts.  If this
    occurs, the ``backport.sh`` script will exit and the conflict should be
    resolved manually, followed by running ``git commit``.  To resume the 
    ``backport.sh`` script after the merge conflict, it is currently necessary
    to edit the script to either remove or comment out the ``git cherry-pick``
    commands that already ran successfully.

    The author of the script hopes to improve it in the future to add
    ``git rebase`` like functionality, such that running
    ``backport.sh --continue`` or ``backport.sh --skip`` will be possible in
    such cases.

.. warning::

    It has also been noted that the ``suggest_backports.py`` script is not
    perfect, and can either miss issues that need to be backported, and in some
    cases can report false positives.

    It's always a good idea before finalizing a bug fix release to look on
    GitHub through the list of closed issues in the release milestone and check
    that each one has a fix in the bug fix branch.  Usually a quick way to do
    this is for each issue to run::

        $ git log --oneline <bugfix-branch> | grep #<issue>

    Most fixes will mention their related issue in the commit message, so this
    tends to be pretty reliable.  Some issues won't show up in the commit log,
    however, as their fix is in a separate pull request.  Usually GitHub makes
    this clear by cross-referencing the issue with its PR.  A future version
    of the ``suggest_backports.py`` script will perform this check
    automatically.

Finally, not all issues assigned to a release milestone need to be fixed before
making that release.  Usually, in the interest of getting a release with
existing fixes out within some schedule, it's best to triage issues that won't
be fixed soon to a new release milestone.  If the upcoming bug fix release is
'v0.2.2', then go ahead and create a 'v0.2.3' milestone and reassign to it any
issues that you don't expect to be fixed in time for 'v0.2.2'.


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

 5. Be sure you're the "master" branch or, for a bug fix release, on the
    appropriate bug fix branch.  For example, if releasing version 0.2.2 make
    sure to::
    
        $ git checkout v0.2.x

 6. Now install Astropy into the virtualenv::

        $ python setup.py install

    This is necessary for two reasons.  First, the entry points for the
    releaser scripts need to be available, and these are in the Astropy
    package. Second, the build process will generate .c files from the
    Cython .pyx files, and the .c files are necessary for the source
    distribution.

 7. Install zest.releaser into the virtualenv; use ``--upgrade --force`` to
    ensure that the latest version is installed in the virtualenv (if you're
    running a csh variant make sure to run ``rehash`` afterwards too)::

        $ pip install zest.releaser --upgrade --force

 8. Ensure that all changes to the code have been committed, then start the
    release by running::

        $ fullrelease

 9. You will be asked to enter the version to be released.  Press enter to
    accept the default (which will normally be correct) or enter a specific
    version string.  A diff will then be shown of CHANGES.rst and setup.py
    showing that a release date has been added to the changelog, and that the
    version has been updated in setup.py.  Enter 'Y' when asked to commit
    these changes.

 10. You will then be shown the command that will be run to tag the release.
     Enter 'Y' to confirm and run the command.

 11. When asked "Check out the tag (for tweaks or pypi/distutils server
     upload)" enter 'N': zest.releaser does not offer enough control yet over
     how the register and upload are performed so we will do this manually
     until the release scripts have been improved.

 12. You will be asked to enter a new development version.  Normally the next
     logical version will be selected--press enter to accept the default, or
     enter a specific version string.  Do not add ".dev" to the version, as
     this will be appended automatically (ignore the message that says ".dev0
     will be appended"--it will actually be ".dev" without the 0).  For
     example, if the just-released version was "0.1" the default next version
     will be "0.2".  If we want the next version to be, say "1.0" then that
     must be entered manually.

 13. You will be shown a diff of CHANGES.rst showing that a new section has
     been added for the new development version, and showing that the version
     has been updated in setup.py.  Enter 'Y' to commit these changes.

 14. When asked to push the changes to a remote repository, enter 'Y'.  This
     should complete the portion of the process that's automated at this point.

 15. Check out the tag of the released version.  For example::

         $ git checkout v0.1

 16. Create the source distribution by doing::

         $ python setup.py sdist

     Copy the produced ``.tar.gz`` somewhere and verify that you can unpack it,
     build it, and get all the tests to pass.  It would be best to create a new
     virtualenv in which to do this.

 17. Register the release on PyPI with::

         $ python setup.py register

 18. Upload the source distribution to PyPI; this is preceded by re-running
     the sdist command, which is necessary for the upload command to know
     which distribution to upload::

         $ python setup.py sdist upload

 19. Update the website to reflect the fact there is now a stable release.

 20. Update Readthedocs so that it builds docs for the corresponding github tag,
     and set the default page to the new release.

 21. If this was a major/minor release (not a bug fix release) create a bug fix
     branch for this line of release.  That is, if the version just released
     was "v<major>.<minor>.0", create bug fix branch with the name
     "v<major>.<minor>.x".  Starting from the commit tagged as the release,
     just checkout a new branch and push it to the remote server.  For example,
     after releasing version 0.3, do::

         $ git checkout -b v0.3.x

     Then edit ``setup.py`` so that the ``VERSION`` variable is
     ``'0.3.1.dev'``, and commit that change. Then, do::

         $ git push upstream v0.3.x

    .. note::

        You may need to replace ``upstream`` here with ``astropy`` or
        whatever remote name you use for the main astropy repository.

     The purpose of this branch is for creating bug fix releases like "0.3.1"
     and "0.3.2", while allowing development of new features to continue in
     the master branch.  Only changesets that fix bugs without making
     significant API changes should be merged to the bug fix branches.

 22. Create a bug fix label on GitHub; this should have the same name as the
     just created bug fix branch prepended with "backport-".  For the previous
     example this would be "backport-0.3.x"  This label should be applied to
     all issues that should be backported to the bug fix branch.  Also create a
     milestone for the next bug fix release if it hasn't been made already.


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


.. _signed tags: http://git-scm.com/book/en/Git-Basics-Tagging#Signed-Tags
.. _zest.releaser: http://pypi.python.org/pypi/zest.releaser
.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _cython: http://www.cython.org/
