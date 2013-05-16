==================
Release Procedures
==================

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
-------------------------------------------

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
^^^^^^^^^^^^
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


Release Procedure
-----------------

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


Creating a MacOS X Installer on a DMG
-------------------------------------

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
