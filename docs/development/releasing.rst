==================
Release Procedures
==================

The current release procedure for Astropy involves a combination of an
automated release script and some manual steps.  Future versions will automate
more of the process, if not all.


.. _release-procedure:

Release Procedure
-----------------

This is the standard release procedure for releasing Astropy (or affiliated
packages that use the full bugfix/maintainence branch approach.)

#. (Optional) You may want to set up a clean environment to build the release.
   For more on setting up virtual environments, see :ref:`virtual_envs`, but
   for the sake of example we will assume you're using `Anaconda`_. This is not
   necessary if you know your normal python environment has what you need, but
   you might want to do something like this for safety's sake::

      $ conda create -n astropy_release_build_v<verson> astropy
      $ source activate astropy_release_build_v<verson>
      $ conda uninstall astropy  # still keeps the dependencies
      $ pip install -r pip-requirements-dev  # any that might be left over

#. Before doing a release of Astropy, you may need to do a release of
   `astropy-helpers`_.  This is not always necessary if there are no significant
   changes in the helpers.  See the section :ref:`helpers-release-info` for
   more on this.

#. Make sure that Travis and any other continuous integration is passing for
   the branch you're going to release.  You may also want to locally run the
   tests in ``remote-data`` mode, as those are not necessarily run
   automatically::

      $ python setup.py test --remote-data

#. Ensure you have a GPG key pair available for when git needs to sign the
   tag you create for the release.  See :ref:`key-signing-info` for more on
   this.

#. Obtain a *clean* version of the Astropy repository.  That is, one
   where you don't have any intermediate build files.  Either use a fresh
   ``git clone`` or do ``git clean -dfx``.

#. Be sure you're on the branch appropriate for the version you're about to
   release.  For example, if releasing version 0.2.2 make sure to::

      $ git checkout v0.2.x

#. Edit the ``setup.py`` file by removing the ``".dev"`` at the end of the
   ``VERSION`` string.

#. Edit the ``CHANGES.rst`` file by changing the date for the version you are
   about to release from "unreleased" to today's date.  Also be sure to remove
   any sections of the changelog for that version that have no entries.

#. Add the changes to ``CHANGES.rst`` and ``setup.py``::

      $ git add CHANGES.rst setup.py

   and commit with message::

      $ git commit -m "Preparing release v<version>"

#. Tag the commit with ``v<version>``, being certain to sign the tag with the
   ``-s`` option::

      $ git tag -s v<version>

#. Edit the ``VERSION`` in ``setup.py`` to be the next version number, but with
   a ``.dev`` suffix at the end (E.g., ``0.2.3.dev``).

#. Also update the ``CHANGES.rst`` file with a new section for the next version.
   You will likely want to use the ``add_to_changelog.py`` script in the
   `astropy-tools`_ repository for this.

#. Add the changes to ``CHANGES.rst`` and ``setup.py``::

      $ git add CHANGES.rst setup.py

   and commit with message::

      $ git commit -m "Back to development: v<next_version>"

#. Now go back and check out the tag of the released version with
   ``git checkout v<version>``.  For example::

      $ git checkout v0.2.2

   Don't forget to remove any non-commited files with::

      $  git clean -dfx

#. Create the source distribution by doing::

         $ python setup.py build sdist

     .. note::

         In the future, the ``build`` command may run automatically as a
         prerequisite for ``sdist``.  But for now, make sure to run it
         whenever running ``sdist`` to ensure that all Cython sources and
         other generated files are built.

#. Run the tests in an environment that mocks up a "typical user" scenario.
   This is not strictly necessary because you ran the tests above, but
   it can sometimes be useful to catch subtle bugs that might come from you
   using a customized developer environment.  For more on setting up virtual
   environments, see :ref:`virtual_envs`, but for the sake of example we will
   assume you're using `Anaconda`_. Do::

      $ conda create -n astropy_release_test_v<verson> numpy
      $ source activate astropy_release_test_v<verson>
      $ cd dist
      $ mkdir test
      $ pip install astropy-<version>.tar.gz
      $ python -c 'import astropy; astropy.test(remote_data=True)'
      $ source deactivate
      $ cd <back to the source directory>

   Assuming everything went smoothly with the last step, we are ready to proceed
   with the release.  If not, you'll have to back up to before starting this
   procedure and start over when things are fixed.

#. Register the release on PyPI with::

      $ python setup.py register

#. Upload the source distribution to PyPI; this is preceded by re-running
   the sdist command, which is necessary for the upload command to know
   which distribution to upload::

      $ python setup.py build sdist upload --sign

#. Go to https://pypi.python.org/pypi?:action=pkg_edit&name=astropy
   and ensure that only the most recent releases in each actively maintained
   release line are *not* marked hidden.  For example, if v0.2.2 was
   just released, v0.2.1 should be hidden.  This is so that users only find
   the latest bugfix releases.

   Do not enabled "Auto-hide old releases" as that may hide bugfix releases
   from older release lines that we may still want to make available.

#. Push up all these changes to the astropy repository::

      $ git push --tags origin v<version>

   .. note::

      You may need to replace ``origin`` here with ``astropy`` or
      whatever remote name you use for the main astropy repository.

#. If this is a release of the current release (not an LTS), update the
   "stable" branch to point to the new release::

      $ git checkout stable
      $ git reset --hard v<version>
      $ git push origin stable --force

#. Update Readthedocs so that it builds docs for the corresponding github tag.
   Also verify that the ``stable`` Readthedocs version builds correctly for
   the new version (it should trigger automatically once you've done the
   previous step.)

   When releasing a patch release, also set the previous version in the
   release history to "protected".  For example when releasing v1.1.2, set
   v1.1.1 to "protected".  This prevents the previous releases from
   cluttering the list of versions that users see in the version dropdown
   (the previous versions are still accessible by their URL though).

#. Update the Astropy web site by editing the ``index.html`` page at
   https://github.com/astropy/astropy.github.com by changing the "current
   version" link and/or updating the list of older versions if this is an LTS
   bugfix or a new major version.

#. In the astropy *master* branch (not just the changelog), be sure to update
   the ``CHANGES.rst`` to reflect the date of the release you just performed and
   to include the new section of the changelog.  Often the easiest way to do
   this is to use ``git cherry-pick`` on the release commit from above, *but* if
   that method is used, be sure to amend that commit and not change the version
   in ``setup.py``. If you aren't sure how to do this, you might be better off
   just copying-and-pasting the relevant parts of the maintainence branch's
   ``CHANGES.rst`` into master.

Modifications for a beta/release candidate release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   For major releases with a lot of changes, we do beta and/or release candiates
   to have a chance to catch significant bugs before the true release.
   If the release you are performing is this kind of pre-release, some
   of the above steps need to be modified.

   The primary modifications to the release procedure are:

   * When entering the new version number, instead of just removing the
     ``.dev``, enter "0.2b1" or "0.2rc1".  It is critical that you follow this
     numbering scheme (``x.yb#`` or ``x.y.zrc#``), as it will ensure the release
     is ordered "before" the main release by various automated tools, and also
     tells PyPI that this is a "pre-release".
   * Do not do step #21 or later, as those are tasks for an actual release.


Performing a Feature Freeze/Branching new Major Versions
--------------------------------------------------------

As outlined in
`APE2 <https://github.com/astropy/astropy-APEs/blob/master/APE2.rst>`_, astropy
releases occur at regular intervals, but feature freezes occur well before the
actual release.  Feature freezes are also the time when the master branch's
development separates from the new major version's maintainence branch.  This
allows new development for the next major version to continue while the
soon-to-be-released version can focus on bug fixes and documentation updates.

The procedure for this is straightforward:

#. Make sure you're on master, and updated to the latest version from github::

      $ git checkout master
      $ git fetch origin
      $ git reset --hard origin/master

#. Create a new branch from master at the point you want the feature freeze to
   occur::

      $ git branch v<version>.x

#. Update the ``CHANGES.rst`` file with a new section at the very top for the
   next major version.  You will likely want to use the ``add_to_changelog.py``
   script in the `astropy-tools`_ repository for this.

#. Update the ``VERSION`` in ``setup.py`` to reflect the new major version. For
   example, if you are about to issue a feature freeze for version ``1.2``, you
   will want to set the new version to ``'1.3.dev'``.

#. Add the changes to ``CHANGES.rst`` and ``setup.py``::

      $ git add CHANGES.rst setup.py

   and commit with message::

      $ git commit -m "Next major version: <next_version>"

#. Push all of these changes up to github::

      $ git push origin v<version>.x:v<version>.x
      $ git push origin master:master

   .. note::

      You may need to replace ``origin`` here with ``astropy`` or
      whatever remote name you use for the main astropy repository.

#. On the github issue tracker, add a new milestone for the next major version.


Maintaining Bug Fix Releases
----------------------------

Astropy releases, as recommended for most Python projects, follows a
<major>.<minor>.<micro> version scheme, where the "micro" version is also
known as a "bug fix" release.  Bug fix releases should not change any user-
visible interfaces.  They should only fix bugs on the previous major/minor
release and may also refactor internal APIs or include omissions from previous
releases--that is, features that were documented to exist but were accidentally
left out of the previous release. They may also include changes to docstrings
that enhance clarity but do not describe new features (e.g., more examples,
typo fixes, etc).

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
``git commit`` and accept the default commit message.  If the fix being
cherry-picked has an associated changelog entry in a separate commit make
sure to backport that as well.

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
release. The rest of the procedure is the same as any other release as
described in :ref:`release-procedure` (although be sure to provide the
right version number).

1. Any existing fixes to the issues assigned to the current bug fix release
   milestone, or labeled with the relevant "backport-x.y.z" label must be
   merged into the bug fix branch.

2. The Astropy changelog must be updated to list all issues--especially
   user-visible issues--fixed for the current release.  The changelog should
   be updated in the master branch, and then merged into the bug fix branch.
   Most issues *should* already have changelog entries for them. But it's
   typical to forget this, so if doesn't exist yet please add one in
   the process of backporting.  See :ref:`changelog-format` for more details.

To aid in this process there is a `suggest_backports.py script in the astropy-tools repository <https://github.com/astropy/astropy-tools/blob/master/suggest_backports.py>`_.
The script is not perfect and still needs a little work, but it will get most of
the work done.  For example, if
the current bug fix branch is called 'v0.2.x' run it like so::

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


.. _helpers-release-info:

Coordinating Astropy and astropy-helpers Releases
-------------------------------------------------

A bit more initial effort is required for an Astropy release that has a
corresponding astropy-helpers release.  The main reason for this more complex
procedure is to allow the Astropy core to be tested aginst the new helpers
before anything is released.  Hence the following procedure should be added
to the beginning of the above procedure when this is required.

#. In `astropy-helpers`, create a new release branch "release-<version>".

#. Create the release commit (updating the version info and changelog) in that
   branch.

#. Push the release branch to github.

#. In astropy master (and/or the relevant maintainence branch), issue a PR
   updating the helpers to the commit described in the last step (i.e., the
   head of the astropy-helpers release branch).

#. Wait for Travis to run to ensure that helpers build works with Astropy.
   If it doesn't, back out the release and fix whatever the problem is before
   trying again.

#. Assuming it does succeed, finish the release of the helpers by doing this in
   the helpers repo::

      git checkout master
      git merge --no-ff release-<version>
      git tag -s "v<version>" -m "Tagging v<version>"

   and then adding one more commit updating back to the next dev version.

This way the commit of the helpers that is tagged as the release is the same
commit that the astropy_helpers submodule will be on when the PR to astropy
testing the release gets merged.


.. _key-signing-info:

Creating a GPG Signing Key and a Signed Tag
-------------------------------------------

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

Generating a public/private key pair
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


.. _signed tags: http://git-scm.com/book/en/Git-Basics-Tagging#Signed-Tags
.. _cython: http://www.cython.org/
.. _astropy-tools: https://github.com/astropy/astropy-tools
.. _Anaconda: http://conda.pydata.org/docs/
.. _astropy-helpers: https://github.com/astropy/astropy-helpers
