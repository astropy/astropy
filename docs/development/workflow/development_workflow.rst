.. _development-workflow:

*******************************
How to make a code contribution
*******************************

This document outlines the process for contributing code to the Astropy
project.

**Already experienced with git? Contributed before?** Jump right to
:ref:`astropy-git`.

Pre-requisites
**************

Before following the steps in this document you need:

* an account on `GitHub`_
* a local copy of the astropy source. Instructions for doing that, including the
  basics you need for setting up git and GitHub, are at :ref:`get_devel`.

Strongly Recommended, but not required
**************************************

You cannot easily work on the development version of astropy in a python
environment in which you also use the stable version. It can be done |emdash|
but can only be done *successfully* if you always remember whether the
development version or stable version is the active one.

:ref:`virtual_envs` offer a better solution and take only a few minutes to set
up. It is well worth your time.

Not sure what your first contribution should be? Take a look at the `Astropy
issue list`_ and grab one labeled `"package-novice" <https://github.com/astropy/astropy/issues?q=is%3Aissue+is%3Aopen+label%3Apackage-novice>`_.
These issues are the most accessible ones if you are not familiar with the
Astropy source code. Issues labeled as `"effort-low" <https://github.com/astropy/astropy/issues?q=is%3Aissue+is%3Aopen+label%3Aeffort-low>`_
are expected to take a few hours (at most) to address, while the
`"effort-medium" <https://github.com/astropy/astropy/issues?q=is%3Aissue+is%3Aopen+label%3Aeffort-medium>`_
ones may take a few days. The developers are friendly and want you to help, so
don't be shy about asking questions on the |astropy-dev mailing list|.

New to `git`_?
**************

Some `git`_ resources
=====================

If you have never used git or have limited experience with it, take a few
minutes to look at `Git Basics`_, part of a much longer `git book`_.

In practice, you need only a handful of `git`_ commands to make contributions
to Astropy. There is a more extensive list of :ref:`git-resources` if you
want more background.

Double check your setup
=======================

Before going further, make sure you have set up astropy as described in
:ref:`get_devel`.

In a terminal window, change directory to the one containing your clone of
Astropy. Then, run ``git remote``; the output should look something like this::

    astropy
    origin

If that works, also run ``git fetch --all``. If it runs without errors then
your installation is working and you have a complete list of all branches in
your clone, ``origin`` and ``astropy``.

About names in `git`_
=====================

`git`_ is designed to be a *distributed* version control system. Each clone of
a repository is, itself, a repository. That can lead to some confusion,
especially for the branch called ``main``. If you list all of the branches
your clone of git knows about with ``git branch -a`` you will see there are
*three* different branches called ``main``::

    * main                              # this is main in your local repo
    remotes/astropy/main                # the official development branch of Astropy
    remotes/origin/main                 # main on your fork of Astropy on GitHub

The naming scheme used by `git`_ will also be used here. A plain branch name,
like ``main`` means a branch in your local copy of Astropy. A branch on a
remote, like ``astropy`` , is labeled by that remote, ``astropy/main``.

This duplication of names can get very confusing when working with pull
requests, especially when the official main branch, ``astropy/main``,
changes due to other contributions before your contributions are merged in.
As a result, you should never do any work in your main
branch, ``main``. Always work on a branch instead.

Essential `git`_ commands
=========================

A full `git`_ tutorial is beyond the scope of this document but this list
describes the few ``git`` commands you are likely to encounter in contributing
to Astropy:

* ``git fetch`` gets the latest development version of Astropy, which you will
  use as the basis for making your changes.
* ``git branch`` makes a logically separate copy of Astropy to keep track of
  your changes.
* ``git add`` stages files you have changed or created for addition to `git`_.
* ``git commit`` adds your staged changes to the repository.
* ``git push`` copies the changes you committed to GitHub
* ``git status`` to see a list of files that have been modified or created.

.. note::
    A good graphical interface to git makes some of these steps much
    easier.
    You might find this
    `list of GUI Clients <https://git-scm.com/downloads/guis/>`_ to be helpful.

If something goes wrong
=======================

`git`_ provides a number of ways to recover from errors. If you end up making a
`git`_ mistake, do not hesitate to ask for help. An additional resource that
walks you through recovering from `git`_ mistakes is the
`git choose-your-own-adventure`_.

.. _astropy-git:

Astropy Guidelines for `git`_
*****************************

.. note::
    It is strongly suggested that you automate the code-style checks using the
    provided pre-commit hook, see :ref:`pre-commit` below for details.

* Don't use your ``main`` branch for anything. Consider :ref:`delete-main`.
* Make a new branch, called a *feature branch*, for each separable set of
  changes: "one task, one branch" (`ipython git workflow`_).
* Start that new *feature branch* from the most current development version
  of astropy (instructions are below).
* Name your branch for the purpose of the changes, for example
  ``bugfix-for-issue-14`` or ``refactor-database-code``.
* Make frequent commits, and always include a commit message. Each commit
  should represent one logical set of changes.
* Ask on the |astropy-dev mailing list| if you get stuck.
* Never merge changes from ``astropy/main`` into your feature branch. If
  changes in the development version require changes to our code you can
  :ref:`rebase`.
* If you need to edit `.mailmap <https://git-scm.com/docs/gitmailmap>`_ and
  know how to do it then you can open a pull request for that. Please run
  `git shortlog -es <https://git-scm.com/docs/git-shortlog>`_ locally
  first with your changes to make sure your
  edit is correct, and you only appear in the list once.

In addition, there is a `git`_ naming convention used in this
document:

* Name the remote that is the primary Astropy repository
  ``astropy``; in prior versions of this documentation it was referred to as
  ``upstream``.


.. _pre-commit:

Pre-commit
**********

All of the coding style checks described in :ref:`code-style` can be performed automatically
when you make a git commit using our provided `pre-commit hook <https://pre-commit.com/>`_
for git, for more information see
`git hooks <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks#_git_hooks>`_.
We encourage you to setup and use these hooks to ensure that your code always meets
our coding style standards. The easiest way to do this is by installing ``pre-commit``::

    python -m pip install pre-commit

Or if you prefer |conda|::

    conda install pre-commit

This next command needs be done by installing ``pre-commit`` in the **root**
of your astropy repository by running::

    pre-commit install

For more detailed instructions on installing ``pre-commit``, see the
`install guide <https://pre-commit.com/#install>`_. Once this installation is
complete, all the coding style checks will be run each time you commit and the
necessary changes will automatically be applied to your code if possible.

.. note::
  The changes made by ``pre-commit`` will not be automatically staged, so you
  will need to review and re-stage any files that ``pre-commit`` has changed.

In general, git will not allow you to commit until the ``pre-commit`` hook has
run successfully. If you need to make a commit which fails the ``pre-commit`` checks,
you can skip these checks by running::

  git commit --no-verify

If you do not want to use ``pre-commit`` as part of your git workflow, you can
still run the checks manually (see, :ref:`code-style`) using::

  tox -e codestyle

Or this will run whether you did ``pre-commit install`` or not::

  pre-commit run

Again, this will automatically apply the necessary changes to your code if possible.

.. note::
  Once you have made a pull-request the ``pre-commit.ci`` bot is available to assist
  you with fixing any issues with your code style, see :ref:`pre-commit_bot` for details
  on how to use this bot.

Workflow
********

These, conceptually, are the steps you will follow in contributing to Astropy:

#. :ref:`fetch-latest`
#. :ref:`make-feature-branch`; you will make your changes on this branch.
#. :ref:`install-branch`
#. Follow :ref:`edit-flow` to write/edit/document/test code - make
   frequent, small commits.
#. :ref:`add-changelog`
#. :ref:`push-to-github`
#. From GitHub, :ref:`pull-request` to let the Astropy maintainers know
   you have contributions to review.
#. :ref:`revise and push` in response to comments on the pull
   request. Pushing those changes to GitHub automatically updates the
   pull request.

This way of working helps to keep work well organized, with readable history.
This in turn makes it easier for project maintainers (that might be you) to
see what you've done, and why you did it.

A worked example that follows these steps for fixing an Astropy issue is at
:ref:`astropy-fix-example`.

Some additional topics related to `git`_ are in :ref:`additional-git`.

.. _delete-main:

Deleting your main branch
=========================

It may sound strange, but deleting your own ``main`` branch can help reduce
confusion about which branch you are on.  See `deleting main on github`_ for
details.

.. _fetch-latest:

Fetch the latest Astropy
************************

From time to time you should fetch the development version (i.e., Astropy
``astropy/main``) changes from GitHub::

   git fetch astropy --tags

This will pull down any commits you don't have, and set the remote branches to
point to the latest commit. For example, 'trunk' is the branch referred to by
``astropy/main``, and if there have been commits since
you last checked, ``astropy/main`` will change after you do the fetch.

.. _make-feature-branch:

Make a new feature branch
*************************

Make the new branch
===================

When you are ready to make some changes to the code, you should start a new
branch. Branches that are for a collection of related edits are often called
'feature branches'.

Making a new branch for each set of related changes will make it easier for
someone reviewing your branch to see what you are doing.

Choose an informative name for the branch to remind yourself and the rest of us
what the changes in the branch are for. Branch names like ``add-ability-to-fly``
or ``bugfix-for-issue-42`` clearly describe the purpose of the branch.

Always make your branch from ``astropy/main`` so that you are basing your
changes on the latest version of Astropy::

    # Update the mirror of trunk
    git fetch astropy --tags

    # Make new feature branch starting at astropy/main
    git branch my-new-feature astropy/main
    git checkout my-new-feature

Connect the branch to GitHub
============================

At this point you have made and checked out a new branch, but `git`_ does not
know it should be connected to your fork on GitHub. You need that connection
for your proposed changes to be managed by the Astropy maintainers on GitHub.

The most convenient way for connecting your local branch to GitHub is to `git
push`_ this new branch up to your GitHub repo with the ``--set-upstream``
option::

   git push --set-upstream origin my-new-feature

From now on git will know that ``my-new-feature`` is related to the
``my-new-feature`` branch in your GitHub fork of Astropy.

You will still need to :ref:`push-to-github` periodically. The
setup in this section will make that easier because any following pushes of
this branch can be performed without having to write out the remote and branch
names, but it never hurts to be explicit in typing out the commands.

.. _install-branch:

Install your branch
*******************

Ideally you should set up a Python virtual environment just for this fix;
instructions for doing to are at :ref:`virtual_envs`. Doing so ensures you
will not corrupt your main ``astropy`` install and makes it very easy to recover
from mistakes, and thus, is recommended before you proceed.

Assuming you have set up and activated this virtual environment, you need to
install the version of ``astropy`` you are working on into it. Do that with:

.. code-block:: bash

    python -m pip install --editable ".[test]"

This will install ``astropy`` itself, along with a few packages which will be
useful for testing the changes you will make down the road.

For more details on building ``astropy`` from source, see
:ref:`dev-build-astropy-subpkg`.

.. _edit-flow:

The editing workflow
********************

Conceptually, you will:

#. Make changes to one or more files and/or add a new file.
#. Check that your changes do not break existing code.
#. Add documentation to your code and, as appropriate, to the Astropy
   documentation.
#. Ideally, also make sure your changes do not break the documentation.
#. Add tests of the code you contribute.
#. Commit your changes in `git`_.
#. Repeat as necessary.


In more detail
==============

#. Make some changes to one or more files. You should follow the Astropy
   :ref:`code-guide`. Each logical set of changes should be treated as one
   commit. For example, if you are fixing a known bug in Astropy and notice
   a different bug while implementing your fix, implement the fix to that new
   bug as a different set of changes.

#. Test that your changes do not lead to *regressions*, i.e. that your
   changes do not break existing code, by running the Astropy tests. You can
   run all of the Astropy tests from ipython with::

     import astropy
     astropy.test()

   If your change involves only a small part of Astropy, e.g. Time, you can
   run just those tests::

     import astropy
     astropy.test(package='time')

   Tests can also be run from the command line while in the package
   root directory, e.g.::

     pytest

   To run the tests in only a single package, e.g. Time, you can do::

     pytest -P time

   For more details on running tests, please see :ref:`testing-guidelines`.

#. Make sure your code includes appropriate docstrings, in the
   `Numpydoc format`_.
   If appropriate, as when you are adding a new feature,
   you should update the appropriate documentation in the ``docs`` directory;
   a detailed description is in :ref:`documentation-guidelines`.

#. If you have sphinx installed, you can also check that
   the documentation builds and looks correct by running, from the
   ``astropy`` directory::

     cd docs
     make html

   The last line should just state ``build succeeded``, and should not mention
   any warnings.  (For more details, see :ref:`documentation-guidelines`.)

#. Add tests of your new code, if appropriate. Some changes (e.g. to
   documentation) do not need tests. Detailed instructions are at
   :ref:`writing-tests`, but if you have no experience writing tests or
   with the |pytest| testing framework submit your changes without adding
   tests, but mention in the pull request that you have not written tests.
   An example of writing a test is in
   :ref:`astropy-fix-add-tests`.

#. Stage your changes using ``git add`` and commit them using ``git commit``.
   An example of doing that, based on the fix for an actual Astropy issue, is
   at :ref:`astropy-fix-example`.

   .. note::
        Make your `git`_ commit messages short and descriptive. If a commit
        fixes an issue, include, on the second or later line of the commit
        message, the issue number in the commit message, like this:
        ``Closes #123``. Doing so will automatically close the issue when the
        pull request is accepted.

#. Some modifications require more than one commit; if in doubt, break
   your changes into a few, smaller, commits rather than one large commit
   that does many things at once. Repeat the steps above as necessary!

.. _add-changelog:

Add a changelog entry
*********************

Add a changelog fragment briefly describing the change you made by creating
a new file in ``docs/changes/<sub-package>/``. The file should be named like
``<PULL REQUEST>.<TYPE>.rst``, where ``<PULL REQUEST>`` is a pull request
number, and ``<TYPE>`` is one of:

* ``feature``: New feature.
* ``api``: API change.
* ``bugfix``: Bug fix.
* ``other``: Other changes and additions.

An example entry, for the changes in `PR 1845
<https://github.com/astropy/astropy/pull/1845>`_, the file would be
``docs/changes/wcs/1845.bugfix.rst`` and would contain::

    ``astropy.wcs.Wcs.printwcs`` will no longer warn that ``cdelt`` is
    being ignored when none was present in the FITS file.

If you are opening a new pull request, you may not know its number yet, but you
can add it *after* you make the pull request. If you're not sure where to
put the changelog entry, wait at least until a maintainer has reviewed your
PR and assigned it to a milestone.

When writing changelog entries, do not attempt to make API reference links
by using single-backticks.  This is because the changelog (in its current
format) runs for the history of the project, and API references you make today
may not be valid in a future version of Astropy.  However, use of
double-backticks for monospace rendering of module/class/function/argument
names and the like is encouraged.

.. _push-to-github:

Push your changes to GitHub
***************************

To push your changes on a branch named ``my-new-feature`` to your
GitHub fork of ``astropy`` with the same branch name::

    git push origin my-new-feature

.. _pull-request:

Ask for your changes to be reviewed
***********************************

A *pull request* on GitHub is a request to merge the changes you have made into
another repository.

When you are ready to ask for someone to review your code and consider merging
it into Astropy:

#. Go to the URL of your fork of Astropy, e.g.,
   ``https://github.com/your-user-name/astropy``.

#. Use the 'Switch Branches' dropdown menu to select the branch with your
   changes:

   .. image:: branch_dropdown.png

#. Click on the 'Pull request' button:

   .. image:: pull_button.png

   Enter a title for the set of changes, and some explanation of what you've
   done. If there is anything you'd like particular attention for, like a
   complicated change or some code you are not happy with, add the details
   here.

   If you don't think your request is ready to be merged, just say so in your
   pull request message.  This is still a good way to start a preliminary
   code review.

   You may also opt to open a work-in-progress pull request.
   If you do so, instead of clicking "Create pull request", click on the small
   down arrow next to it and select "Create draft pull request". This will let
   the maintainers know that your work is not ready for a full review nor to be
   merged yet. In addition, if your commits are not ready for CI testing, you
   should also use ``[ci skip]`` or ``[skip ci]`` directive in your commit message.
   For usage of pre-commit hooks and directives, see :ref:`pre-commit` and :ref:`pre-commit_bot`.

.. _revise and push:

Revise and push as necessary
****************************

You may be asked to make changes in the discussion of the pull request. Make
those changes in your local copy, commit them to your local repo and push them
to GitHub. GitHub will automatically update your pull request.

.. _no git pull:

Do Not Create a Merge Commit
****************************

If your branch associated with the pull request falls behind the ``main``
branch of https://github.com/astropy/astropy, GitHub might offer you the option
to catch up or resolve conflicts via its web interface, but do not use this. Using
the web interface might create a "merge commit" in your commit history, which is
undesirable, as a "merge commit" can introduce maintenance overhead for the
release manager as well as undesirable branch structure complexity. Do not use the ``git pull`` command either.

Instead, in your local checkout, do a ``fetch`` and then a ``rebase``, and
resolve conflicts as necessary. See :ref:`rebase` and :ref:`howto_rebase`
for further information.

.. _rebase:

Rebase if necessary
*******************

Rebasing means taking your changes and applying them to the latest
version of the ``main`` branch of the official ``astropy`` repository as though that was the
version you had originally branched from. Each individual commit remains
visible, but with new commit hashes.

When to rebase
==============

Pull requests **must** be rebased (but not necessarily squashed to a single
commit) if at least one of the following conditions is met:

* There are conflicts with ``main``.
* There are commits in ``main`` (after the PR branch point) needed for continuous
  integration tests to run correctly.
* There are merge commits from ``astropy/main`` in the PR commit history (merge
  commits from PRs to the user's fork are fine).

.. _howto_rebase:

How to rebase
=============

It is easier to make mistakes rebasing than other areas of `git`_, so before you
start make a branch to serve as a backup copy of your work::

    git branch tmp my-new-feature # make temporary branch--will be deleted later

After altering the history, e.g., with ``git rebase``, a normal ``git push``
is prevented, and a ``--force`` option will be required.

.. warning::

    Do not update your branch with ``git pull``. Pulling changes from
    ``astropy/main`` includes merging the branches, which combines them in a
    way that preserves the commit history of both. The purpose of rebasing is
    rewriting the commit history of your branch, not preserving it.

Behind the scenes, `git`_ is deleting the changes and branch you made, making the
changes others made to the development branch of Astropy, then re-making your
branch from the development branch and applying your changes to your branch.

First, fetch the latest development astropy and go to your branch of interest::

    git fetch astropy main
    git checkout my-new-feature

Now, do the rebase::

    git rebase astropy/main

You are more likely to run into *conflicts* here — places where the changes you made
conflict with changes that someone else made — than anywhere else. Ask for help if you
need it. Instructions are available on how to `resolve merge conflicts after a Git
rebase
<https://help.github.com/en/articles/resolving-merge-conflicts-after-a-git-rebase>`_.

.. _squash-if-necessary:

Squash if necessary
*******************

Squashing refers to combining multiple commits into a single commit. This can be done
using the ``git rebase`` command or via :ref:`github-squash-and-merge`.

An Astropy maintainer will be happy to guide you through this process.

When to squash
==============

If a pull request contains commits with large (approximately > 10KB) intermediate
changes which are ultimately removed or modified in the PR, the intermediate diffs
**must** be squashed. An example is adding a large data file in one commit and then
removing it in a subsequent commit. The goal is to avoid an unnecessary increase in the
size of the ``astropy`` repository.

Squashing commits is **encouraged** when there are numerous commits which do not add
value to the commit history. Most small to medium pull requests can be done with a few
commits. Some large or intricate pull requests may require more commits to capture the
logical progression of the changes.

In general, commits that reflect a specific atomic change (e.g., "Fixed bug revealed by
tests for feature X") may have value for the history.

Commits that are good candidates for squashing include but not limited to:

* Content that gets removed later, most commonly changes in the implementation or
  temporary debugging code or, especially, data files (see above).
* Non-specific commits; e.g., "Added more code and fixed stuff."
* Fixes of typos, linting fixes, or other inconsequential changes.
* Commit messages that violate the code of conduct.

.. _howto_squash:

How to squash
=============

In many cases, squashing can be done by a maintainer using the :ref:`github-squash-and-merge` button on
the GitHub pull request page. If this is not possible, we typically squash using `git
rebase --interactive <https://git-scm.com/docs/git-rebase#_interactive_mode>`_. In
particular, you can rebase and squash within the existing branch using::

  git fetch astropy
  git rebase -i astropy/main

The last command will open an editor with all your commits, allowing you to
squash several commits together, rename them, etc. Helpfully, the file you are
editing has the instructions on what to do.

.. _howto_push_force:

How to push
***********

After using ``git rebase`` you will still need to push your changes to
GitHub so that they are visible to others and the pull request can be
updated.  Use of a simple ``git push`` will be prevented because of the
changed history, and will need to be manually overridden using::

    git push origin my-new-feature --force

If you run into any problems, do not hesitate to ask. A more detailed conceptual
discussing of rebasing is at :ref:`rebase-on-trunk`.

Once the modifications and new git history are successfully pushed to GitHub you
can delete any backup branches that may have been created::

    git branch -D tmp

.. include:: links.inc

.. _Git Basics: https://git-scm.com/book/en/Getting-Started-Git-Basics
.. _git book: https://git-scm.com/book/
.. _Astropy issue list: https://github.com/astropy/astropy/issues
.. _git choose-your-own-adventure: http://sethrobertson.github.io/GitFixUm/fixup.html
.. _numpydoc format: https://numpydoc.readthedocs.io/en/latest/format.html
