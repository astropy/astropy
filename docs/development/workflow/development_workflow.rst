.. _development-workflow:

===============================
How to make a code contribution
===============================

This document outlines the process for contributing code to the Astropy
project.

**Already experienced with git? Contributed before?** Jump right to
:ref:`astropy-git`.

Pre-requisites
==============

Before following the steps in this document you need:

+ an account on `GitHub`_
+ a local copy of the astropy source. Instructions for doing that, including the
  basics you need for setting up git and GitHub, are at :ref:`get_devel`.

Strongly Recommended, but not required
======================================

You cannot easily work on the development version of astropy in a python
environment in which you also use the stable version. It can be done |emdash|
but can only be done *successfully* if you always remember whether the
development version or stable version is the active one.

:ref:`virtual_envs` offer a better solution and take only a few minutes to set
up. It is well worth your time.

Not sure what your first contribution should be? Take a look at the `Astropy
issue list`_ and grab one labeled "package-novice". These issues are the
most accessible ones if you are not familiar with the Astropy source
code. Issues labeled as "effort-low" are expected to take a few hours (at
most) to address, while the "effort-medium" ones may take a few days. The
developers are friendly and want you to help, so don't be shy about asking
questions on the `astropy-dev mailing list`_.

New to `git`_?
==============

Some `git`_ resources
---------------------

If you have never used git or have limited experience with it, take a few
minutes to look at these resources:

* `Interactive tutorial`_ that runs in a browser
* `Git Basics`_, part of a much longer `git book`_.

In practice, you need only a handful of `git`_ commands to make contributions
to Astropy. There is a more extensive list of :ref:`git-resources` if you
want more background.

Double check your setup
-----------------------

Before going further, make sure you have set up astropy as described in
:ref:`get_devel`.

In a terminal window, change directory to the one containing your clone of
Astropy. Then, run ``git remote``; the output should look something like this::

    your-github-username
    astropy

If that works, also run ``git fetch --all``. If it runs without errors then
your installation is working and you have a complete list of all branches in
your clone, ``your-github-username`` and ``astropy``.

About names in `git`_
---------------------

`git`_ is designed to be a *distributed* version control system. Each clone of
a repository is, itself, a repository. That can lead to some confusion,
especially for the branch called ``master``. If you list all of the branches
your clone of git knows about with ``git branch -a`` you will see there are
*three* different branches called ``master``::

    * master                              # this is master in your local repo
    remotes/your-github-username/master   # master on your fork of Astropy on GitHub
    remotes/astropy/master                # the official development branch of Astropy

The naming scheme used by `git`_ will also be used here. A plain branch name,
like ``master`` means a branch in your local copy of Astropy. A branch on a
remote, like ``astropy`` , is labeled by that remote, ``astropy/master``.

This duplication of names can get very confusing for maintainers when trying
to merge code contributions into the official master branch,
``astropy/master``. As a result, you should never do any work in your master
branch, ``master``. Always work on a branch instead.

Essential `git`_ commands
-------------------------

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
    easier. Some options are described in :ref:`git_gui_options`.

If something goes wrong
-----------------------

`git`_ provides a number of ways to recover from errors. If you end up making a
`git`_ mistake, do not hesitate to ask for help. An additional resource that
walks you through recovering from `git`_ mistakes is the
`git choose-your-own-adventure`_.

.. _astropy-git:

Astropy Guidelines for `git`_
=============================

* Don't use your ``master`` branch for anything.
* Make a new branch, called a *feature branch*, for each separable set of
  changes: "one task, one branch" (`ipython git workflow`_).
* Start that new *feature branch* from the most current development version
  of astropy (instructions are below).
* Name your branch for the purpose of the changes, for example
  ``bugfix-for-issue-14`` or ``refactor-database-code``.
* Make frequent commits, and always include a commit message. Each commit
  should represent one logical set of changes.
* Ask on the `astropy-dev mailing list`_ if you get stuck.
* Never merge changes from ``astropy/master`` into your feature branch. If
  changes in the development version require changes to our code you can
  :ref:`rebase`.

In addition there are a couple of `git`_ naming conventions used in this
document:

* Change the name of the remote ``origin`` to ``your-github-username``.
* Name the remote that is the primary Astropy repository
  ``astropy``; in prior versions of this documentation it was referred to as
  ``upstream``.

Workflow
========

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

.. _fetch-latest:

Fetch the latest Astropy
========================

From time to time you should fetch the development version (i.e. Astropy
``astropy/master``) changes from GitHub::

   git fetch astropy

This will pull down any commits you don't have, and set the remote branches to
point to the latest commit. For example, 'trunk' is the branch referred to by
``astropy/master``, and if there have been commits since
you last checked, ``astropy/master`` will change after you do the fetch.

.. _make-feature-branch:

Make a new feature branch
=========================

Make the new branch
-------------------

When you are ready to make some changes to the code, you should start a new
branch. Branches that are for a collection of related edits are often called
'feature branches'.

Making a new branch for each set of related changes will make it easier for
someone reviewing your branch to see what you are doing.

Choose an informative name for the branch to remind yourself and the rest of us
what the changes in the branch are for. Branch names like ``add-ability-to-fly``
or ``buxfix-for-issue-42`` clearly describe the purpose of the branch.

Always make your branch from ``astropy/master`` so that you are basing your
changes on the latest version of Astropy::

    # Update the mirror of trunk
    git fetch astropy

    # Make new feature branch starting at astropy/master
    git branch my-new-feature astropy/master
    git checkout my-new-feature

Connect the branch to GitHub
----------------------------

At this point you have made and checked out a new branch, but `git`_ does not
know it should be connected to your fork on GitHub. You need that connection
for your proposed changes to be managed by the Astropy maintainers on GitHub.

To connect your local branch to GitHub, you `git push`_ this new branch up to
your GitHub repo with the ``--set-upstream`` option::

   git push --set-upstream your-github-username my-new-feature

From now on git will know that ``my-new-feature`` is related to the
``your-github-username/my-new-feature`` branch in your GitHub fork of Astropy.

You will still need to ``git push`` your changes to GitHub periodically. The
setup in this section will make that easier.

.. _install-branch:

Install your branch
===================

Ideally you should set up a python virtual environment just for this fix;
instructions for doing to are at :ref:`virtual_envs`. Doing so ensures you
will not corrupt your main astropy install and makes it very easy to recover
from mistakes.

Once you have activated that environment you need to install the version of
Astropy you are working on. Do that with:

.. code-block:: bash

    python setup.py develop  # typically python 2.x, not python 3

or:

.. code-block:: bash

    python3 setup.py install # python 3...
    # ...though python3 may be called python3.3 or just python,
    # depending on your system.

If you are using python 3 you will need to re-install after making changes to
the Astropy source code. Re-installing goes much faster than the initial install
because it typically does not require new compilation.

.. _edit-flow:

The editing workflow
====================

Conceptually, you will:

#. Make changes to one or more files and/or add a new file.
#. Check that your changes do not break existing code.
#. Add documentation to your code and, as appropriate, to the Astropy
   documentation.
#. Ideally, also make sure your changes do not break the documentation.
#. Add tests of the code you contribute.
#. Commit your changes in `git`_
#. Repeat as necessary.


In more detail
--------------

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
     astropy.test('time')

#. Make sure your code includes appropriate docstrings, described at
   :ref:`doc-rules`. If appropriate, as when you are adding a new feature,
   you should update the appropriate documentation in the ``docs`` directory;
   a detailed description is in :ref:`documentation-guidelines`.

#. If you have sphinx installed, you can also check that
   the documentation builds and looks correct by running, from the
   ``astropy`` directory::

     python setup.py build_sphinx

   The last line should just state ``build succeeded``, and should not mention
   any warnings.  (For more details, see :ref:`documentation-guidelines`.)

#. Add tests of your new code, if appropriate. Some changes (e.g. to
   documentation) do not need tests. Detailed instructions are at
   :ref:`testing-guidelines`, but if you have no experience writing tests or
   with the `py.test`_ testing framework submit your changes without adding
   tests, but mention in the pull request that you have not written tests. An
   example of writing a test is in :ref:`astropy-fix-example`.

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
=====================

Add an entry to the file ``CHANGES.rst`` briefly describing the change you
made. Include the pull request number if the change fixes an issue. An
example entry, for the changes which fixed
`issue 1845 <https://github.com/astropy/astropy/pull/1845>`_, is::

  - `astropy.wcs.Wcs.printwcs` will no longer warn that `cdelt` is
    being ignored when none was present in the FITS file. [#1845]

If the change is a new feature, rather than an existing issue, you will not be
able to put in the issue number until *after* you make the pull request.

.. _push-to-github:

Copy your changes to GitHub
===========================

This step is easy because of the way you created the feature branch. Just::

    git push

.. _pull-request:

Ask for your changes to be reviewed
===================================

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

.. _revise and push:

Revise and push as necessary
============================

You may be asked to make changes in the discussion of the pull request. Make
those changes in your local copy, commit them to your local repo and push them
to GitHub. GitHub will automatically update your pull request.

.. _rebase:

Rebase, but only if asked
=========================

Sometimes the maintainers of Astropy will ask you to *rebase* your changes
before they are merged into the main Astropy repository.

Conceptually, rebasing means taking your changes and applying them to the latest
version of the development branch of the official astropy as though that was the
version you had originally branched from.

Behind the scenes, `git`_ is deleting the changes and branch you made, making the
changes others made to the development branch of Astropy, then re-making your
branch from the development branch and applying your changes to your branch.
This results in re-writing the history of commits, which is why you should do it
only if asked.

It is easier to make mistakes rebasing than other areas of `git`_, so before you
start make a branch to serve as a backup copy of your work::

    git branch tmp my-new-feature # make temporary branch--will be deleted later

The actual rebasing is usually easy::

    git fetch astropy/master  # get the latest development astropy
    git rebase astropy/master my-new-feature

You are more likely to run into *conflicts* here--places where the changes you
made conflict with changes that someone else made--than anywhere else. Ask for
help if you need it.

After the rebase you need to push your changes to GitHub; you will need force
the push because `git`_ objects to re-writing the history of the repository
after you have pushed it somewhere::

    git push -f

If you run into any problems, do not hesitate to ask. A more detailed conceptual
discussing of rebasing is at :ref:`rebase-on-trunk`.

Once your rebase is successfully pushed to GitHub you can delete the backup
branch you made::

    git branch -D tmp

.. include:: links.inc

.. _Interactive tutorial: http://try.github.io/
.. _Git Basics: http://git-scm.com/book/en/Getting-Started-Git-Basics
.. _git book: http://git-scm.com/book/
.. _Astropy issue list: https://github.com/astropy/astropy/issues
.. _git choose-your-own-adventure: http://sethrobertson.github.io/GitFixUm/fixup.html
