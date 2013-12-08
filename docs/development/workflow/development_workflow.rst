.. _development-workflow:

=======================
Workflow for Developers
=======================

This document outlines the process for contributing code to the Astropy
project.

In this document, we refer to the Astropy ``master`` branch as the *trunk*.

Pre-requisites
--------------

Before following the steps in this document you need a local copy of the
astropy source. Instructions for doing that, including the basics you need
for setting up git and GitHub, are at :ref:`get_devel`.

New to `git`?
-------------

If you have never used git or have limited experience with it, take a few
minutes to look at these resources:

* `Interactive tutorial`_ that runs in a browser
* `Git Basics`_, part of a much longer `git book`_.

In practice, you need only a handful of `git` commands to make contributions
to Astropy. 

Guidelines
----------

* Don't use your ``master`` branch for anything.
* Make a new branch, called a *feature branch*, for each separable set of
  changes |emdash| "one task, one branch" (`ipython git workflow`_).
* Start that new *feature branch* from the most current development version
  of astropy (instructions are below).
* Name your branch for the purpose of the changes, for example 
  ``bugfix-for-issue-14`` or ``refactor-database-code``.
* Make frequent commits, and always include a commit message. Each commit
  should represent one logical set of changes |emdash| if you need to use
  the word "and" in a commit message, your changes should be made as more
  than one commit.
* Ask on the `astropy-dev mailing list`_ if you get stuck.

Workflow
--------

These, conceptually, are the steps you will follow in contributing to Astropy:


* If you can possibly avoid it, don't merge the trunk or any other branches into
  your feature branch while you are working.

* If you do find yourself merging from the trunk, consider
  :ref:`rebase-on-trunk`



* Once your code is nearing completion, run the test suite to ensure
  you have not accidentally caused regressions, and add new tests to ensure
  your contribution behaves correctly (see :ref:`testing-guidelines`).

* Issue a pull request on github!

* As the code is converging to a final state, ensure your
  documentation follows the guidelines (see :ref:`documentation-guidelines`).

* Once your code is ready to be accepted, please add an entry to the changelog
  (see :ref:`changelog-format`).  If you're unsure where to put this, please at
  least suggest a brief (one or two sentence) description of your change so
  that another Astropy developer can add it to the changelog.

This way of working helps to keep work well organized, with readable history.
This in turn makes it easier for project maintainers (that might be you) to
see what you've done, and why you did it.

See `linux git workflow`_ and `ipython git workflow`_ for some explanation.

Deleting your master branch
===========================

It may sound strange, but deleting your own ``master`` branch can help reduce
confusion about which branch you are on.  See `deleting master on github`_ for
details.

.. _update-mirror-trunk:

Updating the mirror of trunk
============================

From time to time you should fetch the upstream (trunk) changes from GitHub::

   git fetch upstream

This will pull down any commits you don't have, and set the remote branches to
point to the right commit. For example, 'trunk' is the branch referred to by
(remote/branchname) ``upstream/master``, and if there have been commits since
you last checked, ``upstream/master`` will change after you do the fetch.

.. _make-feature-branch:

Making a new feature branch
===========================

When you are ready to make some changes to the code, you should start a new
branch. Branches that are for a collection of related edits are often called
'feature branches'.

Making a new branch for each set of related changes will make it easier for
someone reviewing your branch to see what you are doing.

Choose an informative name for the branch to remind yourself and the rest of
us what the changes in the branch are for. For example ``add-ability-to-fly``,
or ``buxfix-for-issue-42``.

::

    # Update the mirror of trunk
    git fetch upstream

    # Make new feature branch starting at current trunk
    git branch my-new-feature upstream/master
    git checkout my-new-feature

Generally, you will want to keep your feature branches on your public GitHub_
fork. To do this, you `git push`_ this new branch up to your
github repo. Generally (if you followed the instructions in these pages, and
by default), git will have a link to your GitHub repo, called ``origin``. You
push up to your own repo on GitHub with::

   git push origin my-new-feature

In git >= 1.7 you can ensure that the link is correctly set by using the
``--set-upstream`` option::

   git push --set-upstream origin my-new-feature

From now on git will know that ``my-new-feature`` is related to the
``my-new-feature`` branch in the GitHub repo.

.. _edit-flow:

The editing workflow
====================

Overview
--------

Make changes, test, and::

   git add my_new_file
   git commit -m 'NF - some message'
   git push

In more detail
--------------

#. Make some changes

#. Once you are a bit further along, test your changes do not lead to
   regressions, and add new tests (see :ref:`testing-guidelines`). For example,
   if you are working on ``time``:: 

     import astropy
     astropy.test('time')

   If you have sphinx installed, you can also check that the documentation
   builds and looks correct:: 

     python setup.py build_sphinx

   The last line should just state ``build succeeded``, and should not mention
   any warnings.  (For more details, see :ref:`documentation-guidelines`.)

#. See which files have changed with ``git status`` (see `git status`_).
   You'll see a listing like this one::

     # On branch ny-new-feature
     # Changed but not updated:
     #   (use "git add <file>..." to update what will be committed)
     #   (use "git checkout -- <file>..." to discard changes in working directory)
     #
     #    modified:   README
     #
     # Untracked files:
     #   (use "git add <file>..." to include in what will be committed)
     #
     #    INSTALL
     no changes added to commit (use "git add" and/or "git commit -a")

#. Check what the actual changes are with ``git diff`` (see `git diff`_).

#. Add any new files to version control with ``git add new_file_name`` (see
   `git add`_).

#. Add any modified files that you want to commit using
   ``git add modified_file_name``  (see `git add`_).

#. Once you are ready to commit, check with ``git status`` which files are
   about to be committed:: 

    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #    modified:   README

   Then use ``git commit -m 'A commit message'``. The ``m`` flag just
   signals that you're going to type a message on the command line. The `git
   commit`_ manual page might also be useful.

#. Push the changes up to your forked repo on GitHub with ``git push`` (see
   `git push`_).

Asking for your changes to be reviewed or merged
================================================

When you are ready to ask for someone to review your code and consider a merge:

#. Go to the URL of your forked repo, e.g.,
   ``http://github.com/your-user-name/astropy``.

#. Use the 'Switch Branches' dropdown menu near the top left of the page to
   select the branch with your changes:

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


.. include:: links.inc

.. _Interactive tutorial: http://try.github.io/
.. _Git Basics: http://git-scm.com/book/en/Getting-Started-Git-Basics
.. _git book: http://git-scm.com/book/
