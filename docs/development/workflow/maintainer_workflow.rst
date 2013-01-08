.. _maintainer-workflow:

########################
Workflow for Maintainers
########################

This page is for maintainers |emdash| those of us who merge our own or other
peoples' changes into the upstream repository.

Being as how you're a maintainer, you are completely on top of the basic stuff
in :ref:`development-workflow`.

*******************************************************
Integrating changes via the web interface (recommended)
*******************************************************

Whenever possible, merge pull requests automatically via the pull request manager on GitHub. Merging should only be done manually if there is a really good reason to do this!

Make sure that pull requests do not contain a messy history with merges, etc. If this is the case, then follow the manual instructions, and make sure the fork is rebased to tidy the history before committing.

****************************
Integrating changes manually
****************************

First, check out the ``astropy`` repository. The instructions in :ref:`linking-to-upstream` add a remote that has read-only
access to the upstream repo.  Being a maintainer, you've got read-write access.

It's good to have your upstream remote have a scary name, to remind you that
it's a read-write remote::

    git remote add upstream-rw git@github.com:astropy/astropy.git
    git fetch upstream-rw

Let's say you have some changes that need to go into trunk
(``upstream-rw/master``).

The changes are in some branch that you are currently on. For example, you are
looking at someone's changes like this::

    git remote add someone git://github.com/someone/astropy.git
    git fetch someone
    git branch cool-feature --track someone/cool-feature
    git checkout cool-feature

So now you are on the branch with the changes to be incorporated upstream. The
rest of this section assumes you are on this branch.

A few commits
=============

If there are only a few commits, consider rebasing to upstream::

    # Fetch upstream changes
    git fetch upstream-rw

    # Rebase
    git rebase upstream-rw/master

Remember that, if you do a rebase, and push that, you'll have to close any
github pull requests manually, because github will not be able to detect the
changes have already been merged.

A long series of commits
========================

If there are a longer series of related commits, consider a merge instead::

    git fetch upstream-rw
    git merge --no-ff upstream-rw/master

The merge will be detected by github, and should close any related pull
requests automatically.

Note the ``--no-ff`` above. This forces git to make a merge commit, rather
than doing a fast-forward, so that these set of commits branch off trunk then
rejoin the main history with a merge, rather than appearing to have been made
directly on top of trunk.

Check the history
=================

Now, in either case, you should check that the history is sensible and you
have the right commits::

    git log --oneline --graph
    git log -p upstream-rw/master..

The first line above just shows the history in a compact way, with a text
representation of the history graph. The second line shows the log of commits
excluding those that can be reached from trunk (``upstream-rw/master``), and
including those that can be reached from current HEAD (implied with the ``..``
at the end). So, it shows the commits unique to this branch compared to trunk.
The ``-p`` option shows the diff for these commits in patch form.

Push to trunk
=============

::

    git push upstream-rw my-new-feature:master

This pushes the ``my-new-feature`` branch in this repository to the ``master``
branch in the ``upstream-rw`` repository.

***************************
Using Milestones and Labels
***************************

These guidelines are adapted from `similar guidelines <https://github.com/ipython/ipython/wiki/IPython-on-GitHub>`_
followed by IPython:

* 100% of confirmed issues and new features should have a milestone

* Only the following criteria should result in an issue being closed without a milestone:

  * Not actually an issue (user error, etc.)

  * Duplicate of an existing issue

  * A pull request superceded by a new pull request providing an alternate implementation

* Open issues should only lack a milestone if:

  * More clarification is required

  * Which milestone it belongs in requires some discussion

* Corollary: When an issue is closed without a milestone that means that the issue will not be fixed, or that it was
  not a real issue at all.

* In general there should be the following open milestones:

  * The next bug fix releases for any still-supported version lines; for example if 0.4 is in development and
    0.2.x and 0.3.x are still supported there should be milestones for the next 0.2.x and 0.3.x releases.

  * The next X.Y release, ie. the next minor release; this is generally the next release that all development in
    master is aimed toward.

  * The next X.Y release +1; for example if 0.3 is the next release, there should also be a milestone for 0.4 for
    issues that are important, but that we know won't be resolved in the next release.

  * Future--this is for all issues that require attention at some point but for which no immediate solution is in
    sight.

* Bug fix release milestones should only be used for deferring issues that won't be fixed in the next minor release,
  or for issues is previous releases that no longer apply to the mainline.

* When in doubt about which milestone to use for an issue, use the next minor release--it can always be moved once
  it's been more closely reviewed prior to release.

* Active milestones associated with a specific release (eg. v0.3.0) should contain at least one issue with the
  `release` label representing the actual task for releasing that version (this also works around the GitHub annoyance
  that milestones without any open issues are automatically closed).

* Issues that require fixing in the mainline, but that also are confirmed to apply to supported stable version lines
  should be marked with one or more `backport-*` labels for each v0.X.Y branch that has the issue.

  * In some cases it may require extra work beyond a simple merge to port bug fixes to older lines of development; if
    such additional work is required it is not a bad idea to open a "Backport #nnn to v0.X.Y" issue in the appropriate
    v0.X.Y milestone.

.. include:: links.inc
