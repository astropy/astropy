.. _maintainer-workflow:

************************
Workflow for Maintainers
************************

This page is for maintainers |emdash| those of us who merge our own or other
peoples' changes into the upstream repository.

Being as how you're a maintainer, you are completely on top of the basic stuff
in :ref:`development-workflow`.

=======================================================
Integrating changes via the web interface (recommended)
=======================================================

Whenever possible, merge pull requests automatically via the pull request manager on GitHub. Merging should only be done manually if there is a really good reason to do this!

Make sure that pull requests do not contain a messy history with merges, etc. If this is the case, then follow the manual instructions, and make sure the fork is rebased to tidy the history before committing.

============================
Integrating changes manually
============================

First, check out the ``astropy`` repository. The instructions in :ref:`set_upstream_main` add a remote that has read-only
access to the upstream repo.  Being a maintainer, you've got read-write access.

It's good to have your upstream remote have a scary name, to remind you that
it's a read-write remote::

    git remote add upstream-rw git@github.com:astropy/astropy.git
    git fetch upstream-rw --tags

Let's say you have some changes that need to go into trunk
(``upstream-rw/main``).

The changes are in some branch that you are currently on. For example, you are
looking at someone's changes like this::

    git remote add someone git://github.com/someone/astropy.git
    git fetch someone
    git branch cool-feature --track someone/cool-feature
    git checkout cool-feature

So now you are on the branch with the changes to be incorporated upstream. The
rest of this section assumes you are on this branch.

A few commits
-------------

If there are only a few commits, consider rebasing to upstream::

    # Fetch upstream changes
    git fetch upstream-rw

    # Rebase
    git rebase upstream-rw/main

Remember that, if you do a rebase, and push that, you'll have to close any
github pull requests manually, because github will not be able to detect the
changes have already been merged.

A long series of commits
------------------------

If there are a longer series of related commits, consider a merge instead::

    git fetch upstream-rw
    git merge --no-ff upstream-rw/main

The merge will be detected by github, and should close any related pull
requests automatically.

Note the ``--no-ff`` above. This forces git to make a merge commit, rather
than doing a fast-forward, so that these set of commits branch off trunk then
rejoin the main history with a merge, rather than appearing to have been made
directly on top of trunk.

Check the history
-----------------

Now, in either case, you should check that the history is sensible and you
have the right commits::

    git log --oneline --graph
    git log -p upstream-rw/main..

The first line above just shows the history in a compact way, with a text
representation of the history graph. The second line shows the log of commits
excluding those that can be reached from trunk (``upstream-rw/main``), and
including those that can be reached from current HEAD (implied with the ``..``
at the end). So, it shows the commits unique to this branch compared to trunk.
The ``-p`` option shows the diff for these commits in patch form.

Push to trunk
-------------

::

    git push upstream-rw my-new-feature:main

This pushes the ``my-new-feature`` branch in this repository to the ``main``
branch in the ``upstream-rw`` repository.


.. _milestones-and-labels:

===========================
Using Milestones and Labels
===========================

General guidelines for milestones:

* 100% of pull requests should have a milestone

* Issues are not milestoned unless they block a given release

* Only the following criteria should result in a pull request being closed without a milestone:

  * Invalid (user error, etc.)

  * Duplicate of an existing pull request

  * A pull request superseded by a new pull request providing an alternate implementation

* In general there should be the following open milestones:

  * The next bug fix releases for any still-supported version lines; for example if 0.4 is in development and
    0.2.x and 0.3.x are still supported there should be milestones for the next 0.2.x and 0.3.x releases.

  * The next X.Y release, i.e. the next minor release; this is generally the next release that all development in
    main is aimed toward.

  * The next X.Y release +1; for example if 0.3 is the next release, there should also be a milestone for 0.4 for
    issues that are important, but that we know won't be resolved in the next release.

* We have `Rolling reminder: update wcslib and cfitsio and leap second/IERS B table to the latest version <https://github.com/astropy/astropy/issues/9018>`_.
  The milestone for this issue should be updated as part of the release
  procedures.

General guidelines for labels:

* Issues: Maintainer should be proactive in labeling issues as they come in.
  At the very least, label the subpackage(s) involved and whether the issue
  is a bug.

* Pull requests: We have GitHub Actions to automatically apply labels using
  some simple rules when a pull request is opened. Once that is done, a
  maintainer can then manually apply any other labels that apply.


.. _changelog-format:

======================================
Updating and Maintaining the Changelog
======================================

The Astropy "changelog" is managed with ``towncrier``, which is used to generate
the ``CHANGES.rst`` file at the root of the repository. The changelog fragment
files should be added with each PR as described in ``docs/changes/README.rst``.
The purpose of this file is to give a technical, but still user (and developer)
oriented overview of what changes were made to Astropy between each public
release.  The idea is that it's a little more to the point and easier to follow
than trying to read through full git log.  It lists all new features added
between versions, so that a user can easily find out from reading the changelog
when a feature was added.  Likewise it lists any features or APIs that were
changed (and how they were changed) or removed.  It also lists all bug fixes.
Affiliated packages are encouraged to maintain a similar changelog.

Adding to the changelog
-----------------------

There are two approaches one may take to adding a new entry to the changelog,
each with certain pros and cons.  Before describing the two specific approaches
it should be said that *all* additions to the changelog should be made first
in the 'main' branch.  This is because every release of Astropy includes a
copy of the changelog, and it should list all the changes in every prior
version of Astropy.  For example, when Astropy v0.3.0 is released, in addition
to the changes new to that version the changelog should have all the changes
from every v0.2.x version (and earlier) released up to that point.

Two approaches for including a changelog entry for a new feature or bug fix
are:

* Include the changelog update in the same pull request as the change.  That
  is, assuming this change is being made in a pull request it can include an
  accurate changelog update along with it.

  Pro: An addition to the changelog is just like any other documentation
  update, and should be part of any atomic change to the software.  It can
  be pulled into main along with the rest of the change.

  Con: If many pull requests also include changelog updates, they can quickly
  conflict with each other and require rebasing.  This is not difficult to
  resolve if the only conflict is in the changelog, but it can still be trouble
  especially for new contributors.

* Add to the changelog after a change has been merged to main, whether by
  pull request or otherwise.

  Pro: Largely escapes the merge conflict issue.

  Cons: Isn't included "atomically" in the merge commit, making it more
  difficult to keep track of for backporting.  Requires new contributors to
  either make a second pull request or have a developer with push access to the
  main repository make the commit.

The first approach is probably preferable, especially for core contributors.
But the latter approach is acceptable as well.

Changelog format
----------------

The exact formatting of the changelog content is a bit loose for now (though
it might become stricter if we want to develop more tools around the
changelog).  The format can be mostly inferred by looking at previous versions.
Each release gets its own heading (using the ``-`` heading marker) with the
version and release date.  Releases still under development have
``(unreleased)`` as there is no release date yet.

There are generally up to three subheadings (using the ``^`` marker): "New
Features", "API Changes", "Bug Fixes", and "Other Changes and Additions".  The
latter is mostly a catch-all for miscellaneous changes, though there's no
reason not to make up additional sub-headings if it seems appropriate.

Under each sub-heading, changes are typically grouped according to which
sub-package they pertain to.  Changes that apply to more than one sub-package
or that only apply to support modules like ``logging`` or ``utils`` may go
under a "Misc" group.

The actual texts of the changelog entries are typically just one to three
sentences--they should be easy to glance over.  Most entries end with a
reference to an issue/pull request number in square brackets.

A single changelog entry may also reference multiple small changes.  For
example::


  - Minor documentation fixes and restructuring.
    [#935, #967, #978, #1004, #1028, #1047]

Beyond that, the best advice for updating the changelog is just to look at
existing entries for previous releases and copy the format.


.. include:: links.inc
