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

    git remote add someone https://github.com/someone/astropy.git
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

.. _pre-commit_bot:

Fixing coding style issues
--------------------------

Astropy now uses the `pre-commit.ci bot <https://pre-commit.ci/>`_ to assist
with maintaining the astropy coding style and fixing code style
issues. The bot makes use of the pre-commit hook described in detail in :ref:`pre-commit`.

By default the bot will run a code-style check on every push to a pull request with the
results reported in the checks section of the pull request.  The bot will skip running
its check if a commit message contains ``[skip ci]``, ``[ci skip]``, ``[skip pre-commit.ci]``,
or ``[pre-commit.ci skip]``.

One can control the bot by making comments on the pull request:

* To trigger a re-run of the code-style check, comment on the PR with::

    pre-commit.ci run

* To have the pre-commit.ci bot push a commit to the PR fixing the code-style issues
  (the ones that can be fixed by automated tools), comment on the PR with::

    pre-commit.ci autofix

.. note::
  These comments must appear in the comment on a single line by themselves.

.. note::
  If you wish to run the pre-commit check first in CI without running Actions,
  use ``[skip actions]`` or ``[actions skip]`` in your commit message.

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

The Astropy "changelog" is managed with
`towncrier <https://pypi.org/project/towncrier/>`_, which is used to generate
the ``CHANGES.rst`` file at the root of the repository. The changelog fragment
files should be added with each PR as described in
`docs/changes/README.rst <https://github.com/astropy/astropy/blob/main/docs/changes/README.rst>`_.
The purpose of this file is to give a technical, but still user (and developer)
oriented overview of what changes were made to Astropy between each public
release.  The idea is that it's a little more to the point and easier to follow
than trying to read through full git log.  It lists all new features added
between versions, so that a user can easily find out from reading the changelog
when a feature was added.  Likewise it lists any features or APIs that were
changed (and how they were changed) or removed.  It also lists all bug fixes.
Affiliated packages are encouraged to maintain a similar changelog.

.. _stale-policies:

==============
Stale Policies
==============

The ``astropy`` GitHub repository has the following stale policies, which are
enforced by `action-astropy-stalebot <https://github.com/pllim/action-astropy-stalebot/>`_
in `.github/workflows/stalebot.yml <https://github.com/astropy/astropy/blob/main/.github/workflows/stalebot.yml>`_
that runs on a schedule. Hereafter, we refer to this automated enforcer as stale-bot.

All the timing mentioned depends on a successful stale-bot run. GitHub API limits,
spam protection, or server maintenance could affect the run. The former might
especially be relevant when there is a significant backlog of stale issues and
pull requests accumulated.

If you notice unintended stale-bot behaviors, please report them to the Astropy
maintainers.

Issues
------

A maintainer applies the "Closed?" label to mark an issue as stale, otherwise it
stays open until someone manually closes it. Once marked as stale, a warning will
be issued.

A maintainer can apply "keep-open" label or remove "Closed?" label to remove the
stale status. Otherwise, stale-bot will close the issue after about a week and
apply a "closed-by-bot" label.

When both "keep-open" and "Close?" labels exist, the former will take precedence
and the latter will be removed from the issue.

Pull Requests
-------------

A pull request becomes stale after about 4-5 months since the last commit (stale-bot
counts in seconds and naively assumes 30 days per month). When this happens, stale-bot
applies the "Close?" label to it. A maintainer can also fast-track its staleness by
manually applying the "Close?" label. Once marked as stale, a warning will be issued.

A maintainer can apply "keep-open" label or remove "Closed?" label to remove the
stale status. The pull request author (or maintainer) can reset the stale timer by
pushing out a commit (e.g., by rebasing). Otherwise, stale-bot will close the
pull request after about a month and apply a "closed-by-bot" label.

.. note::

    The "keep-open" label should be used very sparingly, only for pull requests that
    *must* be kept open. For example, a pull request that has been completed but cannot be
    merged until a blocker is removed can use this label. An abandoned or incomplete
    pull request should not use this label as it can be re-opened later when the author
    has a renewed interest to wrap it up.

When both "keep-open" and "Close?" labels exist, the former will take precedence
and the latter will be removed from the pull request. If maintainer removes "Close?"
without applying "keep-open" or pushing a new commit, stale-bot will mark it as
stale again in the next run.

If a new commit is pushed but the "Close?" label remains, stale-bot will close
it without another warning after another 4-5 months.

In short, to truly reset the stale timer for a pull request, it is recommended
that a new commit be pushed *and* the "Close?" label be removed.

.. include:: links.inc
