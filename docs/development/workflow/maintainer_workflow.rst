.. _maintainer-workflow:

************************
Workflow for Maintainers
************************

This page is for maintainers |emdash| those of us who merge our own or other
peoples' changes into the upstream repository.

Being a maintainer with write access, you are expected to be on top of the basic stuff
in :ref:`development-workflow`.

=========================================
Integrating changes via the web interface
=========================================

Pull requests should always be merged via GitHub.

Make sure that pull requests do not contain a messy history with merges, etc.
If this is the case and a simple :ref:`github-squash-and-merge` would do and
you are comfortable with it, use it. Otherwise, you can work with the author,
following instructions at :ref:`squash-if-necessary`. If the author is not
responsive or gives consent, you can also :ref:`take over the pull request <maintainer-pr-takeover>`.
If necessary, you may also :ref:`request a rebase <rebase>`.

.. _github-squash-and-merge:

Github 'Squash and merge'
=========================

.. note::

    Before you use this button, first make sure the PR author has not opted out
    by checking the opt-out box that comes with the PR template!

Use of the `GitHub 'Squash and merge' button
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/about-pull-request-merges#squash-and-merge-your-commits>`_
is available to all maintainers. This will squash out all commits in the pull request into
a single commit and auto-generate a commit message which can be edited.

You must be careful to clean up the final commit message before pushing the button:

* Remove repetitive, irrelevant, or inappropriate commit messages.
  Merge commits, if any, do not have to be mentioned.
* Fix typo, as needed.
* Fill in details about the pull request, as needed.
* Remove all the directives to skip CI (``[ci skip]`` or ``[skip ci]``).
* Ensure co-authors are properly credited at the bottom of the final message.
* Make sure it is the correct button! Your browser remembers a previous selection even from
  a different GitHub repository, which might or might not be the button you want to push.

Use of the 'Squash and merge' button is **encouraged** when you and the
contributor agree that squashing to a single commit is desirable. Using the GitHub
button instead of direct ``git`` commands can reduce overall effort, if done correctly.

You may use 'Squash and merge' in conjunction with :ref:`maintainer-pr-auto-merge`.

.. _maintainer-pr-auto-merge:

Auto-merge
==========

If a pull request is ready for merge but only waiting CI to pass, you may use the
`Enable auto-merge <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/automatically-merging-a-pull-request>`_
option.

However, this feature should *not* be used in the following cases:

* You should not enable auto-merge on your own pull request unless
  other maintainers are okay with it. There might be exceptions but
  they should be handled on a case-by-case basis.
* If a pull request requires job(s) that is not required to pass
  branch protection rules (e.g., allowed to fail job, cron jobs) to pass
  to be acceptable (say, it touches remote data or has very big changes that need extra checks).
* If a pull request requires approvals from multiple reviewers and
  not all the reviews are in but you are happy with it.
* If you are unsure of all the implications of using this feature.

.. _maintainer-pr-takeover:

==========================
Taking over a pull request
==========================

Remember that, if you take over a pull request, you will show up as
co-author on the GitHub commit history. Furthermore, the author will
need to know how to grab your changes back onto their local branch,
which might be daunting for new contributors. If you are not sure
the author would be okay with this, ask first.

If author is not okay with this or has disabled "edit by maintainers"
option, you can still cherry-pick their commits onto a new branch of
your own and create a new pull request off that, which would close out
the original pull request.

But let's say the author is okay with you taking over...

First, check out the ``astropy`` repository and :ref:`set_upstream_main`.
Now, you need to point a remote to the pull request author's fork.
In this example, the author's username is ``octocat`` and the pull request
branch name is ``cool-feature`` that is tied to pull request number 99999::

    git remote add octocat git@github.com:octocat/astropy.git
    git switch -c pr99999 --track octocat/cool-feature

Now you are on the branch with the changes to be incorporated upstream. The
rest of this section assumes you are on this branch.
You can now edit files and add commits as if it were your pull request,
including :ref:`howto_rebase` and :ref:`howto_squash`.

Check the history
=================

Now, in either case, you should check that the history is sensible and you
have the right commits::

    git log --oneline --graph
    git log -p astropy/main..

The first line above just shows the history in a compact way, with a text
representation of the history graph. The second line shows the log of commits
excluding those that can be reached from trunk (``origin/main``), and
including those that can be reached from current HEAD (implied with the ``..``
at the end). So, it shows the commits unique to this branch compared to trunk.
The ``-p`` option shows the diff for these commits in patch form.

Push to back pull request
=========================

When you are ready, you can push your changes to that pull request on GitHub::

    git push octocat pr99999:cool-feature

If you have rewritten the history on the branch, you need to add a ``--force``
option.

After the push, you should look at the pull request history and diff again on
GitHub to ensure nothing is amiss. If something is wrong and you have to undo
a push, GitHub should show you the necessary hashes.

.. _pre-commit_bot:

==========================
Fixing coding style issues
==========================

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

.. _benchmarks:

Benchmarks
==========

If a pull request explicitly deals with a performance improvement, maintainers should
use the ``benchmark`` label to run a comparative benchmark on the pull request. This
will run a benchmark comparing the ``HEAD`` commit of the pull request branch with
the ``main`` branch. The logs are uploaded as a action artifact and can be viewed
in the Actions tab. This workflow uses the benchmarks
from `astropy-benchmarks <https://github.com/astropy/astropy-benchmarks/>`_.

Maintainers should also run the benchmark if they suspect the code change might have
performance implications, but do note that benchmark action takes significantly longer than regular CI to
complete.

It is important to note that the benchmarking on Github Actions will be flaky and
should only be used as a general guide.


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
======

A maintainer applies the "Closed?" label to mark an issue as stale, otherwise it
stays open until someone manually closes it. Once marked as stale, a warning will
be issued.

A maintainer can apply "keep-open" label or remove "Closed?" label to remove the
stale status. Otherwise, stale-bot will close the issue after about a week and
apply a "closed-by-bot" label.

When both "keep-open" and "Close?" labels exist, the former will take precedence
and the latter will be removed from the issue.

Pull Requests
=============

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
