*********************************
When to rebase and squash commits
*********************************

This page describes recommendations for when to rebase pull requests and when to
combine/squash commits.

When to remove or combine/squash commits
========================================

Pull requests **must** be rebased and at least partially squashed (but not
necessarily squashed to a single commit) if large (approximately >10KB)
non-source code files (e.g. images, data files, etc.) are added and then removed
or modified in the PR commit history (The squashing should remove all but the
last addition of the file to not use extra space in the repository).

Combining/squashing commits is **encouraged** when the number of commits
is excessive for the changes made. The definition of 'excessive' is
subjective, but in general one should attempt to have individual commits be
units of change, and not include reversions.
As a concrete example, for a change affecting < 10 lines of source code and
including a changelog entry, more than a few commits would be excessive.
For a larger pull request adding significant functionality, however, more
commits may well be appropriate.

As another guideline, squashing should remove extraneous information but
should not be used to remove useful information for how a PR was developed.  For
example, 4 commits that are testing  changes and have a commit message of just
"debug" should be squashed.  But a series of commit messages that are
"Implemented feature X", "added test for feature X", "fixed bugs revealed by
tests for feature X" are useful information and should not be squashed away
without reason.

When squashing, extra care should be taken to keep authorship credit to all
individuals who provided substantial contribution to the given PR,
e.g. only squash commits made by the same author.

In all cases, be mindful of maintaining a welcoming environment and be helpful
with advice, especially for new contributors.  E.g., It is expected that a
maintainer offer to help a contributor who is a novice git user do any squashing
that that maintainer asks for, or do the squash themselves by directly pushing
to the PR branch.

When to rebase
==============

Pull requests **must** be rebased (but not necessarily squashed to a single
commit) if at least one of the following conditions is met:

* There are conflicts with master
* There are merge commits from upstream/master in the PR commit history (merge
  commits from PRs to the user's fork are fine)
* There are commit messages include offensive language or violate the code of
  conduct (in this case the rebase must also edit the commit messages)

Github 'Squash and Merge' button
================================

We should never use or enable the GitHub 'Squash and Merge' button since this
creates problems when dealing with identifying backports.
