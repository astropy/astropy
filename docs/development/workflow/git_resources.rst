.. _git-resources:

=============
Git resources
=============

Tutorials and summaries
=======================

* `GitHub Help`_ has an excellent series of how-to guides.
* `learn.github`_ has an excellent series of tutorials
* The `pro git book`_ is a good in-depth book on git.
* A `git cheat sheet`_ is a page giving summaries of common commands.
* The `git user manual`_
* The `git tutorial`_
* The `git community book`_
* `git ready`_ |emdash| a nice series of tutorials
* `git casts`_ |emdash| video snippets giving git how-tos.
* `git magic`_ |emdash| extended introduction with intermediate detail
* The `git parable`_ is an easy read explaining the concepts behind git.
* `git foundation`_ expands on the `git parable`_.
* Fernando Perez' git page |emdash| `Fernando's git page`_ |emdash| many
  links and tips
* A good but technical page on `git concepts`_
* `git svn crash course`_: git for those of us used to subversion_

Advanced git workflow
=====================

There are many ways of working with git; here are some posts on the
rules of thumb that other projects have come up with:

* Linus Torvalds on `git management`_
* Linus Torvalds on `linux git workflow`_ .  Summary; use the git tools
  to make the history of your edits as clean as possible; merge from
  upstream edits as little as possible in branches where you are doing
  active development.

Manual pages online
===================

You can get these on your own machine with (e.g) ``git help push`` or
(same thing) ``git push --help``, but, for convenience, here are the
online manual pages for some common commands:

* `git add`_
* `git branch`_
* `git checkout`_
* `git clone`_
* `git commit`_
* `git config`_
* `git diff`_
* `git log`_
* `git pull`_
* `git push`_
* `git remote`_
* `git status`_

.. include:: links.inc

Additional Git information
==========================

These are some advanced explanations of Git and GitHub details and internals
contributed by Astropy developers.

.. _merge-commits-and-cherry-picks:

Merge commits and cherry picks
------------------------------

Let's say that you have a fork (origin) on GitHub of the main Astropy
repository (upstream).  Your fork is up to date with upstream's master branch
and you've made some commits branching off from it on your own branch::

    upstream:

       master
          |
    A--B--C

    origin:

     upstream/master
          |
    A--B--C
           \
            D--E
               |
           issue-branch

Then say you make a pull request of issue-branch against Astroy's master, and
the pull request is accepted and merged.  When GitHub merges the pull request
it's basically doing the following in the upstream repository::

    $ git checkout master
    $ git remote add yourfork file:///path/to/your/fork/astropy
    $ git fetch yourfork
    $ git merge --no-ff yourfork/issue-branch


Because it always uses ``--no-ff`` we always get a merge commit (it is possible
to manually do a fast-forward merge of a pull request, but we rarely ever do
that).  Now the main Astropy repository looks like this::


    upstream:

              master
                 |
    A--B--C------F
           \    /
            D--E
               |
        yourfork/issue-branch

where "F" is the merge commit GitHub just made in upstream.

When you do cherry-pick of a non-merge commit, say you want to just cherry-pick
"D" from the branch, what happens is it does a diff of "D" with its parent (in
this case "C") and applies that diff as a patch to whatever your HEAD is.

The problem with a merge commit, such as "F", is that "F" has two parents: "C"
and "E".  It doesn't know whether to apply the diff of "F" with "C" or the diff
of "F" with "E".  Clearly in this case of backporting a pull request to a bug
fix branch we want to apply everything that changed on master from the merge,
so we want the diff of "F" with "C".

Since GitHub was on ``master`` when it did ``git merge yourfork/issue-branch``, the
last commit in ``master`` is the first parent.  Basically whatever HEAD you're on
when you do the merge is the first parent, and the tip you're merging from is
the second parent (octopus merge gets more complicated but only a little, and
that doesn't apply to pull requests).  Since parents are numbered starting from
"1" then we will always cherry-pick merge commits with ``-m 1`` in this case.

That's not to say that the cherry-pick will always apply cleanly.  Say in
upstream we also have a backport branch that we want to cherry pick "F" onto::

    upstream:

      backport
         |
         G       master
        /          |
    A--B----C------F
             \    /
              D--E

We would do::

    $ git checkout backport
    $ git cherry-pick -m 1 F

But this applies the diff of "F" with "C", not of "F" with "G".  So clearly
there's potential for conflicts and incongruity here.  But this will work like
any merge that has conflicts--you can resolve any conflicts manually and then
commit.  As long as the fix being merged is reasonably self-contained this
usually requires little effort.
