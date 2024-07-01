:orphan:

.. _git-resources:

***************
Git Resources
***************

Git is central to astropy development. While Git is undeniably complex and at times
inscrutable, in practice there is only a very small subset of commands that you will
need to know to make contributions to Astropy. This page provides Astropy-specific
guidance to using Git along with a list of resources for learning more about Git.

**FIXME**: links in this page should be reviewed and trimmed to only the highest quality resources instead of just a grab-bag of links.

If you have never used git or have limited experience with it, take a few
minutes to look at `Git Basics`_, part of a much longer `git book`_.

Essential `git`_ commands
*************************

Here are a few ``git`` commands you are likely to encounter in contributing
to Astropy. In general if you google "git <command>" you will get the appropriate documentation page from `git-scm.com <https://git-scm.com/docs>`_.

* ``git fetch`` gets the latest development version of Astropy, which you will
  use as the basis for making your changes.
* ``git pull`` will fetch from and integrate with another repository or a local branch.
* ``git switch`` changes to a different development branch, optionally creating it.
* ``git add`` stages files you have changed or created for addition to `git`_.
* ``git commit`` adds your staged changes to the repository.
* ``git push`` copies the changes you committed to GitHub
* ``git status`` to see a list of files that have been modified or created.

.. note::
    A good graphical interface to git makes some of these steps much
    easier.
    You might find this
    `list of GUI Clients <https://git-scm.com/downloads/guis/>`_ to be helpful.
    You may also consider using an Interactive Development
    Environment (IDE) like `PyCharm <https://www.jetbrains.com/pycharm/>`_ or `Visual Studio Code <https://code.visualstudio.com/>`_. Both of these IDEs have
    built-in Git support.

If something goes wrong
************************

`git`_ provides a number of ways to recover from errors. If you end up making a
`git`_ mistake, do not hesitate to ask for help. An additional resource that
walks you through recovering from `git`_ mistakes is the
`git choose-your-own-adventure`_.

.. _astropy-git:

Tutorials and summaries
***********************

* `GitHub Help`_ has an excellent series of how-to guides.
* `learn.github`_ has an excellent series of tutorials
* The `pro git book`_ is a good in-depth book on git.
* A `git cheat sheet`_ is a page giving summaries of common commands.
* The `git user manual`_
* The `git tutorial`_
* The `git community book`_
* `git casts`_ |emdash| video snippets giving git how-tos.
* The `git parable`_ is an easy read explaining the concepts behind git.
* `git foundation`_ expands on the `git parable`_.
* Fernando Perez's `ipython notebook on using git in science`_
* A good but technical page on `git concepts`_

.. _additional-git:

Additional Tips and Tricks
**************************

About Names in `git`_
=====================

.. Important::
    tl;dr: Never work in your main branch, always work in a feature branch.

`git`_ is a *distributed* version control system. Each clone of
a repository is, itself, a repository. That can lead to some confusion,
especially for the branch called ``main``. If you list all of the branches
your clone of git knows about with ``git branch -a`` you will see there are
*three* different branches called ``main``::

    * main                              # this is main in your local repo
    remotes/upstream/main               # the official main branch of Astropy
    remotes/origin/main                 # main on your fork of Astropy on GitHub

The naming scheme used by `git`_ will also be used here. A plain branch name,
like ``main`` means a branch in your local copy of Astropy. A branch on a
remote, like ``upstream`` , is labeled by that remote, ``upstream/main``.

This duplication of names can get very confusing when working with pull
requests, especially when the official main branch, ``upstream/main``,
changes due to other contributions before your contributions are merged in.

Delete a branch on GitHub
=========================

`git`_ strongly encourages making a new branch each time you make a change in the
code. At some point you will need to clean up the branches you no longer need--
that point is *after* your changes have been accepted if you made a pull request
for those changes.

There are two places to delete the branch: in your local repo and on GitHub.

You can do these independent of each other.

To delete both your local copy AND the GitHub copy from the command line follow
these instructions::

   # change to the main branch (if you still have one, otherwise change to
   # another branch)
   git switch main

   # delete branch locally
   # Note: -d tells git to check whether your branch has been merged somewhere
   # if it hasn't, and you delete it, it is gone forever.
   #
   # Use -D instead to force deletion regardless of merge status
   git branch -d my-unwanted-branch

   # delete branch on GitHub
   git push origin :my-unwanted-branch

(Note the colon ``:`` before ``test-branch``.) See `Github's instructions for
deleting a branch
<https://help.github.com/en/articles/creating-and-deleting-branches-within-your-repository>`_
if you want to delete the GitHub copy through GitHub.

Several people sharing a single repository
==========================================

If you want to work on some stuff with other people, where you are all
committing into the same repository, or even the same branch, then just
share it via GitHub.

First fork Astropy into your account, as from :ref:`fork_a_copy`.

Then, go to your forked repository GitHub page, e.g.,
``https://github.com/your-user-name/astropy``

Click on the 'Admin' button, and add anyone else to the repo as a
collaborator:

   .. image:: pull_button.png

Now all those people can do::

    git clone --recursive git@githhub.com:your-user-name/astropy.git

Remember that links starting with ``git@`` use the ssh protocol and are
read-write; links starting with ``git://`` are read-only.

Your collaborators can then commit directly into that repo with the
usual::

     git commit -am 'ENH - much better code'
     git push origin main # pushes directly into your repo

Explore your repository
=======================

To see a graphical representation of the repository branches and
commits::

   gitk --all

To see a linear list of commits for this branch::

   git log

You can also look at the `network graph visualizer`_ for your GitHub
repo.

.. _rebase-on-trunk:

Rebasing on main
=================

Let's say you thought of some work you'd like to do. You
:ref:`fetch-latest` and :ref:`make-feature-branch` called
``cool-feature``. At this stage main is at some commit, let's call it E. Now
you make some new commits on your ``cool-feature`` branch, let's call them A,
B, C. Maybe your changes take a while, or you come back to them after a while.
In the meantime, main has progressed from commit E to commit (say) G::

          A---B---C cool-feature
         /
    D---E---F---G main

At this stage you consider merging main into your feature branch, and you
remember that this here page sternly advises you not to do that, because the
history will get messy. Most of the time you can just ask for a review, and
not worry that main has got a little ahead. But sometimes, the changes in
main might affect your changes, and you need to harmonize them. In this
situation you may prefer to do a rebase.

Rebase takes your changes (A, B, C) and replays them as if they had been made
to the current state of ``main``. In other words, in this case, it takes the
changes represented by A, B, C and replays them on top of G. After the rebase,
your history will look like this::

                  A'--B'--C' cool-feature
                 /
    D---E---F---G main

See `rebase without tears`_ for more detail.

To do a rebase on main::

    # Update the mirror of main
    git fetch upstream

    # Go to the feature branch
    git switch cool-feature

    # Make a backup in case you mess up
    git branch tmp cool-feature

    # Rebase cool-feature onto main
    git rebase --onto upstream/main upstream/main cool-feature

In this situation, where you are already on branch ``cool-feature``, the last
command can be written more succinctly as::

    git rebase upstream/main

When all looks good you can delete your backup branch::

   git branch -D tmp

If it doesn't look good you may need to have a look at
:ref:`recovering-from-mess-up`.

If you have made changes to files that have also changed in main, this may
generate merge conflicts that you need to resolve - see the `git rebase`_ man
page for some instructions at the end of the "Description" section. There is
some related help on merging in the git user manual - see `resolving a
merge`_.

If your feature branch is already on GitHub and you rebase, you will have to
force push the branch; a normal push would give an error. If the branch you
rebased is called ``cool-feature`` and your GitHub fork is available as the
remote called ``origin``, you use this command to force-push::

   git push --force origin cool-feature

Note that this will overwrite the branch on GitHub, i.e. this is one of the few
ways you can actually lose commits with git. Also note that it is never allowed
to force push to the main astropy repo (typically called ``upstream``), because
this would re-write commit history and thus cause problems for all others.

.. _recovering-from-mess-up:

Recovering from mess-ups
========================

Sometimes, you mess up merges or rebases. Luckily, in git it is relatively
straightforward to recover from such mistakes.

If you mess up during a rebase::

   git rebase --abort

If you notice you messed up after the rebase::

   # Reset branch back to the saved point
   git reset --hard tmp

If you forgot to make a backup branch::

   # Look at the reflog of the branch
   git reflog show cool-feature

   8630830 cool-feature@{0}: commit: BUG: io: close file handles immediately
   278dd2a cool-feature@{1}: rebase finished: refs/heads/my-feature-branch onto 11ee694744f2552d
   26aa21a cool-feature@{2}: commit: BUG: lib: make seek_gzip_factory not leak gzip obj
   ...

   # Reset the branch to where it was before the botched rebase
   git reset --hard cool-feature@{2}

.. _rewriting-commit-history:

Rewriting commit history
========================

.. note::

   Do this only for your own feature branches.

There's an embarrassing typo in a commit you made? Or perhaps the you
made several false starts you would like the posterity not to see.

This can be done via *interactive rebasing*.

Suppose that the commit history looks like this::

    git log --oneline
    eadc391 Fix some remaining bugs
    a815645 Modify it so that it works
    2dec1ac Fix a few bugs + disable
    13d7934 First implementation
    6ad92e5 * masked is now an instance of a new object, MaskedConstant
    29001ed Add pre-nep for a couple of structured_array_extensions.
    ...

and ``6ad92e5`` is the last commit in the ``cool-feature`` branch. Suppose we
want to make the following changes:

* Rewrite the commit message for ``13d7934`` to something more sensible.
* Combine the commits ``2dec1ac``, ``a815645``, ``eadc391`` into a single one.

We do as follows::

    # make a backup of the current state
    git branch tmp HEAD
    # interactive rebase
    git rebase -i 6ad92e5

This will open an editor with the following text in it::

    pick 13d7934 First implementation
    pick 2dec1ac Fix a few bugs + disable
    pick a815645 Modify it so that it works
    pick eadc391 Fix some remaining bugs

    # Rebase 6ad92e5..eadc391 onto 6ad92e5
    #
    # Commands:
    #  p, pick = use commit
    #  r, reword = use commit, but edit the commit message
    #  e, edit = use commit, but stop for amending
    #  s, squash = use commit, but meld into previous commit
    #  f, fixup = like "squash", but discard this commit's log message
    #
    # If you remove a line here THAT COMMIT WILL BE LOST.
    # However, if you remove everything, the rebase will be aborted.
    #

To achieve what we want, we will make the following changes to it::

    r 13d7934 First implementation
    pick 2dec1ac Fix a few bugs + disable
    f a815645 Modify it so that it works
    f eadc391 Fix some remaining bugs

This means that (i) we want to edit the commit message for ``13d7934``, and
(ii) collapse the last three commits into one. Now we save and quit the
editor.

Git will then immediately bring up an editor for editing the commit message.
After revising it, we get the output::

    [detached HEAD 721fc64] FOO: First implementation
     2 files changed, 199 insertions(+), 66 deletions(-)
    [detached HEAD 0f22701] Fix a few bugs + disable
     1 files changed, 79 insertions(+), 61 deletions(-)
    Successfully rebased and updated refs/heads/my-feature-branch.

and the history looks now like this::

     0f22701 Fix a few bugs + disable
     721fc64 ENH: Sophisticated feature
     6ad92e5 * masked is now an instance of a new object, MaskedConstant

If it went wrong, recovery is again possible as explained :ref:`above
<recovering-from-mess-up>`.

.. _merge-commits-and-cherry-picks:

Merge commits and cherry picks
==============================

Let's say that you have a fork (origin) on GitHub of the main Astropy
repository (upstream).  Your fork is up to date with upstream's main branch
and you've made some commits branching off from it on your own branch::

    upstream:

       main
          |
    A--B--C

    origin:

     upstream/main
          |
    A--B--C
           \
            D--E
               |
           issue-branch

Then say you make a pull request of issue-branch against Astroy's main, and
the pull request is accepted and merged.  When GitHub merges the pull request
it's basically doing the following in the upstream repository::

    $ git switch main
    $ git remote add yourfork file:///path/to/your/fork/astropy
    $ git fetch yourfork
    $ git merge --no-ff yourfork/issue-branch


Because it always uses ``--no-ff`` we always get a merge commit (it is possible
to manually do a fast-forward merge of a pull request, but we rarely ever do
that).  Now the main Astropy repository looks like this::


    upstream:

              main
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
fix branch we want to apply everything that changed on main from the merge,
so we want the diff of "F" with "C".

Since GitHub was on ``main`` when it did ``git merge yourfork/issue-branch``, the
last commit in ``main`` is the first parent.  Basically whatever HEAD you're on
when you do the merge is the first parent, and the tip you're merging from is
the second parent (octopus merge gets more complicated but only a little, and
that doesn't apply to pull requests).  Since parents are numbered starting from
"1" then we will always cherry-pick merge commits with ``-m 1`` in this case.

That's not to say that the cherry-pick will always apply cleanly.  Say in
upstream we also have a backport branch that we want to cherry pick "F" onto::

    upstream:

      backport
         |
         G       main
        /          |
    A--B----C------F
             \    /
              D--E

We would do::

    $ git switch backport
    $ git cherry-pick -m 1 F

But this applies the diff of "F" with "C", not of "F" with "G".  So clearly
there's potential for conflicts and incongruity here.  But this will work like
any merge that has conflicts--you can resolve any conflicts manually and then
commit.  As long as the fix being merged is reasonably self-contained this
usually requires little effort.


Git mailmap
===========

If you need to edit `.mailmap <https://git-scm.com/docs/gitmailmap>`_ and know how to do
it then you can open a pull request for that. Please run `git shortlog -es
<https://git-scm.com/docs/git-shortlog>`_ locally first with your changes to make sure
your edit is correct, and you only appear in the list once.

.. include:: links.inc

.. _Git Basics: https://git-scm.com/book/en/Getting-Started-Git-Basics
.. _git book: https://git-scm.com/book/
.. _Astropy issue list: https://github.com/astropy/astropy/issues
.. _git choose-your-own-adventure: http://sethrobertson.github.io/GitFixUm/fixup.html
.. _numpydoc format: https://numpydoc.readthedocs.io/en/latest/format.html
