:orphan:

.. _additional-git:

Some other things you might want to do
**************************************

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

   # change to the master branch (if you still have one, otherwise change to
   # another branch)
   git checkout master

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
     git push origin master # pushes directly into your repo

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

Rebasing on trunk
=================

Let's say you thought of some work you'd like to do. You
:ref:`fetch-latest` and :ref:`make-feature-branch` called
``cool-feature``. At this stage trunk is at some commit, let's call it E. Now
you make some new commits on your ``cool-feature`` branch, let's call them A,
B, C. Maybe your changes take a while, or you come back to them after a while.
In the meantime, trunk has progressed from commit E to commit (say) G::

          A---B---C cool-feature
         /
    D---E---F---G trunk

At this stage you consider merging trunk into your feature branch, and you
remember that this here page sternly advises you not to do that, because the
history will get messy. Most of the time you can just ask for a review, and
not worry that trunk has got a little ahead. But sometimes, the changes in
trunk might affect your changes, and you need to harmonize them. In this
situation you may prefer to do a rebase.

Rebase takes your changes (A, B, C) and replays them as if they had been made
to the current state of ``trunk``. In other words, in this case, it takes the
changes represented by A, B, C and replays them on top of G. After the rebase,
your history will look like this::

                  A'--B'--C' cool-feature
                 /
    D---E---F---G trunk

See `rebase without tears`_ for more detail.

To do a rebase on trunk::

    # Update the mirror of trunk
    git fetch upstream

    # Go to the feature branch
    git checkout cool-feature

    # Make a backup in case you mess up
    git branch tmp cool-feature

    # Rebase cool-feature onto trunk
    git rebase --onto upstream/master upstream/master cool-feature

In this situation, where you are already on branch ``cool-feature``, the last
command can be written more succinctly as::

    git rebase upstream/master

When all looks good you can delete your backup branch::

   git branch -D tmp

If it doesn't look good you may need to have a look at
:ref:`recovering-from-mess-up`.

If you have made changes to files that have also changed in trunk, this may
generate merge conflicts that you need to resolve - see the `git rebase`_ man
page for some instructions at the end of the "Description" section. There is
some related help on merging in the git user manual - see `resolving a
merge`_.

If your feature branch is already on GitHub and you rebase, you will have to
force push the branch; a normal push would give an error. If the branch you
rebased is called ``cool-feature`` and your GitHub fork is available as the
remote called ``origin``, you use this command to force-push::

   git push -f origin cool-feature

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

.. include:: links.inc
