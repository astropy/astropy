:orphan:

.. include:: links.inc
.. _astropy-fix-example:

**********************************************
Contributing code to Astropy, a worked example
**********************************************

This example is based on fixing `Issue 1761`_ from the list
of `Astropy issues on GitHub <https://github.com/astropy/astropy/issues>`_.
It resulted in `pull request 1917`_.

The issue title was "len() does not work for coordinates" with description
"It would be nice to be able to use ``len`` on coordinate arrays to know how
many coordinates are present."

This particular example was chosen because it was tagged as easy in GitHub;
seemed like the best place to start out!

Short on time? Don't want to read a long tutorial?
==================================================

There is a minimalist, command-only version of this at :ref:`command_history`.
You should  have `pull request 1917`_ open as you read the commands so you can
see the edits made to the code.

There is also a very exciting terminal-cast at :ref:`terminal_cast`.

Before you begin
================

Make sure you have a local copy of astropy set up as described in
:ref:`get_devel`. In a nutshell, the output of ``git remote -v``, run in the
directory where your local of Astropy resides, should be something like this::

    astropy   git://github.com/astropy/astropy.git (fetch)
    astropy   git://github.com/astropy/astropy.git (push)
    your-user-name     git@github.com:your-user-name/astropy.git (fetch)
    your-user-name     git@github.com:your-user-name/astropy.git (push)

The precise form of the URLs for ``your-user-name`` depends on the
authentication method you set up with GitHub.

The important point is that ``astropy`` should point to the official Astropy
repo and ``your-user-name`` should point to *your* copy of Astropy on GitHub.


Grab the latest updates to astropy
==================================

A few steps in this tutorial take only a single command. They are broken out
separately to outline the process in words as well as code.

Inform your local copy of Astropy about the latest changes in the development
version with::

    git fetch astropy

Set up an isolated workspace
============================
+ Make a new `git`_ branch for fixing this issue, check it out, and let my
  GitHub account know about this branch::

    git branch fix-1761 astropy/master   # branch based on latest from GitHub
    git checkout fix-1761                   # switch to this branch
    git push --set-upstream origin fix-1761 # tell my github acct about it

+ Make a python environment just for this fix and switch to that environment.
  The example below shows the necessary steps in the Anaconda python
  distribution::

        conda create -n apy-1761 --clone root  # Anaconda distribution only
        source activate apy-1761

  If you are using a different distribution, see :ref:`virtual_envs` for
  instructions for creating and activating a new environment.
+ Install our branch in this environment with::
  
    pip install -e .
  
  Remember also to use the proper version of pip or python in this context.

Do you really have to set up a separate python environment for each fix? No,
but you definitely want to have a python environment for your work on code
contributions. Making new environments is fast, doesn't take much space and
provide a way to keep your work organized.

Test first, please
==================

It would be hard to overstate the importance of testing to Astropy. Tests are
what gives you confidence that new code does what it should and that it
doesn't break old code.

You should at least run the relevant tests before you make any changes to make
sure that your python environment is set up properly.

The first challenge is figuring out where to look for relevant tests. `Issue
1761`_ is a problem in the `~astropy.coordinates` package, so the tests for
it are in ``astropy/coordinates/tests``. The rest of Astropy has a similar
layout, described at :ref:`testing-guidelines`.

Change to that directory and run the current tests with::

    cd astropy/coordinates/test
    pytest

The tests all pass, so I need to write a new test to expose this bug.


There are several files with tests in them, though::

    $ ls
    test_angles.py
    test_angular_separation.py
    test_api.py
    test_arrays.py
    test_distance.py
    test_formatting.py
    test_matching.py
    test_name_resolve.py
    test_transformations.py

`Issue 1761`_ affects arrays of coordinates, so it seems sensible to put the
new test in ``test_arrays.py``. As with all of the steps, if you are not
sure, ask on the `astropy-dev mailing list`_.

The goal at this point may be a little counter-intuitive: write a test that we
know will fail with the current code. This test allows Astropy to check,
in an automated way, whether our fix actually works and to make sure future
changes to code do not break our fix.

Looking over the existing code in ``test_arrays.py``, each test is a function
whose name starts with ``test_``; the last test in the file is
``test_array_indexing`` so an appropriate place to add the test is right after
that.

Give the test a reasonably clear name; I chose: ``test_array_len``. The
easiest way to figure out what you need to import and how to set up the test
is to look at other tests. The full test is in the traceback below and in
`pull request 1917`_

Write the test, then see if it works as expected--remember, in this case we
expect to *fail*. Running ``pytest test_arrays.py`` gives the expected
result; an excerpt from the output is::

    ================= FAILURES =============================
    ______________ test_array_len __________________________

        def test_array_len():
            from .. import ICRS

            input_length = 5
            ra = np.linspace(0, 360, input_length)
            dec = np.linspace(0, 90, input_length)

            c = ICRS(ra, dec, unit=(u.degree, u.degree))

    >       assert len(c) == input_length
    E       TypeError: object of type 'ICRS' has no len()

    test_arrays.py:291: TypeError

Success!

Add this test to your local `git`_ repo
=======================================

Keep `git`_ commits small and focused on one logical piece at a time. The test
we just wrote is one logical change, so we will commit it. You could, if you
prefer, wait and commit this test along with your fix.

For this tutorial I'll commit the test separately. If you aren't sure what to
do, ask on `astropy-dev mailing list`_.

Check what was changed
----------------------

We can see what has changed with ``git status``::

    $ git status
    On branch fix-1761
    Your branch is up-to-date with 'origin/fix-1761'.

    Changes not staged for commit:
      (use "git add <file>..." to update what will be committed)
      (use "git checkout -- <file>..." to discard changes in working directory)

        modified:   test_arrays.py

    no changes added to commit (use "git add" and/or "git commit -a")

There are two bits of information here:

+ one file changed, ``test_arrays.py``
+ We have not added our changes to git yet, so it is listed under ``Changes
  not staged for commit``.

For more extensive changes it can be useful to use ``git diff`` to see what
changes have been made::

    $ git diff
    diff --git a/astropy/coordinates/tests/test_arrays.py b/astropy/coordinates/test
    index 2785b59..7eecfbb 100644
    --- a/astropy/coordinates/tests/test_arrays.py
    +++ b/astropy/coordinates/tests/test_arrays.py
    @@ -278,3 +278,14 @@ def test_array_indexing():
         assert c2.equinox == c1.equinox
         assert c3.equinox == c1.equinox
         assert c4.equinox == c1.equinox
    +
    +def test_array_len():
    +    from .. import ICRS
    +
    +    input_length = 5
    +    ra = np.linspace(0, 360, input_length)
    +    dec = np.linspace(0, 90, input_length)
    +
    +    c = ICRS(ra, dec, unit=(u.degree, u.degree))
    +
    +    assert len(c) == input_length

A graphical interface to git makes keeping track of these sorts of changes
even easier; see :ref:`git_gui_options` if you are interested.

Stage the change
----------------

`git`_ requires you to add changes in two steps:

+ stage the change with ``git add test_arrays.py``; this adds the file to
  the list of items that will be added to the repo when you are ready to
  commit.
+ commit the change with ``git commit``; this actually adds the changes to
  your repo.

These can be combined into one step; the advantage of doing it in two steps
is that it is easier to undo staging than committing. As we will see later,
``git status`` even tells you how to do it.

Staging can be very handy if you are making changes in a couple of different
places that you want to commit at the same time. Make your first changes,
stage it, then make your second change and stage that. Once everything is
staged, commit the changes as one commit.

In this case, first stage the change::

    git add test_arrays.py

You get no notice at the command line that anything has changed, but
``git status`` will let you know::

    $ git status
    On branch fix-1761
    Your branch is up-to-date with 'origin/fix-1761'.

    Changes to be committed:
      (use "git reset HEAD <file>..." to unstage)

        modified:   test_arrays.py

Note that `git`_ helpfully includes the command necessary to unstage the
change if you want to.

Commit your change
------------------

I prefer to make commits frequently, so I'll commit the test without the fix::

    $ git commit -m'Add test for array coordinate length (issue #1761)'
    [fix-1761 dd4ef8c] Add test for array coordinate length (issue #1761)
     1 file changed, 11 insertions(+)

Commit messages should be short and descriptive. Including the GitHub issue
number allows GitHub to automatically create links to the relevant issue.

Use ``git status`` to get a recap of where we are so far::

    $ git status
    On branch fix-1761
    Your branch is ahead of 'origin/fix-1761' by 1 commit.
      (use "git push" to publish your local commits)

    nothing to commit, working directory clean

In other words, we have made a change to our local copy of astropy but we
have not pushed (transferred) that change to our GitHub account.

Fix the issue
=============

Write the code
--------------

Now that we have a test written, we'll fix the issue. A full discussion of
the fix is beyond the scope of this tutorial, but the fix is to add a
``__len__`` method to ``astropy.coordinates.SphericalCoordinatesBase`` in
``coordsystems.py``. All of the spherical coordinate systems inherit from
this base class and it is this base class that implements the
``__getitem__`` method that allows indexing of coordinate arrays.

See `pull request 1917`_ to view the changes to the code.

.. _test_changes:

Test your change
----------------

There are a few levels at which you want to test:

+ Does this code change make the test we wrote succeed now? Check
  by running ``pytest tests/test_arrays.py`` in the ``coordinates``
  directory. In this case, yes!
+ Do the rest of the coordinate tests still pass? Check by running ``pytest``
  in the ``coordinates`` directory. In this case, yes--we have not broken
  anything!
+ Do all of the astropy tests still succeed? Check by moving to the top level
  directory (the one that contains ``setup.py``) and run ``python setup.py
  test``. This may take several minutes depending on the speed of your system.
  Success again!

.. note::
    Tests that are skipped or xfailed are fine. A fail or an error is not
    fine. If you get stuck, ask on `astropy-dev mailing list`_ for help!

Stage and commit your change
----------------------------

Add the file to your git repo in two steps: stage, then commit.

To make this a little different than the commit we did above, make sure you
are still in the top level directory and check the ``git status``::

    $ git status
    On branch fix-1761
    Your branch is ahead of 'origin/fix-1761' by 1 commit.
      (use "git push" to publish your local commits)

    Changes not staged for commit:
      (use "git add <file>..." to update what will be committed)
      (use "git checkout -- <file>..." to discard changes in working directory)

        modified:   astropy/coordinates/coordsystems.py

    no changes added to commit (use "git add" and/or "git commit -a")

Note that git knows what has changed no matter what directory you are in (as
long as you are in one of the directories in the repo, that is).

Stage the change with::

    git add astropy/coordinates/coordsystems.py

For this commit it is helpful to use a multi-line commit message that will
automatically close the issue on GitHub when this change is accepted. The
snippet below accomplishes that in bash (and similar shells)::

    $ git commit -m"
    > Add len() to coordinates
    >
    > Closes #1761"
    [fix-1761 f196771] Add len() to coordinates
     1 file changed, 4 insertions(+)

If this was not a tutorial I would write the commit message in a git gui or
run ``git commit`` without a message and git would put me in an editor.

However you do it, the message after committing should look like this::

    Add len() to coordinates

    Closes #1761

You can check the commit messages by running ``git log``. If the commit
message doesn't look right, ask about fixing it at `astropy-dev mailing list`_.

Push your changes to your GitHub fork of astropy
================================================

This one is easy: ``git push``

This copies the changes made on your computer to your copy of Astropy on
GitHub. At this point none of the Astropy maintainers know anything about
your change.

We'll take care of that in a moment with a "pull request", but first...

.. _astropy-fix-add-tests:

Stop and think: any more tests or other changes?
================================================

It never hurts to pause at this point and review whether your proposed
changes are complete. In this case I realized there were some tests I could
have included but didn't:

+ What happens when ``len()`` is called on a coordinate that is *not* an
  array?
+ Does ``len()`` work when the coordinate is an array with one entry?

Both of these are mentioned in the pull request so it doesn't hurt to check
them. In this case they also provide an opportunity to illustrate a feature
of the `pytest`_ framework.

I'll move back to the directory containing the tests with
``cd astropy/coordinates/tests`` to make it a bit easier to run just the test
I want.

The second case is easier, so I'll handle that one first following the cycle
we used above:

+ Make the change in ``test_arrays.py``
+ Test the change

The test passed; rather than committing this one change I'll also implement
the check for the scalar case.

One could imagine two different desirable outcomes here:

+ ``len(scalar_coordinate)`` behaves just like ``len(scalar_angle)``, raising
  a `TypeError` for a scalar coordinate.
+ ``len(scalar_coordinate)`` returns 1 since there is one coordinate.

If you encounter a case like this and are not sure what to do, ask. The best
place to ask is in GitHub on the page for the issue you are fixing.

Alternatively, make a choice and be clear in your pull request on GitHub what
you chose and why; instructions for that are below.

Testing for an expected error
-----------------------------

In this case I opted for raising a `TypeError`, because
the user needs to know that the coordinate they created is not going to
behave like an array of one coordinate if they try to index it later on. It
also provides an opportunity to demonstrate a test when the desired result
is an error.

The `pytest`_ framework makes testing for an exception relatively
easy; you put the code you expect to fail in a ``with`` block::

    with pytest.raises(TypeError):
        c = ICRS(0, 0, unit=(u.degree, u.degree))
        len(c)

I added this to ``test_array_len`` in ``test_arrays.py`` and re-ran the test
to make sure it works as desired.

Aside: Python lesson--let others do your work
---------------------------------------------

The actual fix to this issue was very, very short. In ``coordsystems.py`` two
lines were added::

    def __len__(self):
        return len(self.lonangle)

``lonangle`` contains the ``Angle``s that represent longitude (sometimes this
is an RA, sometimes a longitude). By simply calling ``len()`` on one of the
angles in the array you get, for free, whatever behavior has been defined in
the ``Angle`` class for handling the case of a scalar.

Adding an explicit check for the case of a scalar here would have the very
big downside of having two things that need to be kept in sync: handling of
scalars in ``Angle`` and in coordinates.

Commit any additional changes
=============================

Follow the cycle you saw above:

+ Check that **all** Astropy tests still pass; see :ref:`test_changes`
+ ``git status`` to see what needs to be staged and committed
+ ``git add`` to stage the changes
+ ``git commit`` to commit the changes

The `git`_ commands, without their output, are::

    git status
    git add astropy/coordinates/tests/test_arrays.py
    git commit -m"Add tests of len() for scalar coordinate and length 1 coordinate"

Edit the changelog
==================

Keeping the list of changes up to date is nearly impossible unless each
contributor makes the appropriate updates as they propose changes.

Changes are in the file ``CHANGES.rst`` in the top-level directory (the
directory where ``setup.py`` is). Put the change under the list that matches
the milestone (aka release) that is set for the issue in GitHub. If you are
proposing a new feature in a pull request you may need to wait on this change
until the pull request is discussed.

This issue was tagged for 0.3.1, as shown in the image below, so the changelog
entry went there.

    .. image:: milestone.png

The entry in ``CHANGES.rst`` should summarize was you did and include the
pull request number. For writing changelog entries you don't need to know
much about the markup language being used (though you can read as much as
you want about it at the `Sphinx primer`_); look at other entries and
imitate.

For this issue the entry was the line that starts ``- Implemented``::

    astropy.coordinates
    ^^^^^^^^^^^^^^^^^^^

    - Implemented ``len()`` for coordinate objects. [#1761]

Starting the line with a ``-`` makes a bulleted list item, indenting it makes
it a sublist of ``astropy.coordinates`` and putting ``len()`` in double
backticks makes that text render in a monospace font.

Commit your changes to the CHANGES.rst
--------------------------------------

You can use ``git status`` as above or jump right to staging and committing::

    git add CHANGES.rst
    git commit -m"Add changelog entry for 1761"


Push your changes to GitHub
===========================

One last push to GitHub with these changes before asking for the changes to
be reviewed::

    git push

Ask for your changes to be merged with a pull request
=====================================================

This stage requires to go to your GitHub account and navigate to *your* copy
of astropy; the url will be something like
``https://github.com/your-user-name/astropy``.

Once there, select the branch that contains your fix from the branches
dropdown:

    .. image:: worked_example_switch_branch.png

After selecting the correct branch click on the "Pull Request" button, like
that in the image below:

    .. image:: pull_button.png

Name your pull request something sensible. Include the issue number with a
leading ``#`` in the description of the pull request so that a link is
created to the original issue.

Please see `pull request 1917`_ for the pull request from this example.

Revise and push as necessary
============================

You may be asked to make changes in the discussion of the pull request. Make
those changes in your local copy, commit them to your local repo and push them
to GitHub. GitHub will automatically update your pull request.

.. _Issue 1761: https://github.com/astropy/astropy/issues/1917
.. _pull request 1917: https://github.com/astropy/astropy/issues/1917
.. _Sphinx primer: http://sphinx-doc.org/rest.html
.. _test commit: https://github.com/mwcraig/astropy/commit/cf7d5ac15d7c63ae28dac638c6484339bac5f8de
