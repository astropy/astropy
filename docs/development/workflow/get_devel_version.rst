.. _get_devel:

***************************
Try the development version
***************************

.. note::
    `git`_ is the name of a source code management system. It is used to keep
    track of changes made to code and to manage contributions coming from
    several different people. If you want to read more about `git`_ right now
    take a look at `Git Basics`_.

    If you have never used `git`_ before, allow one hour the first time you do
    this. If you find this taking more than one hour, post in one of the
    `astropy forums <http://www.astropy.org/help.html>`_ to get help.


Trying out the development version of astropy is useful in three ways:

* More users testing new features helps uncover bugs before the feature is
  released.
* A bug in the most recent stable release might have been fixed in the
  development version. Knowing whether that is the case can make your bug
  reports more useful.
* You will need to go through all of these steps before contributing any
  code to `Astropy`_. Practicing now will save you time later if you plan to
  contribute.

Overview
========

Conceptually, there are several steps to getting a working copy of the latest
version of astropy on your computer:

#. :ref:`fork_a_copy`; this copy is called a *fork* (if you don't have an
   account on `github`_ yet, go there now and make one).
#. :ref:`check_git_install`
#. :ref:`clone_your_fork`; this is called making a *clone* of the repository.
#. :ref:`set_upstream_main`
#. :ref:`make_a_branch`; this is called making a *branch*.
#. :ref:`activate_development_astropy`
#. :ref:`test_installation`
#. :ref:`try_devel`
#. :ref:`deactivate_development`

Step-by-step instructions
=========================

.. _fork_a_copy:

Make your own copy of Astropy on GitHub
---------------------------------------

In the language of `GitHub`_, making a copy of someone's code is called making
a *fork*. A fork is a complete copy of the code and all of its revision
history.

#. Log into your `GitHub`_ account.

#. Go to the `Astropy GitHub`_ home page.

#. Click on the *fork* button:

   .. image:: ../workflow/forking_button.png

   After a short pause and an animation of Octocat scanning a book on a
   flatbed scanner, you should find yourself at the home page for your own
   forked copy of astropy.

.. _check_git_install:

Make sure git is installed and configured on your computer
----------------------------------------------------------

**Check that git is installed:**

Check by typing, in a terminal::

    $ git --version
    # if git is installed, will get something like: git version 2.20.1

If `git`_ is not installed, `get it <https://git-scm.com/downloads>`_.

**Basic git configuration:**

Follow the instructions at `Set Up Git at GitHub`_ to take care of two
essential items:

+ Set your user name and email in your copy of `git`_

+ Set up authentication so you don't have to type your github password every
  time you need to access github from the command line. The default method at
  `Set Up Git at GitHub`_ may require administrative privileges; if that is a
  problem, set up authentication
  `using SSH keys instead <https://help.github.com/en/articles/connecting-to-github-with-ssh>`_

We also recommend setting up `git`_ so that when you copy changes from your
computer to `GitHub`_ only the copy (called a *branch*) of astropy that you are
working on gets pushed up to GitHub.  *If* your version of git is 1.7.11 or,
greater, you can do that with::

    git config --global push.default simple

If you skip this step now it is not a problem; `git`_ will remind you to do it in
those cases when it is relevant.  If your version of git is less than 1.7.11,
you can still continue without this, but it may lead to confusion later, as you
might push up branches you do not intend to push.

.. note::

    Make sure you make a note of which authentication method you set up
    because it affects the command you use to copy your GitHub fork to your
    computer.

    If you set up password caching (the default method) the URLs will look like
    ``https://github.com/your-user-name/astropy.git``.

    If you set up SSH keys the URLs you use for making copies will look
    something like ``git@github.com:your-user-name/astropy.git``.


.. _clone_your_fork:

Copy your fork of Astropy from GitHub to your computer
------------------------------------------------------

One of the commands below will make a complete copy of your `GitHub`_ fork
of `Astropy`_ in a directory called ``astropy``; which form you use depends
on what kind of authentication you set up in the previous step::

    # Use this form if you setup SSH keys...
    $ git clone --recursive git@github.com:your-user-name/astropy.git
    # ...otherwise use this form:
    $ git clone --recursive https://github.com/your-user-name/astropy.git

If there is an error at this stage it is probably an error in setting up
authentication.

.. _set_upstream_main:

Tell git where to look for changes in the development version of Astropy
------------------------------------------------------------------------

Right now your local copy of astropy doesn't know where the development
version of astropy is. There is no easy way to keep your local copy up to
date. In `git`_ the name for another location of the same repository is a
*remote*. The repository that contains the latest "official" development
version is traditionally called the *upstream* remote, but here we use a
more meaningful name for the remote: *astropy*.

Change into the ``astropy`` directory you created in the previous step and
let `git`_ know about about the astropy remote::

    cd astropy
    git remote add astropy https://github.com/astropy/astropy.git

You can check that everything is set up properly so far by asking `git`_ to
show you all of the remotes it knows about for your local repository of
`Astropy`_ with ``git remote -v``, which should display something like::

    astropy   https://github.com/astropy/astropy.git (fetch)
    astropy   https://github.com/astropy/astropy.git (push)
    origin     git@github.com:your-user-name/astropy.git (fetch)
    origin     git@github.com:your-user-name/astropy.git (push)

Note that `git`_ already knew about one remote, called *origin*; that is your
fork of `Astropy`_ on `GitHub`_.

To make more explicit that origin is really *your* fork of `Astropy`_, rename that
remote to your `GitHub`_ user name::

  git remote rename origin your-user-name

.. _make_a_branch:

Create your own private workspace
---------------------------------

One of the nice things about `git`_ is that it is easy to make what is
essentially your own private workspace to try out coding ideas. `git`_
calls these workspaces *branches*.

Your repository already has several branches; see them if you want by running
``git branch -a``. Most of them are on ``remotes/origin``; in other words,
they exist on your remote copy of astropy on GitHub.

There is one special branch, called *main*. Right now it is the one you are
working on; you can tell because it has a marker next to it in your list of
branches: ``* main``.

To make a long story short, you never want to work on main. Always work on a branch.

To avoid potential confusion down the road, make your own branch now; this
one you can call anything you like (when making contributions you should use
a meaningful more name)::

    git branch my-own-astropy

You are *not quite* done yet. Git knows about this new branch; run
``git branch`` and you get::

    * main
      my-own-astropy

The ``*`` indicates you are still working on main. To work on your branch
instead you need to *check out* the branch ``my-own-astropy``. Do that with::

    git checkout my-own-astropy

and you should be rewarded with::

    Switched to branch 'my-own-astropy'

.. _activate_development_astropy:

"Activate" the development version of astropy
---------------------------------------------

Right now you have the development version of astropy, but python will not
see it. Though there are more sophisticated ways of managing multiple versions
of astropy, for now this straightforward way will work (if you want to jump
ahead to the more sophisticated method look at :ref:`virtual_envs`).

.. note::
    If you want to work on C or Cython code in `Astropy`_, this quick method
    of activating your copy of astropy will *not* work -- you need to go
    straight to using :ref:`virtual_envs`.

If you have decided to use the recommended "activation" method with
``pip``, please note the following: Before trying to install,
check that you have the required dependency: "cython".
If not, install it with ``pip``. Note that on some platforms,
the pip command is ``pip3`` instead of ``pip``, so be sure to use
this instead in the examples below if that is the case.
If you have any problem with different versions of ``pip`` installed,
try aliasing to resolve the issue. If you are unsure about which ``pip``
version you are using, try the command ``which pip`` on the terminal.

In the directory where your copy of astropy is type::

    pip install -e .[test]

This command installs astropy itself, along with a few packages which will be
useful for testing the changes you will make down the road. Several pages of
output will follow the first time you do this; this would not be a bad time to
get a fresh cup of coffee. At the end of it you should see something like
``Finished processing dependencies for astropy==3.2.dev6272``.

To make sure it has been activated **change to a different directory outside of
the astropy distribution** and try this in python::

    >>> import astropy
    >>> astropy.__version__  # doctest: +SKIP
    '3.2.dev6272'

The actual version number will be different than in this example, but it
should have ``'dev'`` in the name.

.. warning::
    Right now every time you run Python, the development version of astropy
    will be used. That is fine for testing but you should make sure you change
    back to the stable version unless you are developing astropy. If you want
    to develop astropy, there is a better way of separating the development
    version from the version you do science with. That method, using a
    `virtualenv`_ or a `conda`_ environment, is discussed at :ref:`virtual_envs`.

    For now **remember to change back to your usual version** when you are
    done with this.

.. _test_installation:

Test your development copy
--------------------------

Testing is an important part of making sure astropy produces reliable,
reproducible results. Before you try out a new feature or think you have found
a bug make sure the tests run properly on your system.

Before running your tests, please see :ref:`testing-dependencies`.

If the test *don't* complete successfully, that is itself a bug--please
`report it <https://github.com/astropy/astropy/issues>`_.

To run the tests, navigate back to the directory your copy of astropy is in on
your computer, then, at the shell prompt, type::

    pytest

This is another good time to get some coffee or tea. The number of tests is
large. When the test are done running you will see a message something like
this::

    4741 passed, 85 skipped, 11 xfailed

Skips and xfails are fine, but if there are errors or failures please
`report them <https://github.com/astropy/astropy/issues>`_.

.. _try_devel:

Try out the development version
-------------------------------

If you are going through this to ramp up to making more contributions to
`Astropy`_ you don't actually have to do anything here.

If you are doing this because you have found a bug and are checking that it
still exists in the development version, try running your code.

Or, just for fun, try out one of the :ref:`new features <changelog>` in
the development version.

Either way, once you are done, make sure you do the next step.

.. _deactivate_development:

"Deactivate" the development version
------------------------------------

Be sure to turn the development version off before you go back to doing
science work with astropy.

Navigate to the directory where your local copy of the development version is,
then run::

    pip uninstall astropy

This should remove the development version only. Once again,
it is important to check that you are using the proper version of
``pip`` corresponding to the Python executable desired.

You should really confirm it is deactivated by **changing to a different
directory outside of the astropy distribution** and running this in python::

    >>> import astropy
    >>> astropy.__version__  # doctest: +SKIP
    '3.1.1'

The actual version number you see will likely be different than this example,
but it should not have ``'dev'`` in it.


.. include:: links.inc
.. _Git Basics: https://git-scm.com/book/en/Getting-Started-Git-Basics
.. _Set Up Git at GitHub: https://help.github.com/en/articles/set-up-git#set-up-git
