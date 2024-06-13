.. _contributing.pull_request:

**************************************
Creating and submitting a contribution
**************************************

You can contribute bug fixes, new features, and documentation updates by submitting a
GitHub pull request (PR). This section will guide you through the process. We encourage
you to `ask for help <https://www.astropy.org/help.html>`_ if you get stuck. The astropy
community is welcoming and friendly and will help you!

.. _contributing.version_control:

GitHub and Git
--------------

Astropy is hosted on `GitHub <https://www.github.com/astropy/astropy>`_, and to
contribute, you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_.

We use `Git <https://git-scm.com/>`_ for version control and to allow many people to
work together on the project. If you are new to Git then take a moment to look at the
:ref:`git_resources` page.

If you are new to contributing to projects through forking on GitHub, take a look at the
`GitHub documentation for contributing to projects
<https://docs.github.com/en/get-started/quickstart/contributing-to-projects>`_. GitHub
provides a quick tutorial using a test repository that may help you become more familiar
with forking a repository, cloning a fork, creating a feature branch, pushing changes
and making pull requests.

Below are some useful resources for learning more about forking and pull requests on GitHub:

* the `GitHub documentation for forking a repo <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.
* the `GitHub documentation for collaborating with pull requests <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests>`_.
* the `GitHub documentation for working with forks <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks>`_.

Creating a feature branch
-------------------------

Your local ``main`` branch should always reflect the current state of astropy repository.
First ensure it's up-to-date with the main astropy repository::

    git checkout main
    git pull upstream main --ff-only

Now create a feature branch for making your changes. For example::

    git switch -c shiny-new-feature

This changes your working branch from ``main`` to the ``shiny-new-feature`` branch.
Keep any changes in this branch specific to one bug or feature so it is clear what the
branch brings to astropy. You can have many feature branches and switch in between them
using the ``git switch`` command.

Using a descriptive branch name can help you stay organized. For example
`io-ascii-commented-header-15513` might be a good name for a branch that fixes the
commented header issue `#15513 <https://github.com/astropy/astropy/issues/15513>`_ in
the ``io.ascii`` sub-package.

.. Important:: Never use the ``main`` branch for making changes. Always create a new
   feature branch for your changes.

When you want to update the feature branch with changes in main after
you created the branch, check the section on
:ref:`updating a PR <contributing.update-pr>`.

.. _contributing.commit-code:

Making code or documentation changes
------------------------------------

Before modifying any part of the astropy package, ensure you follow the
:ref:`contributing environment <contributing_environment>` guidelines to set up an
appropriate development environment. Next you will find it helpful to scan through the
:ref:`astropy contribution guidelines <contributing_codebase>`.

Once you have made code or documentation changes you will want to test the code
FIXME:LINK and/or build the documentation FIXME:LINK to ensure that your changes work as
expected.

You can see all the changes you've currently made by running:

.. code-block:: shell

    git status

You can then commit your all your changes to your local repository with an explanatory
commit message:

.. code-block:: shell

    git add files-that-you-changed ...
    git commit -m "your commit message goes here"

You should make frequent commits and always include a commit message. Each commit
should represent one logical set of changes.

.. Important:: Never merge changes from ``upstream/main`` into your feature branch. If
   changes in ``main`` require changes to our code you must :ref:`rebase`.

.. _contributing.push-code:

Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits

.. code-block:: shell

    git push origin shiny-new-feature

Here ``origin`` is the default name given to your remote repository on GitHub.
You can see the remote repositories

.. code-block:: shell

    git remote -v

If you added the upstream repository as described above you will see something
like

.. code-block:: shell

    origin  git@github.com:yourname/astropy.git (fetch)
    origin  git@github.com:yourname/astropy.git (push)
    upstream        https://github.com/astropy/astropy.git (fetch)
    upstream        https://github.com/astropy/astropy.git (push)

Now your code is on GitHub, but it is not yet a part of astropy. For that to
happen, a pull request needs to be submitted on GitHub.

Making a pull request
---------------------

If everything looks good, you are ready to make a pull request. A pull request is how
code from your local repository becomes available to the GitHub community to review and
merged into project to appear the in the next release.

To submit a pull request follow the steps outlined in the GitHub documentation `Creating
a pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_.

This request then goes to the repository maintainers, and they will review the code.

.. _contributing.update-pr:

Updating your pull request
--------------------------

Based on the review you get on your pull request, you will probably need to make
some changes to the code. You can follow the :ref:`code committing steps <contributing.commit-code>`
again to address any feedback and update your pull request.


.. code-block:: shell

    git push origin shiny-new-feature

Any ``git push`` will automatically update your pull request with your branch's changes
and restart the :ref:`Continuous Integration <contributing.ci>` checks.

FIXME: reference docs on rebasing if necessary.

Tips for a successful pull request
==================================

If you have made it to this point and submitted a pull request, one of the core
maintainers will take a look. To make the process as smooth and efficient as possible,
here are some tips:

- **Reference an open issue** for non-trivial changes to clarify the PR's purpose.
- **Ensure you have appropriate tests**.
- **Keep your pull requests as simple as possible** -- larger PRs take longer to review.
- **Ensure that CI is in a green state** -- any required failures should be addressed.
