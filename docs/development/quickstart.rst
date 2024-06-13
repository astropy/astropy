.. _development_quickstart:

=========================
Development Quickstart
=========================

.. _contributing_environment:

Creating a development environment
==================================

To make and test code changes and build the documentation locally you will need to
create a development environment. If you run into problems at any stage do not hesitate
to `ask for help <https://www.astropy.org/help.html>`_.

Set up GitHub and Git
---------------------

Astropy is hosted in a `astropy GitHub repository
<https://www.github.com/astropy/astropy>`_, and to contribute, you will need to sign up
for a `free GitHub account <https://github.com/signup/free>`_.

To contribute to the astropy repository, you will need to ensure that `Git
<https://git-scm.com/>`_ is installed and configured on your machine. `GitHub has
instructions <https://docs.github.com/en/get-started/quickstart/set-up-git>`__ for
installing and configuring git.

More details and further resources are available in the
:ref:`contributing.version_control` section.

Install a C compiler
--------------------

How to do this will depend on your platform.

**Windows**

You will need `Build Tools for Visual Studio 2022
<https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022>`_.

.. note::
        You DO NOT need to install Visual Studio 2022.
        You only need "Build Tools for Visual Studio 2022" found by
        scrolling down to "All downloads" -> "Tools for Visual Studio".
        In the installer, select the "Desktop development with C++" Workloads.

Alternative options include:

- Install the necessary components on the command line using `vs_BuildTools.exe <https://learn.microsoft.com/en-us/visualstudio/install/use-command-line-parameters-to-install-visual-studio?source=recommendations&view=vs-2022>`_.
- Use the `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_.

**MacOS**

Install the Developer Tools using ``xcode-select --install``.

Further details and related information can be found at
https://devguide.python.org/setup/#macos.

**Linux**

For Linux-based :ref:`conda <contributing.conda>` installations, you won't have to
install any additional components.

. _contributing.forking:

Create a fork and clone of astropy
----------------------------------

If you have not done so already, you will need your own copy of ``astropy`` to
work on the code.

First, create a `GitHub Fork
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ by going to the `astropy project page <https://github.com/astropy/astropy>`_
and hitting the ``Fork`` button.

Next, clone your GitHub fork to your machine:

.. code-block:: shell

    git clone https://github.com/YOUR-USER-NAME/astropy.git
    cd astropy
    git remote add upstream https://github.com/astropy/astropy.git
    git fetch upstream

This creates the directory ``astropy`` and connects your repository to the upstream
(main project) `astropy <https://github.com/astropy/astropy>`_ repository.

Create an isolated development environment
------------------------------------------

A key requirement is to have an isolated Python environment, meaning that it is
isolated from both your system Python and any other Python environments you may have
for doing other work. This is important because the development environment will often
be unstable and possibly broken at times, and you don't want to break your other work.

There are *many* good options (see :ref:`virtual_envs` for discussion), but in this
quickstart guide we use the `conda <https://docs.conda.io/en/latest/>`_ package
manager provided by `miniforge <https://github.com/conda-forge/miniforge>`_. This is a
popular choice and generally works well, especially for newcomers. It is easy to install
and use on all platforms and it makes it easy to install different Python versions which
can be useful for testing.

.. _contributing.conda:

Install miniforge and conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you do not already have ``conda`` installed, `download and install miniforge
<https://github.com/conda-forge/miniforge/blob/main/README.md>`_. The details depend on
your system but the end result is to provide a ``conda`` executable that you can use
to create and manage isolated Python environments.

Now create and activate an ``astropy-dev`` conda environment using the following:

.. code-block:: shell

   conda create -n astropy-dev python graphviz
   conda activate astropy-dev

Note the ``graphviz`` package is required for building the documentation.

Install the development version of astropy
------------------------------------------

Now you can install the development version of astropy into your new environment. This
will install the latest version of astropy from your local git repo, along with
all the dependencies needed to build and fully test astropy.

.. code-block:: shell

   python -m pip install --editable ".[dev_all]"

**Checking the build**

At this point you should be able to import astropy from your locally built version::

   python
   >>> import astropy
   >>> astropy.__version__  # note: the exact output will differ
   '7.0.0.dev303+gb394fda545.d20240613'

At this point you may want to try running some or all of the ``astropy`` unit tests.
Running the full test suite can take a while, so you may want to start with a subset
of only the coordinates tests::

.. code-block:: shell

   pytest astropy/coordinates
   pytest

Details on running and writing tests can be found in the :ref:`testing-guidelines`
section.

.. _contributing.pre-commit:

Install pre-commit
------------------

This is optional, but *highly recommended*. Pre-commit is a tool that runs a number of
:ref:`Continuous Integration (CI) <contributing.ci>` checks (e.g. code formatting) on
your code before you commit it. If you skip this step then it is likely that one or more
of those CI checks will fail when you make a pull request, resulting in lost time (yours
and CI resources).

Installation is straightforward. From the root of the astropy repository, run::

    pre-commit install

Now all of the styling checks will be
run each time you commit changes without your needing to run each one manually.

.. tip:: To learn more about pre-commit, see the :ref`pre-commit` section.

.. _contributing.pull_request:

Creating and submitting a pull request
======================================

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

* `GitHub documentation for forking a repo <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.
* `GitHub documentation for collaborating with pull requests <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests>`_.
* `GitHub documentation for working with forks <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks>`_.

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
----------------------------------

If you have made it to this point and submitted a pull request, one of the core
maintainers will take a look. To make the process as smooth and efficient as possible,
here are some tips:

- **Reference an open issue** for non-trivial changes to clarify the PR's purpose.
- **Ensure you have appropriate tests**.
- **Keep your pull requests as simple as possible** -- larger PRs take longer to review.
- **Ensure that CI is in a green state** -- any required failures should be addressed.
