.. _contributing_quickstart:

=========================
Contributing Quickstart
=========================

.. _contributing_environment:

Creating a development environment
==================================

To make and test code changes and build the documentation locally you will need to
create a development environment. If you run into problems at any stage do not hesitate
to `ask for help <https://www.astropy.org/help.html>`_.

Set up GitHub and Git
---------------------

Astropy is hosted on `GitHub <https://www.github.com/astropy/astropy>`_, and to
contribute, you will need a `GitHub account
<https://docs.github.com/en/get-started/start-your-journey/creating-an-account-on-github>`_.

We use `Git <https://git-scm.com/>`_ for version control and to allow many people to
work together on the project. See the `GitHub quickstart instructions
<https://docs.github.com/en/get-started/quickstart/set-up-git>`__ for installing and
configuring git, as well as the :ref:`git-resources` page.

If you are new to contributing to projects through forking on GitHub, see the
`GitHub documentation for contributing to projects
<https://docs.github.com/en/get-started/quickstart/contributing-to-projects>`_.

Install a C compiler if needed
------------------------------

How to do this will depend on your platform.

**Windows**

You will need `Build Tools for Visual Studio
<https://visualstudio.microsoft.com/downloads/?q=build+tools>`_.

.. note::
        You DO NOT need to install Visual Studio.
        You only need "Build Tools for Visual Studio" found by
        scrolling down to "All downloads" -> "Tools for Visual Studio" -> "Build Tools
        for Visual Studio".

Alternative options include:

- Install the necessary components on the command line using `vs_BuildTools.exe <https://learn.microsoft.com/en-us/visualstudio/install/use-command-line-parameters-to-install-visual-studio?source=recommendations&view=vs-2022>`_.
- Use the `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_.

**MacOS**

Install the Developer Tools using ``xcode-select --install``. There is no need to
install the full Xcode application and this command will install only the command line
tools and developer utilities.

Further details and related information can be found at
https://devguide.python.org/setup/#macos.

**Linux**

For Linux-based installations, you won't have to install any additional components.

.. _contributing.forking:

Create a clone of astropy
-------------------------

If you have not done so already, you will need your own copy of ``astropy`` to
build it and/or contribute to the source. Astropy is hosted in the `astropy GitHub repository <https://www.github.com/astropy/astropy>`_ and you need to make a clone.

First, create a `GitHub Fork
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ by going to the `astropy project page <https://github.com/astropy/astropy>`_
and hitting the ``Fork`` button.

Next, `clone <https://git-scm.com/docs/git-clone>`__ your GitHub fork to your machine:

.. code-block:: shell

    git clone https://github.com/YOUR-USER-NAME/astropy.git
    cd astropy
    git remote add upstream https://github.com/astropy/astropy.git
    git fetch upstream --tags

This creates the directory ``astropy`` and connects your repository to the upstream
(main project) `astropy <https://github.com/astropy/astropy>`__ repository.

You can see the remote repositories as follows::

    git remote --verbose

You will see something like::

    origin  git@github.com:YOUR-USER-NAME/astropy.git (fetch)
    origin  git@github.com:YOUR-USER-NAME/astropy.git (push)
    upstream        https://github.com/astropy/astropy.git (fetch)
    upstream        https://github.com/astropy/astropy.git (push)

.. _create-isolated-env:

Create an isolated development environment
------------------------------------------

A key requirement is to have an isolated Python environment, meaning that it is
isolated from both your system Python and any other Python environments you may have
for doing other work. This is important because the development environment will often
be unstable and possibly broken at times, and you don't want to break your other work.

There are many good options for doing this, including a number of virtual environment
managers (e.g., the Python standard library `venv <https://docs.python.org/3/library/venv.html>`_
module). Users who have a preference for a particular virtual environment manager are
encouraged to use it!

For this quickstart guide we use the `conda <https://docs.conda.io/en/latest/>`_ package
manager provided by `miniforge <https://github.com/conda-forge/miniforge>`_. This is a
popular choice and generally works well, especially for newcomers. It is easy to install
and use on all platforms and it makes it easy to install different Python versions which
can be useful for testing.

Install miniforge and conda
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you do not already have ``conda`` installed, `download and install miniforge
<https://github.com/conda-forge/miniforge/blob/main/README.md>`_. The details depend on
your system but the end result is to provide a ``conda`` executable that you can use
to create and manage isolated Python environments.

Now create and activate an ``astropy-dev`` conda environment using the following::

   conda create -n astropy-dev python graphviz
   conda activate astropy-dev

Note the ``graphviz`` package is required for building the documentation.

Install the development version of astropy
------------------------------------------

Now you can install the development version of astropy into your new environment. This
will install the latest version of astropy from your local git repo, along with
all the dependencies needed to build and fully test astropy::

   python -m pip install --editable '.[dev_all]'

**Checking the build**

At this point you should be able to import astropy from your locally built version::

   python -c 'import astropy; astropy.system_info()'

Next you may want to try running some or all of the ``astropy`` unit tests.
Running the full test suite can take a few minutes, so you may want to start with a
single sub-package (e.g. :ref:`astropy-coordinates`)::


   # run a sub set of the test suite
   pytest astropy/coordinates

   # or the whole suite
   pytest

Details on running and writing tests can be found in the :ref:`testing-guidelines`
section.

.. _contributing.pre-commit:

Install pre-commit
------------------

This is optional, but *highly recommended*. `Pre-commit <https://pre-commit.com/>`_ is a
tool that runs a number of :ref:`Continuous Integration (CI) <contributing.ci>` checks
(e.g. code formatting) on your code before you commit it. If you skip this step then it
is likely that one or more of those CI checks will fail when you make a pull request,
resulting in lost time (yours and CI resources).

Installation is straightforward. From the root of the astropy repository, run::

    pre-commit install

Now all of the styling checks will be run each time you commit changes, ensuring that
the CI formatting checks for your :ref:`pull request <quickstart-pull-request>` will
pass.

.. tip:: To learn more about pre-commit, see the :ref:`pre-commit` section.

.. _contributing.pull_request:

Creating and submitting a pull request
======================================

You can contribute bug fixes, new features, and documentation updates by submitting a
GitHub pull request (PR). This section will guide you through the process. We encourage
you to `ask for help <https://www.astropy.org/help.html>`_ if you get stuck. The Astropy
community is welcoming and friendly and will help you!

If you are new to the Astropy Project and interested to submit a large patch
(e.g., a new big feature or significant refactoring), we encourage you to first
discuss your ideas on GitHub to increase the chance of your PR
being accepted.

Creating a branch
-----------------

Your local ``main`` branch should always reflect the current state of astropy repository.
First ensure it's up-to-date with the main astropy repository::

    git switch main
    git pull upstream main --ff-only

Now create a development branch for making your changes. For example::

    git switch -c subpackage-bug-fix

This changes your working branch from ``main`` to the ``subpackage-bug-fix`` branch.
Keep any changes in this branch specific to one bug or feature so it is clear what the
branch brings to astropy. You can have many feature branches and switch in between them
using the `git switch <https://git-scm.com/docs/git-switch>`_ command.

Using a descriptive branch name can help you stay organized. For example
```io-ascii-commented-header``` might be a good name for a branch that fixes the
commented header issue `#15513 <https://github.com/astropy/astropy/issues/15513>`_ in
the ``io.ascii`` sub-package.

When you want to update the feature branch with changes in main after
you created the branch, check the section on
:ref:`updating a PR <contributing.update-pr>`.

.. _contributing.commit-code:

Making code or documentation changes
------------------------------------

Now comes the fun part where you use your favorite editor or IDE to make changes to the
code or documentation! At a high level this breaks into a few parts:

- **Make changes**: Make the changes you want to make. This could be fixing a bug,
  adding a new feature, or updating the documentation.
- **Test changes**: For code changes, ensure that they work as expected following the
  process outlined in the :ref:`testing-guidelines` section.
- **Build documentation**: If you are updating the documentation, you will want to
  :ref:`build the documentation <builddocs>` to ensure that it looks good.
- **Add a changelog entry**: For most code changes you will need to
  :ref:`add-changelog`.

.. tip:: For more information and examples see :ref:`edit-flow` section.

You can see a summary of the changes you've currently made by running:

.. code-block:: shell

    git status

You can then commit your all your changes to your local repository with an explanatory
`commit message <https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html>`_:

.. code-block:: shell

    git add files-that-you-changed ...
    git commit -m "your commit message goes here"

.. Important:: Never merge changes from ``upstream/main`` into your feature branch. If
   changes in ``main`` require changes to our code you must :ref:`rebase`.

.. _contributing.push-code:

Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits::

    git push origin --set-upstream subpackage-bug-fix

Here ``origin`` is the default name given to your fork on GitHub.

Now your code is on GitHub, but it is not visible to the Astropy maintainers. For that
to happen, a pull request needs to be submitted on GitHub.

The first time you push to a new branch on GitHub, you will see a message like below
with a useful link to create a pull request::

  remote: Create a pull request for 'subpackage-bug-fix' on GitHub by visiting:
  remote:      https://github.com/YOUR-USER-NAME/astropy/pull/new/subpackage-bug-fix


.. _quickstart-pull-request:

Making a pull request
---------------------

If everything looks good, you are ready to make a pull request (PR). A PR is how
code from your local repository becomes available to the GitHub community to review and
merged into project to appear the in the next release.

Most of the time you can just follow the link that ``git`` provided when you pushed
your branch and create the PR. If you don't have that link (and for a few more details), you can follow the :ref:`pull-request` instructions.

Follow the instructions in the PR template and fill it out as completely as possible.

If your PR is still a work in progress then instead of clicking "Create pull request",
click on the small down arrow next to it and select "`Create draft pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests#draft-pull-requests>`__".
In addition, if your commits are not ready for CI testing, you
should include ``[ci skip]`` the last commit message – but note that code formatting
checks and documentation building will still be done. Formatting and style errors *should*
already have been fixed before committing if you have locally
:ref:`installed pre-commit<contributing.pre-commit>`; but if you have not,
you can use the :ref:`pre-commit_bot` to fix them automatically in the PR.

Once submitted (and marked as ready), this request goes to the astropy maintainers and
they will review the PR.

.. _contributing.update-pr:

Updating your pull request
--------------------------

Based on the review you get on your pull request, you will probably need to make
some adjustments. You can follow the :ref:`code committing steps <contributing.commit-code>`
again to address any feedback and update your pull request::

    git push origin subpackage-bug-fix

Any ``git push`` will automatically update your pull request with your branch's changes
and restart the :ref:`Continuous Integration <contributing.ci>` checks.

.. Important:: At this point please read (or at least skim) the sections :ref:`revise
    and push`, :ref:`rebase`, and :ref:`squash-if-necessary`. The information here
    covers situations that happen on occasion and can be cause trouble. As always if
    you have questions, ask for help from the maintainer reviewing your PR.

Tips for a successful pull request
----------------------------------

If you have made it to this point and submitted a pull request, one of the core
maintainers will take a look. To make the process as smooth and efficient as possible,
here are some tips:

- **Reference any existing open issue** to `link to that issue <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests#draft-pull-requests>`_ and close the
  issue if the PR is merged.
- **Ensure you have appropriate tests**.
- **Keep your pull requests as simple as possible** -- larger PRs take longer to review.
- **When practical, limit the scope of a PR to one sub-package** -- this means fewer
  required reviewers and a faster review process.
- **Ensure that CI is in a green state** -- any required failures should be addressed.
