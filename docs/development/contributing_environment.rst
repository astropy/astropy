.. _contributing_environment:

==================================
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
In addition, using ``pre-commit`` will also allow you to more easily
remain up-to-date with our code checks as they change.

Note that if needed, you can skip these checks with ``git commit --no-verify``.

Pre-commit DETAILS THAT NEED TO GO SOMEWHERE ELSE
-------------------------------------------------

If you don't want to use ``pre-commit`` as part of your workflow, you can still use it
to run its checks with one of the following::

    pre-commit run --files <files you have modified>
    pre-commit run --from-ref=upstream/main --to-ref=HEAD --all-files

without needing to have done ``pre-commit install`` beforehand.

Finally, we also have some slow pre-commit checks, which don't run on each commit
but which do run during continuous integration. You can trigger them manually with::

    pre-commit run --hook-stage manual --all-files

.. note::

    You may want to periodically run ``pre-commit gc``, to clean up repos
    which are no longer used.

.. note::

    If you have conflicting installations of ``virtualenv``, then you may get an
    error - see `here <https://github.com/pypa/virtualenv/issues/1875>`_.

    Also, due to a `bug in virtualenv <https://github.com/pypa/virtualenv/issues/1986>`_,
    you may run into issues if you're using conda. To solve this, you can downgrade
    ``virtualenv`` to version ``20.0.33``.
