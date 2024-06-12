.. _contributing_environment:

==================================
Creating a development environment
==================================

To make and test code changes and build the documentation locally you will need to create a
development environment. This requires a C/C++ compiler and an isolated Python environment.

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

Alternatively, you can install the necessary components on the commandline using
`vs_BuildTools.exe <https://learn.microsoft.com/en-us/visualstudio/install/use-command-line-parameters-to-install-visual-studio?source=recommendations&view=vs-2022>`_

Alternatively, you could use the `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_
and consult the ``Linux`` instructions below.

**macOS**

To use the :ref:`conda <contributing.conda>`-based compilers, you will need to install the
Developer Tools using ``xcode-select --install``.

If you prefer to use a different compiler, general information can be found here:
https://devguide.python.org/setup/#macos

**Linux**

For Linux-based :ref:`conda <contributing.conda>` installations, you won't have to install any
additional components outside of the conda environment.


FIXME:

Let us know if you have any difficulties by opening an issue or reaching out on our contributor
community :ref:`Slack <community.slack>`.


. _contributing.forking:

Create a fork of astropy
-----------------------

If you have not done so already, you will need your own copy of astropy (aka fork) to
work on the code. Go to the `astropy project page <https://github.com/astropy/astropy>`_
and hit the ``Fork`` button. For more information see the `GitHub Fork documentation
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_.

Next you will want to clone your fork to your machine:

.. code-block:: shell

    git clone https://github.com/your-user-name/astropy.git astropy
    cd astropy
    git remote add upstream https://github.com/astropy/astropy.git
    git fetch upstream

This creates the directory ``astropy`` and connects your repository to
the upstream (main project) *astropy* repository.


Create an isolated development environment
------------------------------------------

A key requirement is to have an isolated Python environment, meaning that it is
isolated from both your system Python and any other Python environments you may have
for doing other work. This is important because the development environment may well
be unstable and possibly broken at times, and you don't want to break your other work.

There are *many* good options (see :ref:`virtual_envs` for discussion), but in this
quickstart guide we will use the `conda <https://docs.conda.io/en/latest/>`_ package
manager provided by `miniforge <https://github.com/conda-forge/miniforge>`_. This is a
popular choice and generally works well, especially for newcomers.

.. _contributing.conda:

Install conda and dev astropy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `Download and install miniforge <https://github.com/conda-forge/miniforge/blob/main/README.md>`_ to install ``conda``.
* Make sure your conda is up to date (``conda update conda``)
* Create and activate the ``astropy-dev`` conda environment using the following commands:

.. code-block:: none

   conda create -n astropy-dev graphviz
   conda activate astropy-dev
   python -m pip install --editable ".[dev_all]"

**Checking the build**

At this point you should be able to import astropy from your locally built version::

   $ python
   >>> import astropy
   >>> print(astropy.__version__)  # note: the exact output may differ
   2.0.0.dev0+880.g2b9e661fbb.dirty


At this point you may want to try
`running the test suite <https://astropy.pydata.org/docs/dev/development/contributing_codebase.html#running-the-test-suite>`_.

**Keeping up to date with the latest build**

When building astropy with meson, importing astropy will automatically trigger a rebuild, even when C/Cython files are modified.
By default, no output will be produced by this rebuild (the import will just take longer). If you would like to see meson's
output when importing astropy, you can set the environment variable ``MESONPY_EDTIABLE_VERBOSE``. For example, this would be::

   # On Linux/macOS
   MESONPY_EDITABLE_VERBOSE=1 python

   # Windows
   set MESONPY_EDITABLE_VERBOSE=1 # Only need to set this once per session
   python

If you would like to see this verbose output every time, you can set the ``editable-verbose`` config setting to ``true`` like so::

   python -m pip install -ve . --config-settings editable-verbose=true

.. tip::
   If you ever find yourself wondering whether setuptools or meson was used to build your astropy,
   you can check the value of ``astropy._built_with_meson``, which will be true if meson was used
   to compile astropy.


Step 4: install pre-commit hooks
---------------------------------

FIXME: This section is not yet complete.

It's recommended to also install the :ref:`pre-commit hooks <contributing.pre-commit>`.
