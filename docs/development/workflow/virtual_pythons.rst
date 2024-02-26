:orphan:

.. include:: links.inc
.. _virtual_envs:

***************************
Python virtual environments
***************************

If you plan to do regular work on Astropy you should do your development in
a Python virtual environment. Conceptually a virtual environment is a
duplicate of the Python environment you normally work in, but sandboxed from
your default Python environment in the sense that packages installed in the
virtual environment do not affect your normal working environment in any way.
This allows you to install, for example, a development version of Astropy
and its dependencies without it conflicting with your day-to-day work with
Astropy and other Python packages.

.. note::

    "Default Python environment" here means whatever Python you are using
    when you log in; i.e. the default Python installation on your system,
    which is not in a Conda environment or virtualenv.

    More specifically, in UNIX-like platforms it creates a parallel root
    "prefix" with its own ``bin/``, ``lib/``, etc. directories.  When you
    :ref:`activate <activate_env>` the virtual environment it places this
    ``bin/`` at the head of your ``$PATH`` environment variable.

    This works similarly on Windows but the details depend on how you
    installed Python and whether or not you're using Miniconda.

There are a few options for using virtual environments; the choice of method
is dictated by the Python distribution you use:

* If you wish to use the `miniconda`_ Python distribution, you must use the
  `conda`_ command to make and manage your virtual environments.

.. note::

    `miniconda`_ is a minimal flavor of the popular Anaconda Python
    distribution, containing only ``conda``, Python, other useful packages
    (such as ``pip``, ``requests``, etc.) and their dependencies in the
    ``base`` environment.
    Further packages and environments can then be bootstrapped from the
    existing ``base`` environment, allowing developers to install *just the
    packages* needed, and nothing more.

* If you do not wish to use Miniconda you can use `virtualenv`_ and the conda-like
  helper commands provided by `virtualenvwrapper`_; you *can not* use this
  with `conda`_. As the name suggests, `virtualenvwrapper`_ is a wrapper
  around `virtualenv`_.

* A third, more recent option which is growing in popularity is `pipenv`_
  which builds on top of `virtualenv`_ to provide project-specific Python
  environments and dependency management.

In both cases you will go through the same basic steps; the commands to
accomplish each step are given for both `conda`_ and `virtualenvwrapper`_:

* :ref:`setup_for_env`
* :ref:`list_env`
* :ref:`create_env`
* :ref:`activate_env`
* :ref:`deactivate_env`
* :ref:`delete_env`

Another well-maintained guide to Python virtualenvs (specifically `pipenv`_
and `virtualenv`_, though it does not discuss `conda`_) which has been
translated into multiple languages is the `Hitchhiker's Guide to Python
<https://docs.python-guide.org/dev/virtualenvs/>`_ chapter on the subject.


.. _setup_for_env:


Set up for virtual environments
===============================

* `conda`_: No setup is necessary beyond installing the Miniconda Python
  distribution . You can find the `conda installation instructions here`_.

* `virtualenvwrapper`_:

  + First, install `virtualenvwrapper`_, which will also install `virtualenv`_,
    with::

        python -m pip install --user virtualenvwrapper

  + From the `documentation for virtualenvwrapper`_, you also need to::

      export WORKON_HOME=$HOME/.virtualenvs
      export PROJECT_HOME=$HOME/
      source /usr/local/bin/virtualenvwrapper.sh

* `pipenv`_: Install the ``pipenv`` command using your default pip (the
  pip in the default Python environment)::

      python -m pip install --user pipenv

.. _list_env:

List virtual environments
=========================

You do not need to list the virtual environments you have created before using
them...but sooner or later you will forget what environments you have defined
and this is the easy way to find out.

* `conda`_: ``conda info -e``
    + you will always have at least one environment, called ``base``.
    + your active environment is indicated by a ``*``

* `virtualenvwrapper`_: ``workon``
    + If this displays nothing you have no virtual environments
    + If this displays ``workon: command not found`` then you haven't done
      the :ref:`setup_for_env`; do that.
    + For more detailed information about installed environments use
      ``lsvirtualenv``.

* `pipenv`_ does not have a concept of listing virtualenvs; it instead
  automatically generates the virtualenv associated with a project directory
  (e.g. the Astropy source repository on your computer).

.. _create_env:

Create a new virtual environment
================================

This needs to be done once for each virtual environment you want. There is one
important choice you need to make when you create a virtual environment:
which, if any, of the packages installed in your default Python environment do
you want in your virtual environment?

Including them in your virtual environment doesn't take much extra space--they
are linked into the virtual environment instead of being copied. Within the
virtual environment you can install new versions of packages like Numpy or
Astropy that override the versions installed in your default Python environment.

The easiest way to get started is to include in your virtual environment the
packages installed in your your default Python environment; the instructions
below do that.

In everything that follows, ``ENV`` represents the name you give your virtual
environment.

**The name you choose cannot have spaces in it.**

* `conda`_:
    + Make an environment called ``ENV`` with all of the packages in your ``base``
      Miniconda environment::

        conda create --name ENV

    + More details, and examples that start with none of the packages from
      your default Python environment, are in the
      `documentation for the conda command`_ and the
      `guide on how to manage environments`_.

    .. note::
      As a general good practice, it is best to keep ``base`` untouched,
      and install any packages you may need into isolated environments.
      This way, even if you create new environments starting from ``base``,
      you control exactly which packages are installed, saving you from
      subtle dependency issues down the road.

    + Next activate the environment ``ENV`` with::

        conda activate ENV

    + Your command-line prompt will contain ``ENV`` in parentheses by default.

    + If Astropy is installed in your ``ENV`` environment, you may need to uninstall it
      in order for the development version to install properly. You can do this
      with the following command::

        conda uninstall astropy

    + Depending on your development use case, you may want to install
      additional packages into this environment in order to carry out tests,
      build documentation, extend specific additional features etc. See
      :ref:`testing-dependencies`, :ref:`builddocs`, and :ref:`Requirements for
      Astropy <astropy-main-req>` respectively to get started according to your
      use case.

* `virtualenvwrapper`_:
    + Make an environment called ``ENV`` with all of the packages in your
      default Python environment::

         mkvirtualenv --system-site-packages ENV

    + Omit the option ``--system-site-packages`` to create an environment
      without the Python packages installed in your default Python environment.
    + Environments created with `virtualenvwrapper`_ always include `pip`_
      and `setuptools <https://setuptools.readthedocs.io>`_ so that you
      can install packages within the virtual environment.
    + More details and examples are in the
      `virtualenvwrapper command documentation`_.

* `pipenv`_:
    + Make sure you are in the Astropy source directory.  See
      :ref:`get_devel` if you are unsure how to get the source code.  After
      running ``git clone <your-astropy-fork>`` run ``cd astropy/`` then::

        pipenv install --editable .

    + This both creates the virtual environment for the project
      automatically, and also installs all of Astropy's dependencies, and
      adds your Astropy repository as the version of Astropy to use in the
      environment.

    + You can activate the environment any time you're in the top-level
      ``astropy/`` directory (cloned from git) by running::

        pipenv shell

      This will open a new shell with the appropriate virtualenv enabled.

      You can also run individual commands from the virtualenv without
      activating it in the shell like::

        pipenv run python

.. _activate_env:

Activate a virtual environment
==============================

To use a new virtual environment you may need to activate it;
`virtualenvwrapper`_ will try to automatically activate your new environment
when you create it. Activation does two things (either of which you could do
manually, though it would be inconvenient):

* Puts the ``bin`` directory for the virtual environment at the front of your
  ``$PATH``.

* Adds the name of the virtual environment to your command prompt. If you
  have successfully switched to a new environment called ``ENV`` your prompt
  should look something like this: ``(ENV)[~] $``

The commands below allow you to switch between virtual environments in
addition to activating new ones.

* `conda`_: Activate the environment ``ENV`` with::

      conda activate ENV

* `virtualenvwrapper`_: Activate the environment ``ENV`` with::

      workon ENV

* `pipenv`_: Activate the environment by changing into the project
  directory (i.e. the copy of the Astropy repository on your computer) and
  running::

      pipenv shell


.. _deactivate_env:

Deactivate a virtual environment
================================

At some point you may want to go back to your default Python environment. Do
that with:

* `conda`_: ``conda deactivate``

* `virtualenvwrapper`_: ``deactivate``
    + Note that in ``virtualenvwrapper 4.1.1`` the output of
      ``mkvirtualenv`` says you should use ``source deactivate``; that does
      not seem to actually work.

* `pipenv`_: ``exit``

  .. note::

    Unlike ``virtualenv`` and ``conda``, ``pipenv`` does not manipulate
    environment variables in your current shell session.  Instead it
    launches a *subshell* which is a copy of your previous shell, in which
    it can then change some environment variables.  Therefore, any
    environment variables you change in the ``pipenv`` shell will be
    restored to their previous value (or lost entirely) when ``exit``-ing
    the subshell.

.. _delete_env:

Delete a virtual environment
============================

In both `virtualenvwrapper`_ and `conda`_ you can simply delete the
directory in which the ``ENV`` is located; both also provide commands to
make that a bit easier.  `pipenv`_ includes a command for deleting the
virtual environment associated with the current directory:

* `conda`_: ``conda remove --all --name ENV``

* `virtualenvwrapper`_: ``rmvirtualenv ENV``

* `pipenv`_: ``pipenv --rm``: As with other ``pipenv`` commands this is
  run from within the project directory.

.. _documentation for virtualenvwrapper: https://virtualenvwrapper.readthedocs.io/en/latest/install.html
.. _virtualenvwrapper command documentation: https://virtualenvwrapper.readthedocs.io/en/latest/command_ref.html
.. _conda installation instructions here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _documentation for the conda command: https://docs.conda.io/projects/conda/en/latest/commands.html
.. _guide on how to manage environments: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
