:orphan:

.. include:: links.inc
.. _virtual_envs:

***************************
Python virtual environments
***************************

If you plan to do regular work on astropy you should do your development in
a python virtual environment. Conceptually a virtual environment is a
duplicate of the python environment you normally work in with as many (or as
few) of the packages from your normal environment included in that virtual
environment. It is sandboxed from your normal python environment in the sense
that packages installed in the virtual environment do not affect your normal
environment in any way.

.. note::
    "Normal python environment" means whatever python you are using when you
    log in.

There are two options for using virtual environments; the choice of method is
dictated by the python distribution you use:

* If you use the anaconda python distribution you must use `conda`_ to make
  and manage your virtual environments.
* If you use any other distribution you use `virtualenvwrapper`_; you *can not*
  use `conda`_. As the name suggests, `virtualenvwrapper`_ is a wrapper around
  `virtualenv`_.

In both cases you will go through the same basic steps; the commands to
accomplish each step are given for both `conda`_ and `virtualenvwrapper`_:

* :ref:`setup_for_env`
* :ref:`list_env`
* :ref:`create_env`
* :ref:`activate_env`
* :ref:`deactivate_env`
* :ref:`delete_env`

.. note::
    + You **cannot** use `virtualenvwrapper`_ or `virtualenv`_ within anaconda.
    + `virtualenvwrapper`_ works with bash and bash-like shells; see
      :ref:`using-virtualenv` for alternatives.

.. _setup_for_env:


Set up for virtual environments
===============================

* `virtualenvwrapper`_:

  + First, install `virtualenvwrapper`_, which will also install `virtualenv`_,
    with ``pip install virtualenvwrapper``.
  + From the `documentation for virtualenvwrapper`_, you also need to::

      export WORKON_HOME=$HOME/.virtualenvs
      export PROJECT_HOME=$HOME/
      source /usr/local/bin/virtualenvwrapper.sh

* `conda`_: No setup is necessary beyond installing the anaconda python
  distribution.

.. _list_env:

List virtual environments
=========================

You do not need to list the virtual environments you have created before using
them...but sooner or later you will forget what environments you have defined
and this is the easy way to find out.

* `virtualenvwrapper`_: ``workon``
    + If this displays nothing you have no virtual environments
    + If this displays ``workon: command not found`` then you haven't done
      the :ref:`setup_for_env`; do that.
    + For more detailed information about installed environments use
      ``lsvirtualenv``.
* `conda`_: ``conda info -e``
    + you will always have at least one environment, called ``root``
    + your active environment is indicated by a ``*``

.. _create_env:

Create a new virtual environment
================================

This needs to be done once for each virtual environment you want. There is one
important choice you need to make when you create a virtual environment:
which, if any, of the packages installed in your normal python environment do
you want in your virtual environment?

Including them in your virtual environment doesn't take much extra space--they
are linked into the virtual environment instead of being copied. Within the
virtual environment you can install new versions of packages like Numpy or
Astropy that override the versions installed in your normal python environment.

The easiest way to get started is to include in your virtual environment the
packages installed in your your normal python environment; the instructions
below do that.

In everything that follows, ``ENV`` represents the name you give your virtual
environment.

**The name you choose cannot have spaces in it.**

* `virtualenvwrapper`_:
    + Make an environment called ``ENV`` with all of the packages in your normal
      python environment::

         mkvirtualenv --system-site-packages ENV

    + Omit the option ``--system-site-packages`` to create an environment
      without the python packages installed in your normal python environment.
    + Environments created with `virtualenvwrapper`_ always include `pip`_
      and `setuptools <https://setuptools.readthedocs.io>`_ so that you
      can install packages within the virtual environment.
    + More details and examples are in the
      `virtualenvwrapper command documentation`_.
* `conda`_:
    + Make an environment called ``ENV`` with all of the packages in your main
      anaconda environment::

        conda create -n ENV anaconda

    + More details, and examples that start with none of the packages from
      your normal python environment, are in the
      `documentation for the conda command`_ and the
      `guide on how to manage environments`_.

    + If astropy is installed in your environment, you may need to uninstall it
      in order for the development version to install properly. You can do this
      with the following command::

        conda uninstall astropy

.. _activate_env:

Activate a virtual environment
==============================

To use a new virtual environment you may need to activate it;
`virtualenvwrapper`_ will try to automatically activate your new environment
when you create it. Activation does two things (either of which you could do
manually, though it would be inconvenient):

* Put the ``bin`` directory for the virtual environment at the front of your
  ``$PATH``.
* Add the name of the virtual environment to your command prompt. If you have
  successfully switched to a new environment called ``ENV`` your prompt should
  look something like this: ``(ENV)[~] $``

The commands below allow you to switch between virtual environments in
addition to activating new ones.

* `virtualenvwrapper`_: Activate the environment ``ENV`` with::

      workon ENV

* ` conda`: Activate the environment ``ENV`` with::

      conda activate ENV


.. _deactivate_env:

Deactivate a virtual environment
================================

At some point you may want to go back to your normal python environment. Do
that with:

* `virtualenvwrapper`_: ``deactivate``
    + Note that in ``virtualenvwrapper 4.1.1`` the output of
      ``mkvirtualenv`` says you should use ``source deactivate``; that does
      not seem to actually work.
* `conda`_: ``conda deactivate``

.. _delete_env:

Delete a virtual environment
============================

In both `virtualenvwrapper`_ and `conda`_ you can simply delete the directory in
which the ``ENV`` is located; both also provide commands to make that a bit easier.

* `virtualenvwrapper`_: ``rmvirtualenv ENV``
* `conda`_: ``conda remove --all -n ENV``

.. _documentation for virtualenvwrapper: https://virtualenvwrapper.readthedocs.io/en/latest/install.html
.. _virtualenvwrapper command documentation: https://virtualenvwrapper.readthedocs.io/en/latest/command_ref.html
.. _documentation for the conda command: https://docs.conda.io/projects/conda/en/latest/commands.html
.. _guide on how to manage environments: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
