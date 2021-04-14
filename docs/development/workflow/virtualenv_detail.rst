:orphan:

.. _using-virtualenv:

Using virtualenv
****************

`virtualenv`_ is a tool for creating and activating isolated Python
environments that allow installing and experimenting with Python packages
without disrupting your production Python environment.  When using commands
such as ``pip install -e .``, for example, it is strongly recommended to
do so within a virtualenv.  This is generally preferable to installing a
development version of Astropy into your system site-packages and having to
keep track of whether or not your environment is in a "known good"
configuration for production/science use.

Using a virtualenv is also a good way to try out new versions of software that
you're not actively doing development work on without disrupting your normal
production environment.

We won't provide a full tutorial on using virtualenv here |emdash| the
virtualenv documentation linked to above is a better place to start.  But here
is a quick overview on how to set up a virtualenv for Astropy development with
your default Python version:

#. Install virtualenv::

       $ pip install virtualenv

   or (on Debian/Ubuntu)::

       $ sudo apt-get install python-virtualenv

   etc.

#. (Recommended) Create a root directory for all your virtualenvs under a path
   you have write access to.  For example::

       $ mkdir ~/.virtualenvs

#. Create the Astropy virtualenv::

       $ virtualenv --distribute --system-site-packages ~/.virtualenvs/astropy-dev

   The ``--system-site-packages`` option inherits all packages already
   installed in your system site-packages directory; this frees you from having
   to reinstall packages like Numpy and Scipy in the virtualenv.  However, if
   you would like your virtualenv to use a development version of Numpy, for
   example, you can still install Numpy into the virtualenv and it will take
   precedence over the version installed in site-packages.

#. Activate the virtualenv::

       $ source ~/.virtualenvs/astropy-dev/bin/activate

   or if you're using a csh-variant::

       $ source ~/.virtualenvs/astropy-dev/bin/activate.csh

   virtualenv works on Windows too |emdash| see the documentation for details.

#. If the virtualenv successfully activated its name should appear in your
   shell prompt::

       (astropy-dev) $

   The virtualenv can be disabled at any time by entering::

       (astropy-dev) $ deactivate

#. Now as long as the virtualenv is activated, packages you install with
   ``pip`` will automatically install into your virtualenv instead of the system
   site-packages.  Consider installing Astropy in develop mode into the
   virtualenv as described :ref:`activate_development_astropy`.

virtualenvwrapper
=================

`virtualenvwrapper`_ is a set of enhancements to virtualenv mostly
implemented through simple shell scripts and aliases.  It automatically
organizes all your virtualenvs under a single directory (as suggested
above). To create a new virtualenv you can just use the ``'mkvirtualenv
<env_name>'`` command and it will automatically create a new virtualenv of
that name in the default location.

To activate a virtualenv with virtualenvwrapper you don't need to think
about the environment's location of the filesystem or which activate script
to run.  Simply run ``'workon <env_name>'``.  You can also list all
virtualenvs with ``lsvirtualenv``.  That just scratches the surface of the
goodies included with virtualenvwrapper.

The one caveat is that it does not support csh-like shells.
There exists `virtualenvwrapper-win`_, which ports virtualenvwrapper to
Windows batch scripts.

venv
====

virtualenv is so commonly used in the Python development community that its
functionality was finally added to the standard library in Python 3.3 under
the name `venv`_.  venv has not gained wide use yet and is not explicitly
supported by tools like virtualenvwrapper, but it is expected to see wider
adoption in the future.

.. include:: links.inc
