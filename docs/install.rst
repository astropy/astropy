************
Installation
************

Requirements
============

Astropy has the following strict requirements:

- `Python <http://www.python.org/>`_ 2.6 (>=2.6.5), 2.7, 3.3, or 3.4

  - Prior to Astropy v1.0 Python 3.1 and 3.2 are also supported.

- `Numpy`_ |minimum_numpy_version| or later

Astropy also depends on other packages for optional features:

- `h5py <http://h5py.org/>`_: To read/write
  :class:`~astropy.table.Table` objects from/to HDF5 files

- `BeautifulSoup <http://www.crummy.com/software/BeautifulSoup/>`_: To read
  :class:`~astropy.table.table.Table` objects from HTML files

- `PyYAML <http://pyyaml.org>`_: To read/write
  :class:`~astropy.table.Table` objects from/to the Enhanced CSV ASCII table format.

- `scipy`_: To power a variety of features (currently
  mainly cosmology-related functionality)

- `xmllint <http://www.xmlsoft.org/>`_: To validate VOTABLE XML files.

- `matplotlib <http://matplotlib.org/>`_: To provide plotting functionality that `astropy.visualization` enhances.

- `WCSAxes <http://wcsaxes.readthedocs.org/en/latest/>`_: To use `astropy.wcs` to define projections in Matplotlib.

However, note that these only need to be installed if those particular features
are needed. Astropy will import even if these dependencies are not installed.

.. TODO: Link to the planned dependency checker/installer tool.

Installing Astropy
==================

Using pip
-------------

To install Astropy with `pip <http://www.pip-installer.org/en/latest/>`_, simply run::

    pip install --no-deps astropy

.. warning::

    Users of the Anaconda python distribution should follow the instructions
    for :ref:`anaconda_install`.

.. note::

    You will need a C compiler (e.g. ``gcc`` or ``clang``) to be installed (see
    `Building from source`_ below) for the installation to succeed.

.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to "help" you
    by upgrading your Numpy installation, which may not always be desired.

.. note::

    If you get a ``PermissionError`` this means that you do not have the
    required administrative access to install new packages to your Python
    installation.  In this case you may consider using the ``--user`` option
    to install the package into your home directory.  You can read more
    about how to do this in the `pip documentation
    <http://www.pip-installer.org/en/1.2.1/other-tools.html#using-pip-with-the-user-scheme>`_.

    Alternatively, if you intend to do development on other software that uses
    Astropy, such as an affiliated package, consider installing Astropy into a
    :ref:`virtualenv<using-virtualenv>`.

    Do **not** install Astropy or other third-party packages using ``sudo``
    unless you are fully aware of the risks.


.. _anaconda_install:

Anaconda python distribution
----------------------------

Astropy is installed by default with Anaconda. To update to the latest version
run::

    conda update astropy

.. note::

    There may be a delay of a day or two between when a new version of Astropy
    is released and when a package is available for Anaconda. You can check
    for the list of available versions with ``conda search astropy``.

.. note::

    Attempting to use ``pip`` to upgrade your installation of Astropy may result
    in a corrupted installation.


Binary installers
-----------------

Binary installers are available on Windows for Python 2.6, 2.7, and >= 3.3
at `PyPI <https://pypi.python.org/pypi/astropy>`_.

.. _testing_installed_astropy:


Testing an installed Astropy
----------------------------

The easiest way to test your installed version of astropy is running
correctly is to use the :ref:`astropy.test()` function::

    import astropy
    astropy.test()

The tests should run and print out any failures, which you can report at
the `Astropy issue tracker <http://github.com/astropy/astropy/issues>`_.

.. note::

    This way of running the tests may not work if you do it in the
    astropy source distribution.  See :ref:`sourcebuildtest` for how to
    run the tests from the source code directory, or :ref:`running-tests`
    for more details.

.. note::

    Running the tests this way is currently disabled in the IPython REPL due
    to conflicts with some common display settings in IPython.  Please run the
    Astropy tests under the standard Python command-line interpreter.



Building from source
====================

Prerequisites
-------------

You will need a compiler suite and the development headers for Python and
Numpy in order to build Astropy. On Linux, using the package manager for your
distribution will usually be the easiest route, while on MacOS X you will
need the XCode command line tools.

The `instructions for building Numpy from source
<http://docs.scipy.org/doc/numpy/user/install.html>`_ are also a good
resource for setting up your environment to build Python packages.

You will also need `Cython <http://cython.org/>`_ (v0.15 or later) and
`jinja2 <http://jinja.pocoo.org/docs/dev/>`_ (v2.7 or later) installed
to build from source, unless you are installing a numbered release. (The
releases packages have the necessary C files packaged with them, and hence do
not require Cython.)

.. note::

    If you are using MacOS X, you will need to the XCode command line tools.
    One way to get them is to install `XCode
    <https://developer.apple.com/xcode/>`_. If you are using OS X 10.7 (Lion)
    or later, you must also explicitly install the command line tools. You can
    do this by opening the XCode application, going to **Preferences**, then
    **Downloads**, and then under **Components**, click on the Install button
    to the right of **Command Line Tools**.  Alternatively, on 10.7 (Lion) or
    later, you do not need to install XCode, you can download just the command
    line tools from https://developer.apple.com/downloads/index.action
    (requires an Apple developer account).


Obtaining the source packages
-----------------------------

Source packages
^^^^^^^^^^^^^^^

The latest stable source package for Astropy can be `downloaded here
<https://pypi.python.org/pypi/astropy>`_.

Development repository
^^^^^^^^^^^^^^^^^^^^^^

The latest development version of Astropy can be cloned from github
using this command::

   git clone git://github.com/astropy/astropy.git

.. note::

   If you wish to participate in the development of Astropy, see
   :ref:`developer-docs`.  This document covers only the basics
   necessary to install Astropy.

Building and Installing
-----------------------

Astropy uses the Python `distutils framework
<http://docs.python.org/install/index.html>`_ for building and
installing and requires the
`distribute <http://pypi.python.org/pypi/distribute>`_ extension--the later is
automatically downloaded when running ``python setup.py`` if it is not already
provided by your system.

If Numpy is not already installed in your Python environment, the
astropy setup process will try to download and install it before
continuing to install astropy.

To build Astropy (from the root of the source tree)::

    python setup.py build

To install Astropy (from the root of the source tree)::

    python setup.py install


Troubleshooting
---------------

If you get an error mentioning that you do not have the correct permissions to
install Astropy into the default ``site-packages`` directory, you can try
installing with::

    python setup.py install --user

which will install into a default directory in your home directory.


External C libraries
^^^^^^^^^^^^^^^^^^^^

The Astropy source ships with the C source code of a number of
libraries.  By default, these internal copies are used to build
Astropy.  However, if you wish to use the system-wide installation of
one of those libraries, you can pass one or more of the
``--use-system-X`` flags to the ``setup.py build`` command.

For example, to build Astropy using the system `libexpat
<http://www.libexpat.org/>`_, use::

    python setup.py build --use-system-expat

To build using all of the system libraries, use::

    python setup.py build --use-system-libraries

To see which system libraries Astropy knows how to build against, use::

    python setup.py build --help

As with all distutils commandline options, they may also be provided in a
``setup.cfg`` in the same directory as ``setup.py``.  For example, to use
the system `libexpat <http://www.libexpat.org/>`_, add the following to the
``setup.cfg`` file::

    [build]
    use_system_expat=1


The required version of setuptools is not available
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If upon running the ``setup.py`` script you get a message like

    The required version of setuptools (>=0.9.8) is not available,
    and can't be installed while this script is running. Please
    install a more recent version first, using
    'easy_install -U setuptools'.

    (Currently using setuptools 0.6c11 (/path/to/setuptools-0.6c11-py2.7.egg))

this is because you have a very outdated version of the `setuptools
<https://pythonhosted.org/setuptools/>`_ package which is used to install
Python packages.  Normally Astropy will bootstrap newer version of
setuptools via the network, but setuptools suggests that you first
*uninstall* the old version (the ``easy_install -U setuptools`` command).

However, in the likely case that your version of setuptools was installed by an
OS system package (on Linux check your package manager like apt or yum for a
package called ``python-setuptools``), trying to uninstall with
``easy_install`` and without using ``sudo`` may not work, or may leave your
system package in an inconsistent state.

As the best course of action at this point depends largely on the individual
system and how it is configured, if you are not sure yourself what do please
ask on the Astropy mailing list.


The Windows installer can't find Python in the registry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a common issue with Windows installers for Python packages that do not
support the new User Access Control (UAC) framework added in Windows Vista and
later.  In particular, when a Python is installed "for all users" (as opposed
to for a single user) it adds entries for that Python installation under the
``HKEY_LOCAL_MACHINE`` (HKLM) hierarchy and *not* under the
``HKEY_CURRENT_USER`` (HKCU) hierarchy.  However, depending on your UAC
settings, if the Astropy installer is not executed with elevated privileges it
will not be able to check in HKLM for the required information about your
Python installation.

In short: If you encounter this problem it's because you need the appropriate
entries in the Windows registry for Python. You can download `this script`__
and execute it with the same Python as the one you want to install Astropy
into.  For example to add the missing registry entries to your Python 2.7::

    C:\>C:\Python27\python.exe C:\Path\To\Downloads\win_register_python.py

__ https://gist.github.com/embray/6042780#file-win_register_python-py

.. _builddocs:

Building documentation
----------------------

.. note::

    Building the documentation is in general not necessary unless you
    are writing new documentation or do not have internet access, because
    the latest (and archive) versions of astropy's documentation should
    be available at `docs.astropy.org <http://docs.astropy.org>`_ .

Building the documentation requires the Astropy source code and some additional
packages:

    - `Sphinx <http://sphinx.pocoo.org>`_ (and its dependencies) 1.0 or later

    - `Graphviz <http://www.graphviz.org>`_

    - `Astropy-helpers <https://github.com/astropy/astropy-helpers>`_ (Astropy
      and most affiliated packages include this as a submodule in the source
      repository, so it does not need to be installed separately.)

    - `WCSAxes <http://wcsaxes.readthedocs.org/en/latest/>`_

.. note::

    Sphinx also requires a reasonably modern LaTeX installation to render
    equations.  Per the `Sphinx documentation
    <http://sphinx-doc.org/builders.html?highlight=latex#sphinx.builders.latex.LaTeXBuilder>`_,
    for the TexLive distribution the following packages are required to be
    installed:

    * latex-recommended
    * latex-extra
    * fonts-recommended

    For other LaTeX distributions your mileage may vary. To build the PDF
    documentation using LaTeX, the ``fonts-extra`` TexLive package or the
    ``inconsolata`` CTAN package are also required.

There are two ways to build the Astropy documentation. The most straightforward
way is to execute the command (from the astropy source directory)::

    python setup.py build_docs

The documentation will be built in the ``build/sphinx/html`` directory, and can
be read by pointing a web browser to ``build/sphinx/html/index.html``.

The LaTeX documentation can be generated by using the command::

    python setup.py build_docs -b latex

The LaTeX file ``Astropy.tex`` will be created in the ``build/sphinx/latex``
directory, and can be compiled using ``pdflatex``.

The above method builds the API documentation from the source code.
Alternatively, you can do::

    cd docs
    make html

And the documentation will be generated in the ``docs/_build/html/``
directory, but using the *installed* version of Astropy.

.. _sourcebuildtest:

Testing a source code build of Astropy
--------------------------------------

The easiest way to test that your Astropy built correctly (without
installing astropy) is to run this from the root of the source tree::

    python setup.py test

There are also alternative methods of :ref:`running-tests`.

.. include:: development/workflow/known_projects.inc
