************
Installation
************

.. _astropy-main-req:

Requirements
============

``astropy`` has the following strict requirements:

- `Python <https://www.python.org/>`_ |minimum_python_version| or later

- `Numpy`_ |minimum_numpy_version| or later

- `pytest`_ 3.1 or later

``astropy`` also depends on other packages for optional features:

- `scipy`_ |minimum_scipy_version| or later: To power a variety of features
  in several modules.

- `h5py <http://www.h5py.org/>`_: To read/write
  :class:`~astropy.table.Table` objects from/to HDF5 files.

- `BeautifulSoup <https://www.crummy.com/software/BeautifulSoup/>`_: To read
  :class:`~astropy.table.table.Table` objects from HTML files.

- `html5lib <https://html5lib.readthedocs.io/en/stable/>`_: To read
  :class:`~astropy.table.table.Table` objects from HTML files using the
  `pandas <https://pandas.pydata.org/>`_ reader.

- `bleach <https://bleach.readthedocs.io/>`_: Used to sanitize text when
  disabling HTML escaping in the :class:`~astropy.table.Table` HTML writer.

- `PyYAML <https://pyyaml.org>`_ |minimum_yaml_version| or later: To read/write
  :class:`~astropy.table.Table` objects from/to the Enhanced CSV ASCII table
  format and to serialize mixins for various formats.

- `xmllint <http://www.xmlsoft.org/>`_: To validate VOTABLE XML files.
  This is a command line tool installed outside of Python.

- `pandas <https://pandas.pydata.org/>`_: To convert
  :class:`~astropy.table.Table` objects from/to pandas DataFrame objects.
  Version 0.14 or higher is required to use the :ref:`table_io_pandas`
  I/O functions to read/write :class:`~astropy.table.Table` objects.

- `sortedcontainers <https://pypi.org/project/sortedcontainers/>`_ for faster
  ``SCEngine`` indexing engine with ``Table``, although this may still be
  slower in some cases than the default indexing engine.

- `bintrees <https://pypi.org/project/bintrees>`_ for faster ``FastRBT`` and
  ``FastBST`` indexing engines with ``Table``, although these will still be
  slower in most cases than the default indexing engine. *This package is
  deprecated because it is no longer maintained.  The ``sortedcontainers``
  package is now recommended.*

- `pytz <https://pythonhosted.org/pytz/>`_: To specify and convert between
  timezones.

- `jplephem <https://pypi.org/project/jplephem/>`_: To retrieve JPL
  ephemeris of Solar System objects.

- `matplotlib <https://matplotlib.org/>`_ 2.0 or later: To provide plotting
  functionality that `astropy.visualization` enhances.

- `scikit-image <https://scikit-image.org/>`_: To downsample a data array in
  `astropy.nddata.utils`.

- `setuptools <https://setuptools.readthedocs.io>`_: Used for discovery of
  entry points which are used to insert fitters into `astropy.modeling.fitting`.

- `mpmath <http://mpmath.org/>`_: Used for the 'kraft-burrows-nousek'
  interval in `~astropy.stats.poisson_conf_interval`.

- `asdf <https://github.com/spacetelescope/asdf>`_ 2.3 or later: Enables the
  serialization of various Astropy classes into a portable, hierarchical,
  human-readable representation.

- `bottleneck <https://pypi.org/project/Bottleneck/>`_: Improves the performance
  of sigma-clipping and other functionality that may require computing
  statistics on arrays with NaN values.

And the following packages can optionally be used when testing:

- `pytest-astropy <https://github.com/astropy/pytest-astropy>`_: See
  :ref:`sourcebuildtest`

- `pytest-xdist <https://pypi.org/project/pytest-xdist/>`_: Used for
  distributed testing.

- `pytest-mpl <https://github.com/matplotlib/pytest-mpl>`_: Used for testing
  with Matplotlib figures.

- `objgraph <https://mg.pov.lt/objgraph/>`_: Used only in tests to test for reference leaks.

- `IPython <https://ipython.org/>`__: Used for testing the notebook interface of
  `~astropy.table.Table`.

- `coverage <https://coverage.readthedocs.io/>`_: Used for code coverage
  measurements.

- `skyfield <https://rhodesmill.org/skyfield/>`_: Used for testing Solar System
  coordinates.

However, note that these packages require installation only if those particular
features are needed. ``astropy`` will import even if these dependencies are not
installed.

Installing ``astropy``
======================

If you are new to Python and/or do not have familiarity with `Python virtual
environments <https://docs.python.org/3/tutorial/venv.html>`_, then we recommend
starting by installing the `Anaconda Distribution
<https://www.anaconda.com/distribution/>`_. This works on all platforms (linux,
Mac, Windows) and installs a full-featured scientific Python in a user directory
without requiring root permissions.

Using pip
---------

.. warning::

    Users of the Anaconda Python distribution should follow the instructions
    for :ref:`anaconda_install`.

To install ``astropy`` with `pip <https://pip.pypa.io>`__, run::

    pip install astropy

If you want to make sure none of your existing dependencies get upgraded, you
can also do::

    pip install astropy --no-deps

On the other hand, if you want to install ``astropy`` along with all of the
available optional dependencies, you can do::

    pip install astropy[all]

Note that you will need a C compiler (e.g. ``gcc`` or ``clang``) to be installed
(see `Building from source`_ below) for the installation to succeed.

If you get a ``PermissionError`` this means that you do not have the required
administrative access to install new packages to your Python installation. In
this case you may consider using the ``--user`` option to install the package
into your home directory. You can read more about how to do this in the `pip
documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_.

Alternatively, if you intend to do development on other software that uses
``astropy``, such as an affiliated package, consider installing ``astropy`` into a
:ref:`virtualenv<using-virtualenv>`.

Do **not** install ``astropy`` or other third-party packages using ``sudo``
unless you are fully aware of the risks.

.. _anaconda_install:

Using Conda
-----------

To install ``astropy`` using conda run::

    conda install astropy

``astropy`` is installed by default with the `Anaconda Distribution
<https://www.anaconda.com/distribution/>`_. To update to the latest version run::

    conda update astropy

There may be a delay of a day or two between when a new version of ``astropy``
is released and when a package is available for conda. You can check
for the list of available versions with ``conda search astropy``.

It is highly recommended that you install all of the optional dependencies with::

    conda install -c astropy -c defaults \
      scipy h5py beautifulsoup4 html5lib bleach pyyaml pandas sortedcontainers \
      pytz matplotlib setuptools mpmath bottleneck jplephem asdf

To also be able to run tests (see below) and support :ref:`builddocs` use the
following. We use ``pip`` for these packages to ensure getting the latest
releases which are compatible with the latest ``pytest`` and ``sphinx`` releases::

    pip install pytest-astropy sphinx-astropy

.. warning::

    Attempting to use `pip <https://pip.pypa.io>`__ to upgrade your installation
    of ``astropy`` itself may result in a corrupted installation.

.. _testing_installed_astropy:

Testing an Installed ``astropy``
--------------------------------

The easiest way to test if your installed version of ``astropy`` is running
correctly is to use the :ref:`astropy.test()` function::

    import astropy
    astropy.test()

The tests should run and print out any failures, which you can report at
the `Astropy issue tracker <https://github.com/astropy/astropy/issues>`_.

This way of running the tests may not work if you do it in the ``astropy`` source
distribution. See :ref:`sourcebuildtest` for how to run the tests from the
source code directory, or :ref:`running-tests` for more details.

Building from Source
====================

Prerequisites
-------------

You will need a compiler suite and the development headers for Python and
NumPy in order to build ``astropy``.

If you are building the latest developer version rather than using a stable
release, you will also need `Cython <https://cython.org/>`_ (v0.29.13 or later) and
`jinja2 <https://jinja.palletsprojects.com/en/master/>`_ (v2.7 or later) installed. The
released packages have the necessary C files packaged with them, and hence do
not require Cython.

Prerequisites for Linux
-----------------------

On Linux, using the package manager for your distribution will usually be the
easiest route to making sure you have the prerequisites to build ``astropy``. In
order to build from source, you will need the Python and Numpy development
package for your Linux distribution, as well as a number of helper packages.

For Debian/Ubuntu::

    sudo apt-get install python3-dev python3-numpy-dev python3-setuptools cython3 python3-jinja2 python3-pytest-astropy

For Fedora/RHEL::

    sudo yum install python3-devel python3-numpy python3-setuptools python3-Cython python3-jinja2 python3-pytest-astropy

.. note:: Building the developer version of ``astropy`` may require
          newer versions of the above packages than are available in
          your distribution's repository.  If so, you could either try
          a more up-to-date distribution (such as Debian ``testing``),
          or install more up-to-date versions of the packages using
          ``pip`` or ``conda`` in a virtual environment.

Prerequisites for Mac OS X
--------------------------

On MacOS X you will need the XCode command line tools which can be installed
using::

    xcode-select --install

Follow the onscreen instructions to install the command line tools required.
Note that you do **not** need to install the full XCode distribution (assuming
you are using MacOS X 10.9 or later).

The `instructions for building NumPy from source
<https://numpy.org/doc/stable/user/building.html>`_ are a good
resource for setting up your environment to build Python packages.

Obtaining the Source Packages
-----------------------------

Source Packages
^^^^^^^^^^^^^^^

The latest stable source package for ``astropy`` can be `downloaded here
<https://pypi.org/project/astropy>`_.

Development Repository
^^^^^^^^^^^^^^^^^^^^^^

The latest development version of ``astropy`` can be cloned from GitHub
using this command::

   git clone --recursive git://github.com/astropy/astropy.git

If you wish to participate in the development of ``astropy``, see
:ref:`developer-docs`. This document covers only the basics necessary to
installing ``astropy``.

Building and Installing
-----------------------

``astropy`` uses the Python built-in `distutils framework
<https://docs.python.org/install/index.html>`_ for building and
installing, and requires the `setuptools`_ package.

If NumPy is not already installed in your Python environment, the ``astropy``
setup process will try to download and install it before continuing to install
``astropy``.

To build and install ``astropy`` (from the root of the source tree)::

    pip install .

If you install in this way and you make changes to the code, you will need to
re-run the install command for changes to be reflected. Alternatively, you can
use::

    pip install -e .

which installs ``astropy`` in develop/editable mode, which means that changes in
the code are immediately reflected in the installed version.

Troubleshooting
---------------

If you get an error mentioning that you do not have the correct permissions to
install ``astropy`` into the default ``site-packages`` directory, you can try
installing with::

    pip install . --user

which will install into a default directory in your home directory.


External C Libraries
^^^^^^^^^^^^^^^^^^^^

The ``astropy`` source ships with the C source code of a number of
libraries. By default, these internal copies are used to build
``astropy``. However, if you wish to use the system-wide installation of
one of those libraries, you can pass one or more of the
``--use-system-X`` flags to the ``setup.py build`` command.

For example, to build ``astropy`` using the system's expat parser
library, use::

    python setup.py build --use-system-expat

To build using all of the system libraries, use::

    python setup.py build --use-system-libraries

To see which system libraries ``astropy`` knows how to build against, use::

    python setup.py build --help

As with all distutils command line options, they may also be provided in a
``setup.cfg`` in the same directory as ``setup.py``. For example, to use
the system `libexpat <http://www.libexpat.org/>`_, add the following to the
``setup.cfg`` file::

    [build]
    use_system_expat=1


The C libraries currently bundled with ``astropy`` include:

- `wcslib <https://www.atnf.csiro.au/people/mcalabre/WCS/>`_ see
  ``cextern/wcslib/README`` for the bundled version.

- `cfitsio <https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html>`_ see
  ``cextern/cfitsio/changes.txt`` for the bundled version.

- `erfa <https://github.com/liberfa>`_ see ``cextern/erfa/README.rst`` for the
  bundled version.

- `expat <https://expat.sourceforge.io/>`_ see ``cextern/expat/README`` for the
  bundled version.


Installing ``astropy`` into CASA
--------------------------------

If you want to be able to use ``astropy`` inside `CASA
<https://casa.nrao.edu/>`_, the easiest way is to do so from inside CASA.

First, we need to make sure `pip <https://pip.pypa.io>`__ is
installed. Start up CASA as normal, and then type::

    CASA <2>: from setuptools.command import easy_install

    CASA <3>: easy_install.main(['--user', 'pip'])

Now, quit CASA and re-open it, then type the following to install ``astropy``::

    CASA <2>: import subprocess, sys

    CASA <3>: subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'astropy'])

Then close CASA again and open it, and you should be able to import ``astropy``::

    CASA <2>: import astropy

Any ``astropy`` affiliated package can be installed the same way (e.g. the
`spectral-cube <https://spectral-cube.readthedocs.io>`_ or other
packages that may be useful for radio astronomy).

.. note:: The above instructions have not been tested on all systems.
   We know of a few examples that do work, but that is not a guarantee
   that this will work on all systems. If you install ``astropy`` and begin to
   encounter issues with CASA, please look at the `known CASA issues
   <https://github.com/astropy/astropy/issues?q=+label%3ACASA-Installation+>`_
   and if you do not encounter your issue there, please post a new one.

.. _builddocs:

Building Documentation
----------------------

.. note::

    Building the documentation is in general not necessary unless you are
    writing new documentation or do not have internet access, because
    the latest (and archive) versions of Astropy's documentation should
    be available at `docs.astropy.org <https://docs.astropy.org>`_ .

Dependencies
^^^^^^^^^^^^

Building the documentation requires the ``astropy`` source code and some
additional packages, including those in :ref:`astropy-main-req`. The easiest
way to install the extra dependencies for documentation is to install
the `sphinx-astropy <https://github.com/astropy/sphinx-astropy>`_ package,
either with pip::

    pip install sphinx-astropy

or with Conda::

    conda install -c astropy sphinx-astropy

In addition to providing configuration common to packages in the Astropy
ecosystem, this package also serves as a way to automatically get the main
dependencies, including:

* `Sphinx <http://www.sphinx-doc.org/>`_ - the main package we use to build
  the documentation
* `astropy-sphinx-theme <https://github.com/astropy/astropy-sphinx-theme>`_ -
  the default 'bootstrap' theme used by ``astropy`` and a number of affiliated
  packages
* `sphinx-automodapi <https://sphinx-automodapi.readthedocs.io>`_ - an extension
  that makes it easy to automatically generate API documentation
* `sphinx-gallery <https://sphinx-gallery.readthedocs.io/en/latest/>`_ - an
  extension to generate example galleries
* `numpydoc <https://numpydoc.readthedocs.io>`_ - an extension to parse
  docstrings in NumPyDoc format
* `pillow <https://pillow.readthedocs.io>`_ - used in one of the examples

In addition, if you want inheritance graphs to be generated, you will need to
make sure that `Graphviz <http://www.graphviz.org>`_ is installed. If you
install sphinx-astropy with Conda, Graphviz will automatically get installed,
but if you use pip, you will need to install Graphviz separately as it is not
a Python package.

.. _astropy-doc-building:

Building
^^^^^^^^

There are two ways to build the Astropy documentation. The first way is to
execute the command (from the ``astropy`` source directory)::

    python setup.py build_docs

The documentation will be built in the ``docs/_build/html`` directory, and can
be read by pointing a web browser to ``docs/_build/html/index.html``.

In the second way, LaTeX documentation can be generated by using the command::

    python setup.py build_docs -b latex

The LaTeX file ``Astropy.tex`` will be created in the ``docs/_build/latex``
directory, and can be compiled using ``pdflatex``.

The above method builds the API documentation from the source code.
Alternatively, you can do::

    cd docs
    make html

And the documentation will be generated in the same location, but using the
*installed* version of ``astropy``.

Reporting Issues/Requesting Features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned above, building the documentation depends on a number of sphinx
extensions and other packages. Since it is not always possible to know which
package is causing issues or would need to have a new feature implemented, you
can open an issue in the `core astropy package issue
tracker <https://github.com/astropy/astropy/issues>`_. However, if you wish, you
can also open issues in the repositories for some of the dependencies:

* For requests/issues related to the appearance of the docs (e.g. related to
  the CSS), you can open an issue in the `astropy-sphinx-theme issue tracker
  <https://github.com/astropy/astropy-sphinx-theme/issues>`_.

* For requests/issues related to the auto-generated API docs which appear to
  be general issues rather than an issue with a specific docstring, you can use
  the `sphinx-automodapi issue tracker
  <https://github.com/astropy/sphinx-automodapi/issues>`_.

* For issues related to the default configuration (e.g which extensions are
  enabled by default), you can use the `sphinx-astropy issue tracker
  <https://github.com/astropy/sphinx-astropy/issues>`_.

.. _sourcebuildtest:

Testing a Source Code Build of ``astropy``
------------------------------------------

Before running tests, it is necessary to make sure that Astropy's test
dependencies are installed. This can be done with the following command::

    pip install pytest-astropy

More information on what the ``pytest-astropy`` package provides can be found
in :ref:`testing-dependencies`.

The most convenient way to test that your Astropy built correctly
(without installing ``astropy``) is to run this from the root of the source tree::

    python setup.py test

There are also alternative methods of :ref:`running-tests`. Note that you will
need `pytest`_ to be installed for this to work.

.. include:: development/workflow/known_projects.inc
