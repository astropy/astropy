************
Installation
************

.. _astropy-main-req:

Requirements
============

Astropy has the following strict requirements:

- `Python <https://www.python.org/>`_ |minimum_python_version| or later

- `Numpy`_ |minimum_numpy_version| or later

- `pytest`_ 3.1 or later

Astropy also depends on other packages for optional features:

- `scipy`_: To power a variety of features in several modules.

- `h5py <http://www.h5py.org/>`_: To read/write
  :class:`~astropy.table.Table` objects from/to HDF5 files.

- `BeautifulSoup <https://www.crummy.com/software/BeautifulSoup/>`_: To read
  :class:`~astropy.table.table.Table` objects from HTML files.

- `bleach <https://bleach.readthedocs.io/>`_: Used to sanitize text when
  disabling HTML escaping in the :class:`~astropy.table.Table` HTML writer.

- `PyYAML <http://pyyaml.org>`_: To read/write
  :class:`~astropy.table.Table` objects from/to the Enhanced CSV ASCII table
  format and to serialize mixins for various formats.

- `xmllint <http://www.xmlsoft.org/>`_: To validate VOTABLE XML files.
  This is a command line tool installed outside of Python.

- `pandas <http://pandas.pydata.org/>`_: To read/write
  :class:`~astropy.table.Table` objects from/to pandas DataFrame objects.

- `bintrees <https://pypi.python.org/pypi/bintrees>`_ for faster ``FastRBT`` and
  ``FastBST`` indexing engines with ``Table``, although these will still be
  slower in most cases than the default indexing engine.

- `sortedcontainers <https://pypi.org/project/sortedcontainers/>`_ for faster
  ``SCEngine`` indexing engine with ``Table``, although this may still be
  slower in some cases than the default indexing engine.

- `pytz <http://pythonhosted.org/pytz/>`_: To specify and convert between timezones.

- `jplephem <https://pypi.org/project/jplephem/>`_: To retrieve JPL
  ephemeris of Solar System objects.

- `matplotlib <http://matplotlib.org/>`_ 2.0 or later: To provide plotting
  functionality that `astropy.visualization` enhances.

- `scikit-image <http://scikit-image.org/>`_: To downsample a data array in `astropy.nddata.utils`.

- `setuptools <https://setuptools.readthedocs.io>`_: Used for discovery of
  entry points which are used to insert fitters into `astropy.modeling.fitting`.

- `mpmath <http://mpmath.org/>`_: Used for the 'kraft-burrows-nousek'
  interval in `~astropy.stats.poisson_conf_interval`.

- `asdf <https://github.com/spacetelescope/asdf>`_ 2.3 or later: Enables the
  serialization of various Astropy classes into a portable, hierarchical,
  human-readable representation.

- `bottleneck <https://pypi.org/project/Bottleneck/>`_: Improves the performance
  of sigma-clipping and other functionality that may required computing statistics
  on arrays with NaN values.

And the following packages can optionally be used when testing:

- `pytest-astropy <https://github.com/astropy/pytest-astropy>`_: See
  :ref:`sourcebuildtest`

- `pytest-xdist <https://docs.pytest.org/en/3.0.0/xdist.html>`_: Used for
  distributed testing.

- `pytest-mpl <https://github.com/matplotlib/pytest-mpl>`_: Used for testing
  with Matplotlib figures.

- `objgraph <https://mg.pov.lt/objgraph/>`_: Used only in tests to test for reference leaks.

- `IPython <https://ipython.org/>`__: Used for testing notebook interface of
  `~astropy.table.Table`.

- `coverage <https://coverage.readthedocs.io/>`_: Used for code coverage
  measurements.

- `skyfield <https://rhodesmill.org/skyfield/>`_: Used for testing Solar System
  coordinates.

However, note that these only need to be installed if those particular features
are needed. Astropy will import even if these dependencies are not installed.

Installing Astropy
==================

Using pip
---------

.. warning::

    Users of the Anaconda python distribution should follow the instructions
    for :ref:`anaconda_install`.

To install Astropy with `pip <https://pip.pypa.io>`__, simply run::

    pip install astropy

If you want to make sure none of your existing dependencies get upgraded, you
can also do::

    pip install astropy --no-deps

On the other hand, if you want to install Astropy along with all the available
optional depdendencies, you can do::

    pip install astropy[all]

Note that you will need a C compiler (e.g. ``gcc`` or ``clang``) to be installed
(see `Building from source`_ below) for the installation to succeed.

If you get a ``PermissionError`` this means that you do not have the required
administrative access to install new packages to your Python installation.  In
this case you may consider using the ``--user`` option to install the package
into your home directory.  You can read more about how to do this in the `pip
documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_.

Alternatively, if you intend to do development on other software that uses
Astropy, such as an affiliated package, consider installing Astropy into a
:ref:`virtualenv<using-virtualenv>`.

Do **not** install Astropy or other third-party packages using ``sudo``
unless you are fully aware of the risks.

.. _anaconda_install:

Using conda
-----------

Astropy is installed by default with the `Anaconda Distribution
<https://www.anaconda.com/download/>`_. To update to the latest version run::

    conda update astropy

There may be a delay of a day or two between when a new version of Astropy
is released and when a package is available for Anaconda. You can check
for the list of available versions with ``conda search astropy``.

.. warning::

    Attempting to use `pip <https://pip.pypa.io>`__ to upgrade your installation
    of Astropy may result in a corrupted installation.

.. _testing_installed_astropy:

Testing an installed Astropy
----------------------------

The easiest way to test your installed version of astropy is running
correctly is to use the :ref:`astropy.test()` function::

    import astropy
    astropy.test()

The tests should run and print out any failures, which you can report at
the `Astropy issue tracker <https://github.com/astropy/astropy/issues>`_.

This way of running the tests may not work if you do it in the astropy source
distribution.  See :ref:`sourcebuildtest` for how to run the tests from the
source code directory, or :ref:`running-tests` for more details.

Building from source
====================

Prerequisites
-------------

You will need a compiler suite and the development headers for Python and
Numpy in order to build Astropy.

If you are building the latest developer version rather than using a stable
release, you will also need `Cython <http://cython.org/>`_ (v0.21 or later) and
`jinja2 <http://jinja.pocoo.org/docs/dev/>`_ (v2.7 or later) installed (the
released packages have the necessary C files packaged with them, and hence do
not require Cython).

Prerequisites for Linux
-----------------------

On Linux, using the package manager for your distribution will usually be the
easiest route to making sure you have the pre-requisites to build Astropy. In
order to build from source, you'll need the python development package for your
Linux distribution.

For Debian/Ubuntu::

    sudo apt-get install python-dev

For Fedora/RHEL::

    sudo yum install python-devel

Prerequisites for Mac OS X
--------------------------

On MacOS X you will need the XCode command line tools which can be installed
using::

    xcode-select --install

Follow the onscreen instructions to install the command line tools required.
Note that you do **not** need to install the full XCode distribution (assuming
you are using MacOS X 10.9 or later).

The `instructions for building Numpy from source
<https://docs.scipy.org/doc/numpy/user/building.html>`_ are a good
resource for setting up your environment to build Python packages.

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

   git clone --recursive git://github.com/astropy/astropy.git

If you wish to participate in the development of Astropy, see
:ref:`developer-docs`.  This document covers only the basics necessary to
install Astropy.

Building and Installing
-----------------------

Astropy uses the Python built-in `distutils framework
<http://docs.python.org/install/index.html>`_ for building and
installing and requires the `setuptools`_ package.

If Numpy is not already installed in your Python environment, the astropy setup
process will try to download and install it before continuing to install
astropy.

To build and install Astropy (from the root of the source tree)::

    pip install .

If you install in this way, if you make changes to the code, you will need to
re-run the install command for changes to be reflected. Alternatively, you can
use::

    pip install -e .

Which installs astropy in develop/editable mode, which means that changes in
the code are immediately reflected in the installed version.

Troubleshooting
---------------

If you get an error mentioning that you do not have the correct permissions to
install Astropy into the default ``site-packages`` directory, you can try
installing with::

    pip install . --user

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


The C libraries currently bundled with Astropy include:

- `wcslib <http://www.atnf.csiro.au/people/mcalabre/WCS/>`_ see
  ``cextern/wcslib/README`` for the bundled version.

- `cfitsio <https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html>`_ see
  ``cextern/cfitsio/changes.txt`` for the bundled version.

- `erfa <https://github.com/liberfa>`_ see ``cextern/erfa/README.rst`` for the
  bundled version.

- `expat <http://expat.sourceforge.net/>`_ see ``cextern/expat/README`` for the
  bundled version.


Installing Astropy into CASA
----------------------------

If you want to be able to use Astropy inside `CASA
<https://casa.nrao.edu/>`_, the easiest way is to do so from inside CASA.

First, we need to make sure `pip <https://pip.pypa.io>`__ is
installed. Start up CASA as normal, and type::

    CASA <2>: from setuptools.command import easy_install

    CASA <3>: easy_install.main(['--user', 'pip'])

Now, quit CASA and re-open it, then type the following to install Astropy::

    CASA <2>: import subprocess, sys

    CASA <3>: subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'astropy'])

Then close CASA again and open it, and you should be able to import Astropy::

    CASA <2>: import astropy

Any astropy affiliated package can be installed the same way (e.g. the
`spectral-cube <http://spectral-cube.readthedocs.io/en/latest/>`_ or other
packages that may be useful for radio astronomy).

.. note:: The above instructions have not been tested on all systems.
   We know of a few examples that do work, but that is not a guarantee
   that this will work on all systems.  If you install astropy and begin to
   encounter issues with CASA, please look at the `known CASA issues
   <https://github.com/astropy/astropy/issues?q=+label%3ACASA-Installation+>`_
   and, if you don't encounter your issue there, post a new one.

.. _builddocs:

Building documentation
----------------------

.. note::

    Building the documentation is in general not necessary unless you
    are writing new documentation or do not have internet access, because
    the latest (and archive) versions of astropy's documentation should
    be available at `docs.astropy.org <http://docs.astropy.org>`_ .

Dependencies
^^^^^^^^^^^^

Building the documentation requires the Astropy source code and some additional
packages, including those in :ref:`astropy-main-req`. The easiest way to
install the extra dependencies for documentation is to install
the `sphinx-astropy <https://github.com/astropy/sphinx-astropy>`_ package,
either with pip::

    pip install sphinx-astropy

or with conda::

    conda install -c astropy sphinx-astropy

In addition to providing configuration common to packages in the Astropy
ecosystem, this package also serves as a way to automatically get the main
dependencies, including:

* `Sphinx <http://sphinx.pocoo.org>`_ - the main package we use to build
  the documentation.
* `astropy-sphinx-theme <https://github.com/astropy/astropy-sphinx-theme>`_ -
  the default 'bootstrap' theme use by Astropy and a number of affiliated
  packages.
* `sphinx-automodapi <http://sphinx-automodapi.readthedocs.io>`_ - an extension
  that makes it easy to automatically generate API documentation.
* `sphinx-gallery <https://sphinx-gallery.readthedocs.io/en/latest/>`_ - an
  extension to generate example galleries
* `numpydoc <https://numpydoc.readthedocs.io>`_ - an extension to parse
  docstrings in NumpyDoc format
* `pillow <https://pillow.readthedocs.io>`_ - used in one of the examples

In addition, if you want inheritance graphs to be generated, you will need to
make sure that `Graphviz <http://www.graphviz.org>`_ is installed. If you
install sphinx-astropy with conda, graphviz will automatically get installed,
but if you use pip, you will need to install Graphviz separately as it isn't
a Python package.

Building
^^^^^^^^

There are two ways to build the Astropy documentation. The most straightforward
way is to execute the command (from the astropy source directory)::

    python setup.py build_docs

The documentation will be built in the ``docs/_build/html`` directory, and can
be read by pointing a web browser to ``docs/_build/html/index.html``.

The LaTeX documentation can be generated by using the command::

    python setup.py build_docs -b latex

The LaTeX file ``Astropy.tex`` will be created in the ``docs/_build/latex``
directory, and can be compiled using ``pdflatex``.

The above method builds the API documentation from the source code.
Alternatively, you can do::

    cd docs
    make html

And the documentation will be generated in the same location, but using the
*installed* version of Astropy.

Reporting issues/requesting features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned above, building the documentation depends on a number of sphinx
extensions and other packages. Since it isn't necessarily straightforward to
know which package is causing issues or would need to have a new feature
implemented, you can simply open an issue in the `core astropy package issue
tracker <https://github.com/astropy/astropy/issues>`_. However, if you wish you
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

Testing a source code build of Astropy
--------------------------------------

Before running tests, it is necessary to make sure that Astropy's test
dependencies are installed. This can be done with the following command::

    pip install pytest-astropy

More information on what the ``pytest-astropy`` package provides can be found
in :ref:`testing-dependencies`.

The easiest way to test that your Astropy built correctly (without
installing astropy) is to run this from the root of the source tree::

    python setup.py test

There are also alternative methods of :ref:`running-tests`. Note that you will
need `pytest <http://pytest.org>`_ to be installed for this to work.

.. include:: development/workflow/known_projects.inc
