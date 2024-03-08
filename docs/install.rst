************
Installation
************

.. _installing-astropy:

Installing ``astropy``
**********************

If you are new to Python and/or do not have familiarity with `Python virtual
environments <https://docs.python.org/3/tutorial/venv.html>`_, then we recommend
starting by installing the `Anaconda Distribution
<https://www.anaconda.com/download/>`_. This works on all platforms (linux,
Mac, Windows) and installs a full-featured scientific Python in a user directory
without requiring root permissions.

Using pip
=========

.. warning::

    Users of the Anaconda Python distribution should follow the instructions
    for :ref:`anaconda_install`.

To install ``astropy`` with |pip|, run::

    python -m pip install astropy

If you want to make sure none of your existing dependencies get upgraded, you
can also do::

    python -m pip install astropy --no-deps

On the other hand, if you want to install ``astropy`` along with recommended
or even all of the available optional :ref:`dependencies <astropy-main-req>`,
you can do::

    python -m pip install "astropy[recommended]"

or::

    python -m pip install "astropy[all]"

In most cases, this will install a pre-compiled version (called a *wheel*) of
astropy, but if you are using a very recent version of Python, if a new version
of astropy has just been released, or if you are building astropy for a platform
that is not common, astropy will be installed from a source file. Note that in
this case you will need a C compiler (e.g., ``gcc`` or ``clang``) to be installed
(see `Building from source`_ below) for the installation to succeed.

If you get a ``PermissionError`` this means that you do not have the required
administrative access to install new packages to your Python installation. In
this case you may consider using the ``--user`` option to install the package
into your home directory. You can read more about how to do this in the `pip
documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_.

Alternatively, if you intend to do development on other software that uses
``astropy``, such as an affiliated package, consider installing ``astropy``
into a :ref:`virtualenv <astropy-dev:virtual_envs>`.

Do **not** install ``astropy`` or other third-party packages using ``sudo``
unless you are fully aware of the risks.

.. _anaconda_install:

Using Conda
===========

To install ``astropy`` using conda run::

    conda install astropy

``astropy`` is installed by default with the `Anaconda Distribution
<https://www.anaconda.com/download/>`_. To update to the latest version run::

    conda update astropy

There may be a delay of a day or two between when a new version of ``astropy``
is released and when a package is available for conda. You can check
for the list of available versions with ``conda search astropy``.

If you want to install ``astropy`` along with recommended or all of the
available optional :ref:`dependencies <astropy-main-req>`, you can do::

    conda install --channel conda-forge --channel defaults scipy matplotlib

or::

    conda install --channel conda-forge --channel defaults scipy matplotlib \
      h5py beautifulsoup4 html5lib bleach pandas sortedcontainers \
      pytz setuptools mpmath bottleneck jplephem asdf-astropy pyarrow

To also be able to run tests (see below) and support :ref:`builddocs` use the
following. We use ``pip`` for these packages to ensure getting the latest
releases which are compatible with the latest ``pytest`` and ``sphinx`` releases::

    python -m pip install pytest-astropy sphinx-astropy

.. warning::

    Attempting to use `pip <https://pip.pypa.io>`__ to upgrade your installation
    of ``astropy`` itself may result in a corrupted installation.

.. _testing_installed_astropy:

Testing an Installed ``astropy``
================================

{% if is_development %}

The easiest way to test if your installed version of ``astropy`` is running
correctly is to use the :ref:`astropy.test()` function::

    import astropy
    astropy.test()

The tests should run and print out any failures, which you can report at
the `Astropy issue tracker <https://github.com/astropy/astropy/issues>`_.

This way of running the tests may not work if you do it in the ``astropy`` source
distribution. See :ref:`sourcebuildtest` for how to run the tests from the
source code directory, or :ref:`running-tests` for more details.

{%else%}

See the :ref:`latest documentation on how to test your installed version of
astropy <astropy-dev:testing_installed_astropy>`.

{%endif%}

.. _astropy-main-req:

Requirements
************

``astropy`` has the following strict requirements:

- |Python| |minimum_python_version| or later

- |NumPy| |minimum_numpy_version| or later

- |PyERFA| |minimum_pyerfa_version| or later

- `PyYAML <https://pyyaml.org>`_ |minimum_pyyaml_version| or later

- |packaging| |minimum_packaging_version| or later

``astropy`` also depends on a number of other packages for optional features.
The following are particularly recommended:

- |SciPy| |minimum_scipy_version| or later: To power a variety of features
  in several modules.

- |Matplotlib| |minimum_matplotlib_version| or later: To provide plotting
  functionality that `astropy.visualization` enhances.

The further dependencies provide more specific features:

- `h5py <http://www.h5py.org/>`_: To read/write
  :class:`~astropy.table.Table` objects from/to HDF5 files.

- `BeautifulSoup <https://www.crummy.com/software/BeautifulSoup/>`_: To read
  :class:`~astropy.table.table.Table` objects from HTML files.

- `html5lib <https://html5lib.readthedocs.io/en/stable/>`_: To read
  :class:`~astropy.table.table.Table` objects from HTML files using the
  `pandas <https://pandas.pydata.org/>`_ reader.

- `bleach <https://bleach.readthedocs.io/>`_: Used to sanitize text when
  disabling HTML escaping in the :class:`~astropy.table.Table` HTML writer.

- `xmllint <http://www.xmlsoft.org/>`_: To validate VOTABLE XML files.
  This is a command line tool installed outside of Python.

- `pandas <https://pandas.pydata.org/>`_: To convert
  :class:`~astropy.table.Table` objects from/to pandas DataFrame objects.
  Version 0.14 or higher is required to use the :ref:`table_io_pandas`
  I/O functions to read/write :class:`~astropy.table.Table` objects.

- `sortedcontainers <https://pypi.org/project/sortedcontainers/>`_ for faster
  ``SCEngine`` indexing engine with ``Table``, although this may still be
  slower in some cases than the default indexing engine.

- `pytz <https://pythonhosted.org/pytz/>`_: To specify and convert between
  timezones.

- `jplephem <https://pypi.org/project/jplephem/>`_: To retrieve JPL
  ephemeris of Solar System objects.

- `setuptools <https://setuptools.readthedocs.io>`_: Used for discovery of
  entry points which are used to insert fitters into `astropy.modeling.fitting`.

- `mpmath <https://mpmath.org/>`_: Used for the 'kraft-burrows-nousek'
  interval in `~astropy.stats.poisson_conf_interval`.

- `asdf-astropy <https://github.com/astropy/asdf-astropy>`_ |minimum_asdf_astropy_version| or later: Enables the
  serialization of various Astropy classes into a portable, hierarchical,
  human-readable representation.

- `bottleneck <https://pypi.org/project/Bottleneck/>`_: Improves the performance
  of sigma-clipping and other functionality that may require computing
  statistics on arrays with NaN values.

- `certifi <https://pypi.org/project/certifi/>`_: Useful when downloading
  files from HTTPS or FTP+TLS sites in case Python is not able to locate
  up-to-date root CA certificates on your system; this package is usually
  already included in many Python installations (e.g., as a dependency of
  the ``requests`` package).

- `pyarrow <https://arrow.apache.org/docs/python/>`_ |minimum_pyarrow_version| or later:
  To read/write :class:`~astropy.table.Table` objects from/to Parquet files.

- |fsspec| |minimum_fsspec_version| or later: Enables access to :ref:`subsets
  of remote FITS files <fits_io_cloud>` without having to download the entire file.

- |s3fs| |minimum_s3fs_version| or later: Enables access to files hosted in
  AWS S3 cloud storage.

However, note that these packages require installation only if those particular
features are needed. ``astropy`` will import even if these dependencies are not
installed.

The following packages can optionally be used when testing:

- |pytest-astropy|: See :ref:`sourcebuildtest`

- `pytest-xdist <https://pypi.org/project/pytest-xdist/>`_: Used for
  distributed testing.

- `pytest-mpl <https://github.com/matplotlib/pytest-mpl>`_: Used for testing
  with Matplotlib figures.

- `objgraph <https://mg.pov.lt/objgraph/>`_: Used only in tests to test for reference leaks.

- |IPython| |minimum_ipython_version| or later:
  Used for testing the notebook interface of `~astropy.table.Table`.

- `coverage <https://coverage.readthedocs.io/>`_: Used for code coverage
  measurements.

- `skyfield <https://rhodesmill.org/skyfield/>`_: Used for testing Solar System
  coordinates.

- `sgp4 <https://pypi.org/project/sgp4/>`_: Used for testing satellite positions.

- `tox <https://tox.readthedocs.io/en/latest/>`_: Used to automate testing
  and documentation builds.

Building from Source
********************

Prerequisites
=============

You will need a compiler suite and the development headers for Python in order
to build ``astropy``. You do not need to install any other specific build
dependencies (such as `Cython <https://cython.org/>`_) since these are
declared in the ``pyproject.toml`` file and will be automatically installed into
a temporary build environment by pip.

Prerequisites for Linux
=======================

On Linux, using the package manager for your distribution will usually be the
easiest route to making sure you have the prerequisites to build ``astropy``. In
order to build from source, you will need the Python development
package for your Linux distribution, as well as pip.

For Debian/Ubuntu::

    sudo apt-get install python3-dev python3-numpy-dev python3-setuptools cython3 python3-pytest-astropy

For Fedora/RHEL::

    sudo yum install python3-devel python3-numpy python3-setuptools python3-Cython python3-pytest-astropy

.. note:: Building the developer version of ``astropy`` may require
          newer versions of the above packages than are available in
          your distribution's repository.  If so, you could either try
          a more up-to-date distribution (such as Debian ``testing``),
          or install more up-to-date versions of the packages using
          ``pip`` or ``conda`` in a virtual environment.

Prerequisites for Mac OS X
==========================

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
=============================

Source Packages
---------------

The latest stable source package for ``astropy`` can be `downloaded here
<https://pypi.org/project/astropy>`_.

Development Repository
----------------------

The latest development version of ``astropy`` can be cloned from GitHub
using this command::

   git clone https://github.com/astropy/astropy.git

If you wish to participate in the development of ``astropy``, see the
:ref:`developer-docs`. The present document covers only the basics necessary to
installing ``astropy``.

Building and Installing
=======================

To build and install ``astropy`` (from the root of the source tree)::

    python -m pip install .

If you install in this way and you make changes to the code, you will need to
re-run the install command for changes to be reflected. Alternatively, you can
use::

    python -m pip install --editable .

which installs ``astropy`` in develop/editable mode -- this then means that
changes in the code are immediately reflected in the installed version.

Troubleshooting
===============

If you get an error mentioning that you do not have the correct permissions to
install ``astropy`` into the default ``site-packages`` directory, you can try
installing with::

    python -m pip install --user .

which will install into a default directory in your home directory.

If you get an error directing use of option ``-std=c99`` or ``-std=gnu99``, you can try
installing with::

    CFLAGS='-std=c99' python -m pip install --editable .

This is necessary for use with CentOS7.

.. _external_c_libraries:

External C Libraries
--------------------

The ``astropy`` source ships with the C source code of a number of
libraries. By default, these internal copies are used to build
``astropy``. However, if you wish to use the system-wide installation of
one of those libraries, you can set environment variables with the
pattern ``ASTROPY_USE_SYSTEM_???`` to ``1`` when building/installing
the package.

For example, to build ``astropy`` using the system's expat parser
library, use::

    ASTROPY_USE_SYSTEM_EXPAT=1 python -m pip install --editable .

To build using all of the system libraries, use::

    ASTROPY_USE_SYSTEM_ALL=1 python -m pip install --editable .

The C libraries currently bundled with ``astropy`` include:

- `wcslib <https://www.atnf.csiro.au/people/mcalabre/WCS/>`_ see
  ``cextern/wcslib/README`` for the bundled version. To use the
  system version, set ``ASTROPY_USE_SYSTEM_WCSLIB=1``.

- `expat <https://libexpat.github.io/>`_ see ``cextern/expat/README`` for the
  bundled version. To use the system version, set ``ASTROPY_USE_SYSTEM_EXPAT=1``.

.. _install_astropy_nightly:

Installing pre-built Development Versions of ``astropy``
========================================================

Most nights a development snapshot of ``astropy`` will be compiled.
This is useful if you want to test against a development version of astropy but
do not want to have to build it yourselves. You can see the
`available astropy dev snapshots page <https://anaconda.org/astropy/astropy/files?type=pypi>`_
to find out what is currently being offered.

Installing these "nightlies" of ``astropy`` can be achieved by using ``pip``::

  python -m pip install --upgrade --index-url https://pypi.anaconda.org/astropy/simple astropy --pre

The extra index URL tells ``pip`` to check the ``pip`` index on
pypi.anaconda.org, where the nightlies are stored, and the ``--pre`` command
tells ``pip`` to install pre-release versions (in this case ``.dev`` releases).

.. _builddocs:

Building Documentation
======================

.. note::

    Building the documentation is in general not necessary unless you are
    writing new documentation or do not have internet access, because
    the latest (and archive) versions of Astropy's documentation should
    be available at `docs.astropy.org <https://docs.astropy.org>`_ .

Dependencies
------------

Building the documentation requires the ``astropy`` source code and some
additional packages. The easiest way to build the documentation is to use `tox
<https://tox.readthedocs.io/en/latest/>`_ as detailed in
:ref:`astropy-doc-building`. If you are happy to do this, you can skip the rest
of this section.

On the other hand, if you wish to call Sphinx manually to build the
documentation, you will need to make sure that a number of dependencies are
installed. If you use conda, the easiest way to install the dependencies is
with::

    conda install --channel conda-forge sphinx-astropy

Without conda, you install the dependencies by specifying ``[docs]`` when
installing ``astropy`` with pip::

    python -m pip install --editable ".[docs]"

You can alternatively install the `sphinx-astropy
<https://github.com/astropy/sphinx-astropy>`_ package with pip::

    python -m pip install sphinx-astropy

In addition to providing configuration common to packages in the Astropy
ecosystem, this package also serves as a way to automatically get the main
dependencies, including:

* `Sphinx <http://www.sphinx-doc.org>`_ - the main package we use to build
  the documentation
* `astropy-sphinx-theme <https://github.com/astropy/astropy-sphinx-theme>`_ -
  the default 'bootstrap' theme used by ``astropy`` and a number of affiliated
  packages
* `sphinx-automodapi <https://sphinx-automodapi.readthedocs.io>`_ - an extension
  that makes it easy to automatically generate API documentation
* `sphinx-gallery <https://sphinx-gallery.readthedocs.io/en/latest/>`_ - an
  extension to generate example galleries
* |numpydoc| - an extension to parse
  docstrings in NumPyDoc format
* `pillow <https://pillow.readthedocs.io>`_ - used in one of the examples
* `Graphviz <http://www.graphviz.org>`_ - generate inheritance graphs (available
  as a conda package or a system install but not in pip)

.. Note::
    Both of the ``pip`` install methods above do not include `Graphviz
    <http://www.graphviz.org>`_.  If you do not install this package separately
    then the documentation build process will produce a very large number of
    lengthy warnings (which can obscure bona fide warnings) and also not
    generate inheritance graphs.

.. _astropy-doc-building:

Building
--------

There are two ways to build the Astropy documentation. The easiest way is to
execute the following tox command (from the ``astropy`` source directory)::

    tox -e build_docs

If you do this, you do not need to install any of the documentation dependencies
as this will be done automatically. The documentation will be built in the
``docs/_build/html`` directory, and can be read by pointing a web browser to
``docs/_build/html/index.html``.

Alternatively, you can do::

    cd docs
    make html

.. note::
   If you have a multi-core processor, and wish to leverage this for building
   documentation, you can do so as follows::

       cd docs
       SPHINXOPTS="-j N" make html

   where ``N`` is the number of processes over which to distribute the build, as
   described in the `sphinx-build Documentation
   <https://www.sphinx-doc.org/en/master/man/sphinx-build.html#cmdoption-sphinx-build-j>`_.

The documentation will be generated in the same location. Note that
this uses the installed version of astropy, so if you want to make sure
the current repository version is used, you will need to install it with
e.g.::

    python -m pip install --editable ".[docs]"

before changing to the ``docs`` directory.

In the second way, LaTeX documentation can be generated by using the command::

    make latex

The LaTeX file ``Astropy.tex`` will be created in the ``docs/_build/latex``
directory, and can be compiled using ``pdflatex``.

Reporting Issues/Requesting Features
------------------------------------

As mentioned above, building the documentation depends on a number of Sphinx
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
==========================================

{% if is_development %}

The easiest way to run the tests in a source checkout of ``astropy``
is to use `tox <https://tox.readthedocs.io/en/latest/>`_::

    tox -e test-alldeps

There are also alternative methods of :ref:`running-tests` if you
would like more control over the testing process.

{%else%}

See the :ref:`latest documentation on how to run the tests in a source
checkout of astropy <astropy-dev:sourcebuildtest>`

{%endif%}
