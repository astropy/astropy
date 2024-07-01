************
Installation
************

.. _installing-astropy:

Installing ``astropy``
**********************

If you are new to Python and/or do not have familiarity with `Python virtual
environments <https://docs.python.org/3/tutorial/venv.html>`__, then we recommend
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
this case you will need a C compiler to be installed
(see `Building from source`_ below) for the installation to succeed.

If you get a ``PermissionError`` this means that you do not have the required
administrative access to install new packages to your Python installation. In
this case you should first create and activate a Python environment using either
:ref:`Conda <anaconda_install>` or a `Python virtual
environment <https://docs.python.org/3/tutorial/venv.html>`__. Both of these options
will also allow you to do development on other software that uses
``astropy``, such as an affiliated package.

.. warning:: Do **not** install ``astropy`` or other third-party packages using ``sudo``.

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

See the :ref:`documentation on how to test your installed version of
astropy <astropy-dev:running-tests-installed-astropy>`.

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

If you want to build the code from source, follow the instructions for
:ref:`contributing_environment`. Note that instead of cloning from your fork, you can
choose to clone from the main repository::

    git clone https://github.com/astropy/astropy.git
    cd astropy

Building the documentation is typically not necessary unless you are
developing code or documentation or do not have internet access, because
the stable, latest, and archived versions of Astropy's documentation are
available at `docs.astropy.org <https://docs.astropy.org>`_ . The process
is described in :ref:`astropy-dev:builddocs`.

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

You can test this installation by running the tests as described in the section
:ref:`astropy-dev:running-tests-installed-astropy`.
