.. _installing-astropy:

************
Installation
************

Overview
========

The first step to installing ``astropy`` is to ensure that you have a Python
environment which is **isolated** from your system Python installation. This is
important because ``astropy`` has many dependencies, and you do not want to accidentally
break your system by installing incompatible versions of these dependencies.

For this installation guide we use the `conda <https://docs.conda.io/en/latest/>`_
package manager provided by `miniforge <https://github.com/conda-forge/miniforge>`_.
This is a popular choice and works well, especially for newcomers. It is easy to install
and use on all platforms and it makes it easy to install the latest Python version. If
you already have a ``miniforge``-based Python environment then you can skip to
:ref:`installing-astropy-with-pip`.

Another option for more experienced users is a virtual environment manager such as the
Python standard library `venv <https://docs.python.org/3/library/venv.html>`_ module.
There are numerous resources available to help you set up a virtual environment in this
manner if you choose this option.

.. note::
   We **do not recommend** using ``astropy`` with an existing `miniconda
   <https://docs.anaconda.com/miniconda/>`_ or `Anaconda Python
   <https://www.anaconda.com/download/>`_ distribution. The ``astropy`` package provided
   by Anaconda Inc. in the ``defaults`` channel can be outdated and these distributions
   can require a license for use at a large organisation. Instead, use ``miniforge`` as
   described below.

Once you have a Python environment set up, you will install ``astropy`` using |pip| or
|conda|. Here we document using |pip| because it is easier to install the optional
dependencies, but feel free to use |conda| if you prefer.

Install ``miniforge``
=====================

You will install Python by first installing `miniforge
<https://github.com/conda-forge/miniforge/#miniforge>`__. This provides the `conda
package manager <https://docs.conda.io/en/latest/>`_ with the default remote package
repository set to the community-led `conda-forge <https://conda-forge.org>`_ channel.

In a new terminal (miniforge Prompt on Windows) run ``conda list`` to test that the
install has worked.

Create Python Environment
=========================

To create a new Python environment for ``astropy`` and other packages, start by
launching a terminal (under a UNIX-like system) or the miniforge Prompt (under Windows).
Now we will create and activate a new virtual environment to install ``astropy`` into:

.. code-block:: bash

    $ conda create --channel conda-forge  --name astropy python
    $ conda activate astropy

In this case the environment we have created is named ``astropy`` but you can use any
name you like.

In the future when you make a new terminal, you will need to run ``conda activate
astropy`` to activate this environment.

.. _installing-astropy-with-pip:

Install ``astropy``
===================

You can install ``astropy`` and the rest of your dependencies using either |pip| or
|conda|. Both methods are fully supported and will work well.

.. warning::
   Once you have created your base Python environment with |conda|, you should try to
   stick with one method for installing new packages in your environment. In particular,
   |conda| is not aware of packages installed with |pip| and may overwrite them.

Using pip
---------
To install ``astropy`` and your choice of :ref:`dependencies <astropy-main-req>`, run
one of the following commands::

    python -m pip install astropy                # Minimum required dependencies
    python -m pip install "astropy[recommended]" # Recommended dependencies
    python -m pip install "astropy[all]"         # All optional dependencies
    python -m pip install "astropy[dev_all]"     # All optional and test dependencies

In most cases, this will install a pre-compiled version of ``astropy`` (called a
*wheel*). However, if you are installing astropy on an uncommon platform, astropy will be
installed from a source file. In this unusual case you will need a C compiler to be
installed (see `Build from source`_ below) for the installation to succeed.

.. warning:: Do **not** install ``astropy`` or other packages using ``sudo`` or any
   elevated privilege.

Using conda
-----------
To install ``astropy`` and the minimal set of required dependencies, run::

  conda install --channel conda-forge astropy

Install the recommended dependencies with::

  conda install --channel conda-forge scipy matplotlib

Install the optional dependencies with::

  conda install --channel conda-forge ipython jupyter dask h5py pyarrow \
     beautifulsoup4 html5lib bleach pandas sortedcontainers pytz jplephem mpmath \
     asdf-astropy bottleneck fsspec s3fs certifi

Testing
-------
You can test that your newly installed version of ``astropy`` is working via the
`documentation on how to test your installed version of astropy
<https://docs.astropy.org/en/latest/development/testguide.html#running-tests-installed-astropy>`_.

.. _astropy-main-req:

Requirements
============

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

- `ipydatagrid <https://pypi.org/project/ipydatagrid/>`_: Used in
  :meth:`astropy.table.Table.show_in_notebook` to display the Astropy table
  in Jupyter notebook for ``backend="ipydatagrid"``.

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

Build from Source
=================

If you want to build the code from source, follow the instructions for
:ref:`contributing_environment`. Note that instead of cloning from your fork, you can
choose to clone from the main repository::

    git clone https://github.com/astropy/astropy.git
    cd astropy

Building the documentation is typically not necessary unless you are
developing code or documentation or do not have internet access, because
the stable, latest, and archived versions of Astropy's documentation are
available at `docs.astropy.org <https://docs.astropy.org>`_ . The process
is described in `Building the Documentation from Source <https://docs.astropy.org/en/latest/development/docguide.html#builddocs>`_.

.. _sourcebuildtest:

Test Source Code Build
----------------------

{% if is_development %}

The easiest way to run the tests in a source checkout of ``astropy``
is to use `tox <https://tox.readthedocs.io/en/latest/>`_::

    tox -e test-alldeps

There are also alternative methods of :ref:`running-tests` if you
would like more control over the testing process.

{%else%}

See the `latest documentation on how to run the tests in a source
checkout of astropy <https://docs.astropy.org/en/latest/install.html#testing-a-source-code-build-of-astropy>`_.

{%endif%}


.. _install_astropy_nightly:

Install Pre-built Development Version
=====================================

Most nights a development snapshot of ``astropy`` will be compiled.
This is useful if you want to test against a development version of astropy but
do not want to have to build it yourselves. You can see the
`available astropy dev snapshots page <https://anaconda.org/astropy/astropy/files?type=pypi>`_
to find out what is currently being offered.

Installing these "nightlies" of ``astropy`` can be achieved by using ``pip``::

  python -m pip install --upgrade --extra-index-url https://pypi.anaconda.org/astropy/simple astropy --pre

The extra index URL tells ``pip`` to check the ``pip`` index on
pypi.anaconda.org, where the nightlies are stored, and the ``--pre`` command
tells ``pip`` to install pre-release versions (in this case ``.dev`` releases).

You can test this installation by running the tests as described in the section
`Running tests on an installed astropy <https://docs.astropy.org/en/latest/development/testguide.html#running-tests-installed-astropy>`_.
