************
Installation
************

Requirements
============

Astropy has the following strict requirements:

- `Python <http://www.python.org/>`_ 2.6, 2.7, 3.1 or 3.2

- `Numpy <http://www.numpy.org/>`_ 1.4 or later

Astropy also depends on other packages for optional features:

- `h5py <http://alfven.org/wp/hdf5-for-python/>`_: To read/write
  :class:`~astropy.table.table.Table` objects from/to HDF5 files

- `scipy <http://www.scipy.org/>`_: To power a variety of features (currently
  mainly cosmology-related functionality)

- `xmllint <http://www.xmlsoft.org/>`_: To validate VOTABLE XML files.

However, note that these only need to be installed if those particular features
are needed. Astropy will import even if these dependencies are not installed.

.. TODO: Link to the planned dependency checker/installer tool.

Installing Astropy
==================

Using `pip`
-----------

To install Astropy with `pip`, simply run::

    pip install --no-deps astropy

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
    to install the package into your home directory.  You can read more about
    how to do this in the `pip documentation <http://www.pip-installer.org/en/1.2.1/other-tools.html#using-pip-with-the-user-scheme>`_.

    Alternatively, if you intend to do development on other software that uses
    Astropy, such as an affiliated package, consider installing Astropy into a
    :ref:`virtualenv<using-virtualenv>`.

    Do **not** install Astropy or other third-party packages using ``sudo``
    unless you are fully aware of the risks.


Binary installers
-----------------

Binary installers are available on Windows for Python 2.6, 2.7, 3.1, and 3.2
at `PyPI <https://pypi.python.org/pypi/astropy>`_.

.. _testing_installed_astropy:

Testing an installed Astropy
----------------------------

The easiest way to test your installed version of astropy is running
correctly is to use the :func:`astropy.test` function::

    import astropy
    astropy.test()

The tests should run and print out any failures, which you can report at
the `Astropy issue tracker <http://github.com/astropy/astropy/issues>`_.

.. note::
    
    This way of running the tests may not work if you do it in the
    astropy source distribution.  See :ref:`sourcebuildtest` for how to
    run the tests from the source code directory, or :ref:`running-tests`
    for more details.



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

You will also need `Cython <http://cython.org/>`_ installed to build
from source, unless you are installing a numbered release. (The releases
packages have the necessary C files packaged with them, and hence do not
require Cython.)

.. note:: If you are using MacOS X, you will need to the XCode command line
          tools.  One way to get them is to install `XCode
          <https://developer.apple.com/xcode/>`_. If you are using OS X 10.7
          (Lion) or later, you must also explicitly install the command line
          tools. You can do this by opening the XCode application, going to
          **Preferences**, then **Downloads**, and then under **Components**,
          click on the Install button to the right of **Command Line Tools**.
          Alternatively, on 10.7 (Lion) or later, you do not need to install
          XCode, you can download just the command line tools from
          https://developer.apple.com/downloads/index.action (requires an Apple
          developer account).

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
`--use-system-X` flags to the `setup.py build` command.

For example, to build Astropy using the system `libexpat`, use::

    python setup.py build --use-system-expat

To build using all of the system libraries, use::

    python setup.py build --use-system-libraries

To see which system libraries Astropy knows how to build against, use::

    python setup.py build --help

As with all distutils commandline options, they may also be provided
in a `setup.cfg` in the same directory as `setup.py`.  For example, to
use the system `libexpat`, add the following to the `setup.cfg` file::

    [build]
    use_system_expat=1

Compatibility packages
^^^^^^^^^^^^^^^^^^^^^^

.. warning:: This feature is still experimental, and you may run into
             unexpected issues with other packages, so we strongly
             recommend simply updating your code to use Astropy if
             possible, rather than rely on these compatibility packages.

Optionally, it is possible to install 'compatibility' packages that
emulate the behavior of previous packages that have now been
incorporated into Astropy. These are:

* `PyFITS <http://www.stsci.edu/institute/software_hardware/pyfits/>`_
* `vo <https://trac.assembla.com/astrolib/>`_
* `PyWCS <https://trac.assembla.com/astrolib/>`_

If you build Astropy with::

    python setup.py build --enable-legacy
    python setup.py install

or simply::

    python setup.py install --enable-legacy

then you will be able to import these modules from your scripts as if
the original packages had been installed. Using::

    import pyfits
    import vo
    import pywcs

will then be equivalent to::

    from astropy.io import fits as pyfits
    from astropy.io import vo
    from astropy import wcs as pywcs

In order to install the compatibility packages none of the
original packages should be present.

.. note:: If you are interested in testing out existing code with Astropy
          without modifying the import statements, but don't want to
          uninstall existing packages, you can use `virtualenv
          <http://www.virtualenv.org/>`_ to set up a clean environment.

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

There are two ways to build the Astropy documentation. The most straightforward
way is to execute the command (from the astropy source directory)::

    python setup.py build_sphinx

The documentation will be built in the ``docs/_build/html`` directory, and can
be read by pointing a web browser to ``docs/_build/html/index.html``.

The above method builds the API documentation from the source code.
Alternatively, you can do::

    cd docs
    make html

And the documentation will be generated in the same location, but using the
*installed* version of Astropy.

.. _sourcebuildtest:

Testing a source code build of Astropy
--------------------------------------

The easiest way to test that your Astropy built correctly (without
installing astropy) is to run this from the root of the source tree::

    python setup.py test

There are also alternative methods of :ref:`running-tests`.
