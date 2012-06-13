************
Installation
************

Requirements
============

Astropy has the following strict requirements:

- `Python <http://www.python.org/>`_ 2.6, 2.7, 3.1 or 3.2

- `Numpy <http://www.numpy.org/>`_ 1.4 or later

Astropy also depends on other projects for optional features.

- `xmllint <http://www.xmlsoft.org/>`_: To validate VOTABLE XML files.

TODO: Link to the planned dependency checker/installer tool.

Installing Astropy
==================

Using `pip`
-----------

TODO: Write once we’re up on PyPI.

Binary installers
-----------------

TODO: Write about where to obtain binary packages (.dmg, .msi etc.)

Building from source
====================

Prerequisites
-------------

You will need a compiler suite and the development headers for Python
and Numpy in order to build Astropy.  Using the package manager for
your platform will usually be the easiest route.

The `instructions for building Numpy from source
<http://docs.scipy.org/doc/numpy/user/install.html>`_ are also a good
resource for setting up your environment to build Python packages.

Obtaining the source
--------------------

Source packages
^^^^^^^^^^^^^^^

Source tarballs of past releases and the current development branch of
astropy can be downloaded from here:

   https://github.com/astropy/astropy/downloads

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
^^^^^^^^^^^^^^^^^^^^^^^

Astropy uses the Python `distutils framework
<http://docs.python.org/install/index.html>`_ for building and
installing.

To build Astropy (from the root of the source tree)::

    python setup.py build

To install Astropy (from the root of the source tree)::

    python setup.py install

.. _builddocs:

Building documentation
^^^^^^^^^^^^^^^^^^^^^^

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

Testing Astropy
^^^^^^^^^^^^^^^

The easiest way to test that your Astropy built correctly is to run
(from the root of the source tree)::

    python setup.py test

There are also alternative methods of :ref:`running-tests`.
