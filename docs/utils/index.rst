.. _utils:

************************************************
Astropy Core Package Utilities (`astropy.utils`)
************************************************

Introduction
============

The `astropy.utils` package contains general-purpose utility functions and
classes.  Examples include data structures, tools for downloading and caching
from URLs, and version intercompatibility functions.

This functionality is not astronomy-specific, but is intended primarily for
use by Astropy developers. It is all safe for users to use, but the functions
and classes are typically more complicated or specific to a particular need of
Astropy.

Because of the mostly standalone and grab-bag nature of these utilities, they
are generally best understood through their docstrings, and hence this
documentation generally does not have detailed sections like the other packages.
The exceptions are below:

.. toctree::
   :maxdepth: 1

   data
   masked/index

:ref:`utils-iers`

.. note:: The ``astropy.utils.compat`` subpackage is not included in this
    documentation. It contains utility modules for compatibility with
    older/newer versions of python and numpy, as well as including some
    bugfixes for the stdlib that are important for ``astropy``. It is recommended
    that developers at least glance over the source code for this subpackage,
    but most of it cannot be reliably included here because of the large
    amount of version-specific code it contains. Its content is solely for
    internal use of ``astropy`` and subject to changes without deprecations.
    Do not use it in external packages or code.

Reference/API
=============

.. toctree::
   :maxdepth: 2

   ref_api
