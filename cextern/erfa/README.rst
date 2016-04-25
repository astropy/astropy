This is the source code repository for ERFA (Essential Routines for
Fundamental Astronomy).  ERFA is a C library containing key algorithms for
astronomy, and is based on the `SOFA library <http://www.iausofa.org/>`_ published by the International
Astronomical Union (IAU).

ERFA is intended to replicate the functionality of SOFA (aside from possible
bugfixes in ERFA that have not yet been included in SOFA), but is licensed
under a three-clause BSD license to enable its compatibility with a wide
range of open source licenses. Permission for this release has been
obtained from the SOFA board, and is avilable in the ``LICENSE`` file included
in this source distribution.

Differences from SOFA
---------------------

This version of ERFA (v1.2.0) is based on SOFA version "20150209_a", with the
differences outlined below.

ERFA branding
^^^^^^^^^^^^^

All references to "SOFA" in the source code have been changed to ERFA, and
functions have the prefix ``era`` instead of ``iau``.

C macro prefixes
^^^^^^^^^^^^^^^^

All C macros used in ERFA are the same as their SOFA equivalents, but with an
``ERFA_`` prefix to prevent namespace collisions.

Bugfixes
^^^^^^^^

ERFA includes smaller changes that may or may not eventually make it into SOFA,
addressing localized bugs or similar smaller issues:

* Typos have been corrected in the documentation of atco13 and atio13 (see https://github.com/liberfa/erfa/issues/29).

Note: Issues identified in ERFA should always be reported upstream to SOFA
at sofa@ukho.gov.uk.

Building and installing ERFA
----------------------------

To build and install a released version of ERFA in your OS's standard
location, simply do::

    ./configure
    make
    make install

If you want to run the tests to make sure ERFA built correctly, before
installing do::

    make check


For developers
^^^^^^^^^^^^^^

If you are using a developer version from github, you will need to first do
``./bootstrap.sh`` before the above commands. This requires ``autoconf`` and
``libtool``.

If you wish to build against the ERFA static library without installing, you
will find it in ``$ERFAROOT/src/.libs/liberfa.a`` after running ``make``.

Creating a single-file version of the source code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, if you wish to bundle the ERFA source code with a separate
package, you can use the ``source_flattener.py`` script from the
`erfa-fetch repository`_ to combine
the ERFA source code into just two files: a ``erfa.c`` source file, and an
``erfa.h`` include file.  You should run this script like this::

    cd /path/to/erfa-source-code
    python /path/to/erfa-fetch/source_flattener.py src -n erfa

If possible, however, it is recommended that you provide an option to use any
copy of the ERFA library that is already installed on the system.

Travis build status
-------------------
.. image:: https://travis-ci.org/liberfa/erfa.png
    :target: https://travis-ci.org/liberfa/erfa

.. _erfa-fetch repository: https://github.com/liberfa/erfa-fetch
