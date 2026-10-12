astropy.extern
==============

This sub-package contains third-party Python packages/modules that are
required for some of Astropy's core functionality.

In particular, this currently includes for Python:

- PLY_: This is a parser generator providing lex/yacc-like tools in Python.
  It is used for Astropy's unit parsing and angle/coordinate string parsing.

Notes for third-party packagers
-------------------------------

Packagers preparing Astropy for inclusion in packaging frameworks have
different options for how to handle these third-party extern packages, if they
would prefer to use their system packages rather than the bundled versions.

To replace any of the other Python modules included in this package,
remove them and update any imports in Astropy to import the system versions
rather than the bundled copies.

.. _PLY: http://www.dabeaz.com/ply/
