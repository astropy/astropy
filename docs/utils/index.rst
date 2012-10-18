.. _utils:

************************************************
Astropy Core Package Utilities (`astropy.utils`)
************************************************

Introduction
============

The `astropy.utils` package contains general-purpose utilities
functions and classes that such as general-purpose data structures,
version intercompatibility functions. Basically, if it's not
astronomy-related, but likely useful for other developers, it probably
lives here. These are safe for users to make use of, but they are
typically more complicated or esoteric enough that they are mostly of
interest to developers.

Because of the mostly standalone and grab-bag nature of these utilities,
they are generally best understood through their docstrings, and hence
only the reference section is currently provided for this subpackage.

.. note::
    The `astropy.utils.compat` subpackage is not included in this
    documentation. It contains utility modules for compatibility with
    older/newer versions of python, as well as including some bugfixes
    for the stdlib that are important for Astropy. It is recommended
    that developers at least glance over the source code for this
    subpackage, but it cannot be reliably included here because of the
    large amount of version-specific code it contains.


Reference/API
=============

.. automodapi:: astropy.utils.misc
    :no-inheritance-diagram:

.. automodapi:: astropy.utils.collections
    :no-inheritance-diagram:

.. automodapi:: astropy.utils.console
    :no-inheritance-diagram:

XML
---
The ``astropy.utils.xml.*`` modules provide various
`XML <http://www.w3.org/XML/>`_ processing tools.

.. automodapi:: astropy.utils.xml.check
    :no-inheritance-diagram:
    :headings: ^"

.. automodapi:: astropy.utils.xml.iterparser
    :no-inheritance-diagram:
    :headings: ^"

.. automodapi:: astropy.utils.xml.validate
    :no-inheritance-diagram:
    :headings: ^"

.. automodapi:: astropy.utils.xml.writer
    :no-inheritance-diagram:
    :headings: ^"




