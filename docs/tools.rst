.. _tools:

Astronomy-Oriented Tools (`astropy.tools`)
==========================================

Introduction
------------

The `astropy.tools` package holds smallish general astronomy functions
or algorithms that are likely of use to users, but either not related to
functionality in an existing package or of general use across multiple
packages.



Getting Started
---------------

The current tools are fairly self-contained, and include relevant examples in
their docstrings.  For example, see `~astropy.tools.misc.sigma_clip`.



.. NOTE TO FUTURE DEVS: When this subpackage gets more populated, it will be
.. wise to add a section of the form shown below.  Be sure to move this file
.. to docs/tools/index.rst and update docs/index.rst to tools/index when
.. that happens.

.. Using packagename
.. -----------------

.. For more complicated packages that require multiple documents, this
.. should just be a table of contents referencing those documents:

.. .. toctree::
..     packagename/subdoc1
..     packagename/subdoc2
..     packagename/subdoc3


See Also
--------

* :ref:`utils`
    The subpackage for utilities that are oriented towards developers, rather than
    users.  These utilities are more general-purposes, while `astropy.tools` is
    more astronomy-focused.


Reference/API
-------------

.. automodapi:: astropy.tools
