astropy.extern
==============

This sub-package contains third-party Python packages/modules that are
required for some of Astropy's core functionality.  It also contains third-
party JavaScript libraries used for browser-based features.

In particular, this currently includes for Python:

- ConfigObj_: This provides the core config file handling for Astropy's
  configuration system.

- PLY_: This is a parser generator providing lex/yacc-like tools in Python.
  It is used for Astropy's unit parsing and angle/coordinate string parsing.

Notes for third-party packagers
-------------------------------

Packagers preparing Astropy for inclusion in packaging frameworks have
different options for how to handle these third-party extern packages, if they
would prefer to use their system packages rather than the bundled versions.

jQuery/DataTables
^^^^^^^^^^^^^^^^^

It is possible to change the default urls for the remote versions of these
files by using the Astropy
`Configuration system <https://docs.astropy.org/en/stable/config/>`_. The default
configuration file (``$HOME/.astropy/config``) contains a commented section
``[table.jsviewer]`` with two items for jQuery_ and DataTables_. It is also
possible to display the default value and modify it by importing the
configuration module::

    In [1]: from astropy.table.jsviewer import conf

    In [2]: conf.jquery_url
    Out[2]: u'https://code.jquery.com/jquery-1.11.3.min.js'

    In [3]: conf.jquery_url = '...'

Third-party packagers can override the defaults for these configuration items
by modifying the configuration objects in ``astropy/table/jsviewer.py``, or
provide astropy config files that include the overrides appropriate for the
packaged version.


Other
^^^^^

To replace any of the other Python modules included in this package, simply
remove them and update any imports in Astropy to import the system versions
rather than the bundled copies.


.. _ConfigObj: https://github.com/DiffSK/configobj
.. _PLY: http://www.dabeaz.com/ply/
.. _jQuery: http://jquery.com/
.. _DataTables: http://www.datatables.net/
