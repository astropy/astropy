astropy.extern
==============

This sub-package contains third-party Python packages/modules that are
required for some of Astropy's core functionality.  It also contains third-
party JavaScript libraries used for browser-based features.

In particular, this currently includes for Python:

- ConfigObj_ (Python 2 and 3 versions): This provides the core config file
  handling for Astropy's configuration system.

- PLY_: This is a parser generator providing lex/yacc-like tools in Python.
  It is used for Astropy's unit parsing and angle/coordinate string parsing.

- pytest_: This is the test framework used to collect and run Astropy's tests.
  It is a fairly large package, so it is bundled as a compressed base64 string
  that is unpacked and loaded as a Python module.  This packed bundled version
  is officially provided by the pytest project.

- Six_: This is a Python 2/3 compatibility library used to ease development
  of a simultaneously Python 2 and 3-compatible code base.

And for JavaScript:

- jQuery_: This is used currently for the browser-based table viewer feature.

- DataTables_: This is a plug-in for jQuery used also for the browser-based
  table viewer.

Notes for developers
--------------------

jQuery/DataTables
^^^^^^^^^^^^^^^^^
the minified files are the ones that are used in the table viewer feature, but
the non-minified versions are also present in the ``js/`` sub-directory for
packaging reasons. These files must also be distributed, to provide the source
files from which the minified ones can be compiled. This is a requirement for
Linux distributions such as Debian and Fedora.


Notes for third-party packagers
-------------------------------

Packagers preparing Astropy for inclusion in packaging frameworks have
different options for how to handle these third-party extern packages, if they
would prefer to use their system packages rather than the bundled versions.

Six
^^^

Because use of the ``six`` module is pervasive throughout the Astropy package
(it is crucial for providing Python 2/3 compatibility at all levels) it has
been made easier to replace.  The actual ``six`` module is included under
``astropy/extern/bundled``.  Importing ``astropy.extern.six`` will first try
to load the ``six`` module from that location, and then will fall back on
trying to import an installed ``six`` module from the standard ``sys.path``
locations (if it is a new enough version).

Packagers may wish to prioritize the system version and/or remove the bundled
version altogether.  To do this, edit the ``_SIX_SEARCH_PATH`` list in
``astropy.extern.six``.  If removing the bundled copy altogether simply delete
``astropy/extern/bundled/six.py`` and update this variable to read::

    _SIX_SEARCH_PATH = ['six']

No other imports in Astropy need to be updated.  Imports from
``astropy.extern.six`` are automatically updated to point to the system version
of the module.


jQuery/DataTables
^^^^^^^^^^^^^^^^^

Packagers may either use system copies of these JavaScript modules, or require
use of online versions (perhaps via URLs of cloud-hosted versions of these
modules).

It is possible to change the default urls for the remote versions of these
files by using the Astropy
`Configuration system <http://docs.astropy.org/en/stable/config/>`_. The default
configuration file (``$HOME/.astropy/config``) contains a commented section
``[table.jsviewer]`` with two items for jQuery and DataTables. It is also
possible to display the default value and modify it by importing the
configuration module::

    In [1]: from astropy.table.jsviewer import conf

    In [2]: conf.jquery_url
    Out[2]: u'https://code.jquery.com/jquery-1.11.3.min.js'

    In [3]: conf.jquery_url = '...'

Third-party packagers can override the defaults for these configuration items
(by modifying the configuration objects in ``astropy/table/jsviewer.py``, or
provide astropy config files that include the overrides appropriate for the
packaged version.  They would *also* need to set the default
``use_local_files`` option to ``False`` for these settings to be read.


Other
^^^^^

To replace any of the other Python modules included in this package, simply
remove them and update any imports in Astropy to import the system versions
rather than the bundled copies.  If requested we could add a mechanism similar
to that used by ``six`` for the other modules.


.. _ConfigObj: https://github.com/DiffSK/configobj
.. _PLY: http://www.dabeaz.com/ply/
.. _pytest: http://pytest.org/latest/
.. _Six: http://pypi.python.org/pypi/six/
.. _jQuery: http://jquery.com/
.. _DataTables: http://www.datatables.net/
