=============================
`astropy.tools` documentation
=============================

The `~astropy.tools` package holds general astronomy functions or algorithms
that are likely of use to users, but either not related to functionality in
an existing package or of general use across multiple packages.

.. note::
    For functions and classes that are more developer-oriented, the correct
    package is `astropy.utils`.  `astropy.tools` is intended primarily for
    functionality that is astronomy-specific and/or of use to users.


Reference/API
-------------
Below are the reference documentation for the tools sub-packages.  All public
functions and classes in these packages will be imported into the
`astropy.tools` package, so the recommended usage is e.g.
``from astropy.tools import sigma_clip`` or ``import astropy.tools`` instead of
``from astropy.tools.misc import sigma_clip`` or similar.

`astropy.tools.misc`
^^^^^^^^^^^^^^^^^^^^^^

.. automodule::  astropy.tools.misc
    :members:
