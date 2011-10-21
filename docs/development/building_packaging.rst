==========================================
Building, Packaging and Distribution Guide
==========================================

There is a central `setup.py`.  It defines which Python packages to
install.  Each package does not have its own standalone `setup.py`.

Each package that needs to build C extensions has a module
`setup_package.py` that contains a function `get_extensions()` which
returns a list of `distutils.core.Extension` objects defining any
extensions to be built.

There are a set of helper functions for commonly occurring things when
building C extensions (e.g. finding the Numpy headers and library) in
`astropy.setup_helpers`.

Future directions
-----------------

We plan to use packaging/distutils2 when it is ready, but in the
meantime will be using distutils.
