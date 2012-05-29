"""
This subpackage contains implementations of command-line scripts that are
included with PyFITS.

It lives under the pyfits package primarily for Python 3 support--2to3 does not
convert these scripts unless they are actually found in the packages (as
opposed to the top-level scripts/ directory.

The actual scripts that are installed in bin/ are simple wrappers for these
modules that will run in any Python version.
"""
