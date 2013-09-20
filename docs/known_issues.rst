============
Known Issues
============

While most bugs and issues are managed using the `astropy issue
tracker <https://github.com/astropy/astropy/issues>`_, this document
lists issues that are too difficult to fix, may require some
intervention from the user to workaround, or are due to bugs in other
projects or packages.

Floating point precision issues on Python 2.6 on Microsoft Windows
------------------------------------------------------------------

When converting floating point numbers to strings on Python 2.6 on a
Microsoft Windows platform, some of the requested precision may be
lost.

The easiest workaround is to install Python 2.7.

The Python issue: http://bugs.python.org/issue7117

Failing logging tests when running the tests in IPython
-------------------------------------------------------

When running the Astropy tests using ``astropy.test()`` in an IPython
interpreter some of the tests in the ``astropy/tests/test_logger.py`` fail.
This is due to mutually incompatible behaviors in IPython and py.test, and is
not due to a problem with the test itself or the feature being tested.

See: https://github.com/astropy/astropy/issues/717

Numpy 1.4.x unreliable on 64-bit Ubuntu
---------------------------------------

As of Ubuntu 12.04 (and possibly earlier), the 1.4.x versions of numpy sometimes
cause segmentation faults.  This problem is not unique to Astropy, as the numpy
tests themselves do not pass, but it does cause some Astropy functionality to
fail.

The solution is to use a more recent version of Numpy.

mmap support for ``astropy.io.fits`` on GNU Hurd
------------------------------------------------

On Hurd and possibly other platforms ``flush()`` on memory-mapped files is not
implemented, so writing changes to a mmap'd FITS file may not be reliable and is
thus disabled.  Attempting to open a FITS file in writeable mode with mmap will
result in a warning (and mmap will be disabled on the file automatically).

See: https://github.com/astropy/astropy/issues/968

Crash on upgrading from Astropy 0.2 to a newer version
------------------------------------------------------

It is possible for installation of a new version of Astropy, or upgrading of an
existing installation to crash due to not having permissions on the
``~/.astropy/`` directory (in your home directory) or some file or subdirectory
in that directory.  In particular this can occur if you installed Astropy as
the root user (such as with ``sudo``) at any point.  This can manifest in
several ways, but the most common is a traceback ending with ``ImportError:
cannot import name config``.  To resolve this issue either run ``sudo chown -R
<your_username> ~/.astropy`` or, if you don't need anything in it you can blow
it away with ``sudo rm -rf ~/.astropy``.

See for example: https://github.com/astropy/astropy/issues/987

Color printing on Windows
-------------------------

Colored printing of log messages and other colored text does work in Windows
but only when running in the IPython console.  Colors are not currently
supported in the basic Python command-line interpreter on Windows.
