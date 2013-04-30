============
Known Issues
============

While most bugs and issues are managed using the `astropy issue
tracker <https://github.com/astropy/astropy/issues>`_, this document
lists issues that are too difficult to fix, may require some
intervention from the user to workaround, or are due to bugs in other
projects or packages.

Some docstrings can not be displayed in IPython < 0.13.2
--------------------------------------------------------

Displaying long docstrings that contain Unicode characters may fail on
some platforms in the IPython console (prior to IPython version
0.13.2)::

    In [1]: import astropy.units as u

    In [2]: u.Angstrom?
    ERROR: UnicodeEncodeError: 'ascii' codec can't encode character u'\xe5' in
    position 184: ordinal not in range(128) [IPython.core.page]

This can be worked around by changing the default encoding to `utf-8`
by adding the following to your `sitecustomize.py` file::

    import sys
    sys.setdefaultencoding('utf-8')

Note that in general, `this is not recommended
<http://ziade.org/2008/01/08/syssetdefaultencoding-is-evil/>`_,
because it can hide other Unicode encoding bugs in your application.
However, in general if your application does not deal with text
processing and you just want docstrings to work, this may be
acceptable.

The IPython issue: https://github.com/ipython/ipython/pull/2738

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

mmap support on GNU Hurd
------------------------

On Hurd and possibly other platforms ``flush()`` on memory-mapped files is not
implemented, so writing changes to a mmap'd file may not be reliable and is
this disabled.  Attempting to open a file in writeable mode with mmap will
result in a warning (and mmap will be disabled on the file automatically).

See: https://github.com/astropy/astropy/issues/968
