.. doctest-skip-all

============
Known Issues
============

While most bugs and issues are managed using the `astropy issue
tracker <https://github.com/astropy/astropy/issues>`_, this document
lists issues that are too difficult to fix, may require some
intervention from the user to workaround, or are due to bugs in other
projects or packages.

.. _quantity_issues:

Quantities lose their units with some operations
------------------------------------------------

Quantities are subclassed from numpy's `~numpy.ndarray` and in some numpy operations
(and in scipy operations using numpy internally) the subclass is ignored, which
means that either a plain array is returned, or a `~astropy.units.quantity.Quantity` without units.
E.g.::

    >>> import astropy.units as u

    >>> import numpy as np

    >>> q = u.Quantity(np.arange(10.), u.m)

    >>> np.dot(q,q)
    285.0

    >>> np.hstack((q,q))
    <Quantity [ 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 0., 1., 2., 3., 4.,
                5., 6., 7., 8., 9.] (Unit not initialised)>

Work-arounds are available for some cases.  For the above::

    >>> q.dot(q)
    <Quantity 285.0 m2>

    >>> u.Quantity([q, q]).flatten()
    <Quantity [ 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 0., 1., 2., 3., 4.,
                5., 6., 7., 8., 9.] m>

An incomplete list of specific functions which are known to exhibit this behavior follows.

* `numpy.dot`
* `numpy.hstack`, `numpy.vstack`, ``numpy.c_``, ``numpy.r_``, `numpy.append`
* `numpy.where`
* `numpy.choose`
* `numpy.vectorize`
* pandas DataFrame(s)


See: https://github.com/astropy/astropy/issues/1274

Quantities float comparison with np.isclose fails
-------------------------------------------------

Comparing Quantities floats using the numpy function `~numpy.isclose` fails on
numpy 1.9 as the comparison between ``a`` and ``b`` is made using the formula

.. math::

    |a - b| \le (a_\textrm{tol} + r_\textrm{tol} \times |b|)

This will result in the following traceback when using this with Quantities::

    >>> from astropy import units as u, constants as const
    >>> import numpy as np
    >>> np.isclose(500* u.km/u.s, 300 * u.km / u.s)
    UnitsError: Can only apply 'add' function to dimensionless quantities when
    other argument is not a quantity (unless the latter is all zero/infinity/nan)

An easy solution is::

    >>> np.isclose(500* u.km/u.s, 300 * u.km / u.s, atol=1e-8 * u.mm / u.s)
    array([False], dtype=bool)



Some docstrings can not be displayed in IPython < 0.13.2
--------------------------------------------------------

Displaying long docstrings that contain Unicode characters may fail on
some platforms in the IPython console (prior to IPython version
0.13.2)::

    >>> import astropy.units as u

    >>> u.Angstrom?
    ERROR: UnicodeEncodeError: 'ascii' codec can't encode character u'\xe5' in
    position 184: ordinal not in range(128) [IPython.core.page]

This can be worked around by changing the default encoding to ``utf-8``
by adding the following to your ``sitecustomize.py`` file::

    import sys
    sys.setdefaultencoding('utf-8')

Note that in general, `this is not recommended
<http://ziade.org/2008/01/08/syssetdefaultencoding-is-evil/>`_,
because it can hide other Unicode encoding bugs in your application.
However, in general if your application does not deal with text
processing and you just want docstrings to work, this may be
acceptable.

The IPython issue: https://github.com/ipython/ipython/pull/2738

Locale errors
-------------

On MacOS X, you may see the following error when running ``setup.py``::

      ...
    ValueError: unknown locale: UTF-8

You may also (on MacOS X or other platforms) see errors such as::

      ...
      stderr = stderr.decode(stdio_encoding)
    TypeError: decode() argument 1 must be str, not None

This is due to the ``LC_CTYPE`` environment variable being incorrectly set to
``UTF-8`` by default, which is not a valid locale setting. To fix this, set
this environment variable, as well as the ``LANG`` and ``LC_ALL`` environment
variables to e.g. ``en_US.UTF-8`` using, in the case of ``bash``::

    export LANG="en_US.UTF-8"
    export LC_ALL="en_US.UTF-8"
    export LC_CTYPE="en_US.UTF-8"

To avoid any issues in future, you should add this line to your e.g.
``~/.bash_profile`` or ``.bashrc`` file.

To test these changes, open a new terminal and type ``locale``, and you should
see something like::

    $ locale
    LANG="en_US.UTF-8"
    LC_COLLATE="en_US.UTF-8"
    LC_CTYPE="en_US.UTF-8"
    LC_MESSAGES="en_US.UTF-8"
    LC_MONETARY="en_US.UTF-8"
    LC_NUMERIC="en_US.UTF-8"
    LC_TIME="en_US.UTF-8"
    LC_ALL="en_US.UTF-8"

If so, you can go ahead and try running ``setup.py`` again (in the new
terminal).

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

Table sorting can silently fail on MacOS X or Windows with Python 3 and Numpy < 1.6.2
-------------------------------------------------------------------------------------

In Python 3, prior to Numpy 1.6.2, there was a bug (in Numpy) that caused
sorting of structured arrays to silently fail under certain circumstances (for
example if the Table contains string columns) on MacOS X, Windows, and possibly
other platforms other than Linux.  Since ``Table.sort`` relies on Numpy to
internally sort the data, it is also affected by this bug.  If you are using
Python 3, and need the sorting functionality for tables, we recommend updating
to a more recent version of Numpy.

Anaconda users should upgrade with ``conda``, not ``pip``
---------------------------------------------------------

Upgrading Astropy in the anaconda python distribution using ``pip`` can result
in a corrupted install with a mix of files from the old version and the new
version. Anaconda users should update with ``conda update astropy``. There
may be a brief delay between the release of Astropy on PyPI and its release
via the ``conda`` package manager; users can check the availability of new
versions with ``conda search astropy``.

Installation fails on Mageia-2 or Mageia-3 distributions
--------------------------------------------------------

Building may fail with warning messages such as::

    unable to find 'pow' or 'sincos'

at the linking phase. Upgrading the OS packages for Python should
fix the issue, though an immediate workaround is to edit the file::

    /usr/lib/python2.7/config/Makefile

and search for the line that adds the option ``-Wl,--no-undefined`` to the
``LDFLAGS`` variable and remove that option.


Remote data utilities in `astropy.utils.data` fail on some Python distributions
-------------------------------------------------------------------------------

The remote data utilities in `astropy.utils.data` depend on the Python
standard library `shelve` module, which in some cases depends on the
standard library `bsddb` module. Some Python distributions, including but
not limited to

* OS X, Python 2.7.5 via homebrew
* Linux, Python 2.7.6 via conda [#]_
* Linux, Python 2.6.9 via conda

are built without support for the ``bsddb`` module, resulting in an error
such as::

    ImportError: No module named _bsddb

One workaround is to install the ``bsddb3`` module.

.. [#] Continuum `says
       <https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/mCQL6fVx55A>`_
       this will be fixed in their next Python build.
