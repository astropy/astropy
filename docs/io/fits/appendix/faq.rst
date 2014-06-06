.. doctest-skip-all

PyFITS FAQ
----------

.. contents::

General Questions
^^^^^^^^^^^^^^^^^

What is PyFITS?
"""""""""""""""

PyFITS_ is a library written in, and for use with the Python_ programming
language for reading, writing, and manipulating FITS_ formatted files.  It
includes a high-level interface to FITS headers with the ability for high and
low-level manipulation of headers, and it supports reading image and table
data as Numpy_ arrays.  It also supports more obscure and non-standard formats
found in some FITS files.

PyFITS includes two command-line utilities for working with FITS files:
fitscheck, which can verify and write FITS checksums; and fitsdiff, which can
analyze and display the differences between two FITS files.

Although PyFITS is written mostly in Python, it includes an optional module
written in C that's required to read/write compressed image data.  However,
the rest of PyFITS functions without this extension module.

.. _PyFITS: http://www.stsci.edu/institute/software_hardware/pyfits
.. _Python: http://www.python.org
.. _FITS: http://fits.gsfc.nasa.gov/
.. _Numpy: http://numpy.scipy.org/

What is the development status of PyFITS?
"""""""""""""""""""""""""""""""""""""""""

PyFITS is written and maintained by the Science Software Branch at the `Space
Telescope Science Institute`_, and is licensed by AURA_ under a `3-clause BSD
license`_ (see `LICENSE.txt`_ in the PyFITS source code).

PyFITS' current primary developer and active maintainer is `Erik Bray`_, though
patch submissions are welcome from anyone.  It has a `Trac site`_ where the
source code can be browsed, and where bug reports may be submitted.  The source
code resides primarily in an `SVN repository`_ which allows anonymous
checkouts, though a `Git mirror`_ also exists.  PyFITS also has a `GitHub
site`_.

The current stable release series is 3.0.x.  Each 3.0.x release tries to
contain only bug fixes, and to not introduce any significant behavioral or API
changes (though this isn't guaranteed to be perfect).  The upcoming 3.1 release
will contain new features and some API changes, though will try maintain as
much backwards-compatibility as possible.  After the 3.1 release there may be
further 3.0.x releases for bug fixes only where possible.  Older versions of
PyFITS (2.4 and earlier) are no longer actively supported.

PyFITS is also included as a major component of upcoming Astropy_ project as
the :mod:`astropy.io.fits` module.  The goal is for Astropy to eventually serve
as a drop-in replacement for PyFITS. However, for the time being PyFITS will
still be released as an independent product as well, until such time that the
Astropy project proves successful and widely-adopted.

.. _Space Telescope Science Institute: http://www.stsci.edu/
.. _AURA: http://www.aura-astronomy.org/
.. _3-clause BSD license: http://en.wikipedia.org/wiki/BSD_licenses#3-clause_license_.28.22New_BSD_License.22_or_.22Modified_BSD_License.22.29
.. _LICENSE.txt: https://trac.assembla.com/pyfits/browser/trunk/LICENSE.txt
.. _Erik Bray: mailto:embray@stsci.edu
.. _Trac site: https://trac.assembla.com/pyfits/
.. _SVN repository: https://subversion.assembla.com/svn/pyfits/
.. _Git mirror: git://github.com/spacetelescope/PyFITS.git
.. _GitHub site: https://github.com/spacetelescope/PyFITS


Build and Installation Questions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Is PyFITS available on Windows?
"""""""""""""""""""""""""""""""

Yes--the majority of PyFITS is pure Python, and can be installed and used on
any platform that supports Python (>=2.5).  However, PyFITS includes an
optional C extension module for reading/writing compressed image HDUs.  As most
Windows systems are not configured to compile C source code, binary installers
are also available for Windows.  Though PyFITS can also be installed from
source even on Windows systems without a compiler by disabling the compression
module.  See `How do I install PyFITS from source on Windows?`_ for more
details.

Where is the Windows installer for version X of PyFITS on version Y of Python?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Every official PyFITS build for Windows is eventually uploaded to PyPI_.  This
includes builds for every major Python release from 2.5.x and up, except for
3.0 as there is no official Numpy release for Python 3.0 on Windows.  The one
binary module included in these builds was linked with Numpy 1.6.1, though it
should work with other recent Numpy versions.

Sometimes the Windows binary installers don't go up immediately after every
PyFITS release.  But if they appear missing they should go up within another
business day or two.  This has gotten better with recent releases thanks to
some automation.

.. _PyPI: http://pypi.python.org/pypi/pyfits

Why is the PyFITS installation failing on Windows?
""""""""""""""""""""""""""""""""""""""""""""""""""

The most likely cause of installation failure on Windows is if building/
installing from source fails due to the lack of a compiler for the optional C
extension module.  Such a failure would produce an error that looks something
like::

    building 'pyfits.compression' extension
    error: Unable to find vcvarsall.bat

Your best bet in cases like this is to install from one of the binary
executable installers available for Windows on PyPI.  However, there are still
cases where you may need to install from source: For example, it's difficult to
use the binary installers with virtualenv.  See `How do I install PyFITS from
source on Windows?`_ for more detailed instructions on building on Windows.

See below for a few answers to other specific questions about Windows
installation. For other installation errors not mentioned by this FAQ, please
contact help@stsci.edu with a description of the problem.

On Windows Vista or later why can't the installer find Python in the registry?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This is a common issue with Windows installers for Python packages that do not
support the new User Access Control (UAC) framework added in Windows Vista and
later.  In particular, when a Python is installed "for all users" (as opposed
to for a single user) it adds entries for that Python installation under the
``HKEY_LOCAL_MACHINE`` (HKLM) hierarchy and *not* under the
``HKEY_CURRENT_USER`` (HKCU) hierarchy.  However, depending on your UAC
settings, if the PyFITS installer is not executed with elevated privileges
it will not be able to check in HKLM for the required information about your
Python installation.

In short: If you encounter this problem it's because you need the appropriate
entries in the Windows registry for Python. You can download `this script`__
and execute it with the same Python as the one you want to install PyFITS
into.  For example to add the missing registry entries to your Python 2.7::

    C:\>C:\Python27\python.exe C:\Path\To\Downloads\win_register_python.py

__ https://gist.github.com/embray/6042780#file-win_register_python-py

How do I install PyFITS from source on Windows?
"""""""""""""""""""""""""""""""""""""""""""""""

There are a few options for building/installing PyFITS from source on Windows.

First of all, as mentioned elsewhere, most of PyFITS is pure-Python.  Only the
C extension module for reading/writing compressed images needs to be compiled.
If you don't need compressed image support, PyFITS can be installed without it.

In future releases this will hopefully be even easier, but for now it's
necessary to edit one file in order to disable the extension module.  Locate
the `setup.cfg`_ file at the root of the PyFITS source code.  This is the file
that describes what needs to be installed for PyFITS.  Find the line that reads
``[extension=pyfits.compression]``.  This is the section that lists what needs
to be compiled for the extension module.  Comment out every line in the
extension section by prepending it with a ``#`` character (stopping at the
``[build_ext]`` line).  It should look like this::

    ...
    scripts = scripts/fitscheck

    #[extension=pyfits.compression]
    #sources =
    #    src/compress.c
    #    src/fits_hcompress.c
    #    src/fits_hdecompress.c
    #    src/fitsio.c
    #    src/pliocomp.c
    #    src/compressionmodule.c
    #    src/quantize.c
    #    src/ricecomp.c
    #    src/zlib.c
    #    src/inffast.c
    #    src/inftrees.c
    #    src/trees.c
    #include_dirs = numpy
    # Squelch a handful of warnings (which actually cause pip to break in tox and
    # other environments due to gcc outputting non-ASCII characters in some
    # terminals; see python issue6135)
    #extra_compile_args =
    #    -Wno-unused-function
    #    -Wno-strict-prototypes

    [build_ext]
    ...

With these lines properly commented out, rerun ``python setup.py install``, and
it should skip building/installing the compression module.  PyFITS will work
fine with out it, but will issue warnings when encountering a compressed image
that it can't read.

If you do need to compile the compression module, this can still be done on
Windows with just a little extra work.  By default, Python tries to compile
extension modules with the same compiler that Python itself was compiled with.

To check what compiler Python was built with, the easiest way is to run::

    python -c "import platform; print platform.python_compiler()"

For the official builds of recent Python versions this should be something
like::

    MSC v.1500 32 bit (Intel)

For unofficial Windows distributions of Python, such as ActiveState, EPD, or
Cygwin, your mileage may vary.

As it so happens, MSC v.15xx is the compiler version included with Visual
C++ 2008.  Luckily, Microsoft distributes a free version of this as `Visual C++
Express Edition`_.  So for building Python extension modules on Windows this is
one of the simpler routes.  Just install the free VC++ 2008.  It should install
a link to the Start Menu at All Programs->Microsoft Visual C++ Express
Edition->Visual Studio Tools->Visual Studio 2008 Command Prompt.

If you run that link, it should launch a command prompt with reasonable
environment variables set up for using Visual C++.  Then change directories to
your copy of the PyFITS source code and re-run ``python setup.py install``.
You may also need to comment out the ``extra_compile_args`` option in the
``setup.cfg`` file (its value is the two lines under it after the equal sign).
Though the need to manually disable this option for MSC will be fixed in a
future PyFITS version.

Another option is to use gcc through `MinGW`_, which is in fact how the PyFITS
releases for Windows are currently built.  This article provides a good
overview of how to set this up: http://seewhatever.de/blog/?p=217

.. _setup.cfg: https://trac.assembla.com/pyfits/browser/trunk/setup.cfg
.. _Visual C++ Express Edition: http://www.microsoft.com/visualstudio/en-us/products/2008-editions/express
.. _MinGW: http://www.mingw.org/

Is PyFITS available for Mac OSX?
""""""""""""""""""""""""""""""""

Yes, but there is no binary package specifically for OSX (such as a .dmg, for
example).  For OSX just download, build, and install the source package.  This
is generally easier on OSX than it is on Windows, thanks to the more
developer-friendly environment.

The only major problem with building on OSX seems to occur for some users of
10.7 Lion, with misconfigured systems.  See the next question for details on
that.

To build PyFITS without the optional compression module, follow the
instructions in `How do I install PyFITS from source on Windows?`_.

Why is the PyFITS installation failing on OSX Lion (10.7)?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

There is a common problem that affects all Python packages with C extension
modules (not just PyFITS) for some users of OSX 10.7.  What usually occurs is
that when building the package several errors will be output, ending with
something like::

    unable to execute gcc-4.2: No such file or directory
    error: command 'gcc-4.2' failed with exit status 1

There are a few other errors like it that can occur.  The problem is that when
you build a C extension, by default it will try to use the same compiler that
your Python was built with. In this case, since you're using the 32-bit
universal build of Python it's trying to use the older gcc-4.2 and is trying
to build with PPC support, which is no longer supported in Xcode.

In this case the best solution is to install the x86-64 build of Python for
OSX (http://www.python.org/ftp/python/2.7.2/python-2.7.2-macosx10.6.dmg for
2.7.2).  In fact, this is the build version officially supported for use on
Lion.  Other, unofficial Python builds such as from `MacPorts`_ may also work.

.. _MacPorts: http://astrofrog.github.com/macports-python/

How do I find out what version of PyFITS I have installed?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To output the PyFITS version from the command line, run::

    $ python -c 'import pyfits; print(pyfits.__version__)'

When PyFITS is installed with stsci_python, it is also possible to check the
installed SVN revision by importing ``pyfits.svn_version``.  Then use
``dir(pyfits.svn_version)`` to view a list of available attributes.  A
feature like this will be available soon in standalone versions of PyFITS as
well.

How do I run the tests for PyFITS?
""""""""""""""""""""""""""""""""""

Currently the best way to run the PyFITS tests is to download the source code,
either from a source release or from version control, and to run the tests out
of the source.  It is not necessary to install PyFITS to run the tests out of
the source code.

The PyFITS tests require `nose`_ to run.  nose can be installed on any Python
version using pip or easy_install.  See the nose documentation for more
details.

With nose installed, it is simple to run the tests on Python 2.x::

    $ python setup.py nosetests

If PyFITS has not already been built, this will build it automatically, then
run the tests.  This does not cause PyFITS to be installed.

On Python 3.x the situation is a little more complicated.  This is due to the
fact that PyFITS' source code is not Python 3-compatible out of the box, but
has to be run through the 2to3 converter.  Normally when you build/install
PyFITS on Python 3.x, the 2to3 conversion is performed automatically.
Unfortunately, nose does not know to use the 2to3'd source code, and will
instead try to import and test the unconverted source code.

To work around this, it is necessary to first build PyFITS (which will run the
source through 2to3)::

    $ python setup.py build

Then run the ``nosetests`` command, but pointing it to the ``build`` tree
where the 2to3'd source code and tests reside, using the ``-w`` switch::

    $ python setup.py nosetests -w build/lib.linux-x86_64-3.2

where the exact path of the ``build/lib.*`` directory will vary depending on
your platform and Python version.

.. _nose: http://readthedocs.org/docs/nose/en/latest/

How can I build a copy of the PyFITS documentation for my own use?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

First of all, it's worth pointing out that the documentation for the latest
version of PyFITS can always be downloaded in `PDF form
<http://stsdas.stsci.edu/download/docs/The_PyFITS_Handbook.pdf>`_ or browsed
online in `HTML <http://packages.python.org/pyfits>`_.  There are also plans
to make the docs for older versions of PyFITS, as well as up-to-date
development docs available online.

Otherwise, to build your own version of the docs either for offline use, or to
build the development version of the docs there are a few requirements.  The
most import requirement is `Sphinx`_, which is the toolkit used to generate
the documentation.  Use ``pip install sphinx`` or ``easy_install sphinx`` to
install Sphinx.  Using pip or easy_install will install the correct versions
of Sphinx's basic dependencies, which include docutils, Jinja2, and Pygments.

Next, the docs require STScI's custom Sphinx theme, `stsci.sphinxext`_.  It's
a simple pure-Python package and can be installed with pip or easy_install.

The next requirement is `numpydoc`_, which is not normally installed with
Numpy itself.  Install it with pip or easy_install.  Numpy is also required,
though it is of course a requirement of PyFITS itself.

Finally, it is necessary to have `matplotlib`_, specifically for
matplotlib.sphinxext.  This is perhaps the most onerous requirement if you do
not already have it installed. Please refer to the matplotlib documentation for
details on downloading and installing matplotlib.

It is also necessary to install PyFITS itself in order to generate the API
documentation.  For this reason, it is a good idea to install Sphinx and
PyFITS into a `virtualenv`_ in order to build the development version of the
docs (see below).

With all the requirements installed, change directories into the ``docs/``
directory in the PyFITS source code, and run::

    $ make html

to build the HTML docs, which will be output to ``build/html``.  To build the
docs in other formats, please refer to the Sphinx documentation.

To summarize, assuming that you already have Numpy and Matplotlib on your
Python installation, perform the following steps from within the PyFITS source
code::

    $ virtualenv --system-site-packages pyfits-docs
    $ source pyfits-docs/bin/activate
    $ pip install sphinx
    $ pip install numpydoc
    $ pip install stsci.sphinxext
    $ python setup.py install pyfits
    $ cd docs/
    $ make html


.. _Sphinx: http://sphinx.pocoo.org/
.. _stsci.sphinxext: http://pypi.python.org/pypi/stsci.sphinxext
.. _numpydoc: http://pypi.python.org/pypi/numpydoc
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _virtualenv: http://pypi.python.org/pypi/virtualenv


Usage Questions
^^^^^^^^^^^^^^^

Something didn't work as I expected.  Did I do something wrong?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Possibly.  But if you followed the documentation and things still did not work
as expected, it is entirely possible that there is a mistake in the
documentation, a bug in the code, or both.  So feel free to report it as a
bug.  There are also many, many corner cases in FITS files, with new ones
discovered almost every week.  PyFITS is always improving, but does not
support all cases perfectly.  There are some features of the FITS format
(scaled data, for example) that are difficult to support correctly and can
sometimes cause unexpected behavior.

For the most common cases, however, such as reading and updating FITS headers,
images, and tables, PyFITS should be very stable and well-tested.  Before
every PyFITS release it is ensured that all its tests pass on a variety of
platforms, and those tests cover the majority of use-cases (until new
corner cases are discovered).

PyFITS crashed and output a long string of code.  What do I do?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This listing of code is what is knows as a `stack trace`_ (or in Python
parlance a "traceback").  When an unhandled exception occurs in the code,
causing the program to end, this is a way of displaying where the exception
occurred and the path through the code that led to it.

As PyFITS is meant to be used as a piece in other software projects, some
exceptions raised by PyFITS are by design.  For example, one of the most
common exceptions is a `~.exceptions.KeyError` when an attempt is made to
read the value of a non-existent keyword in a header::

    >>> import pyfits
    >>> h = pyfits.Header()
    >>> h['NAXIS']
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/path/to/pyfits/header.py", line 125, in __getitem__
        return self._cards[self._cardindex(key)].value
      File "/path/to/pyfits/header.py", line 1535, in _cardindex
        raise KeyError("Keyword %r not found." % keyword)
    KeyError: "Keyword 'NAXIS' not found."

This indicates that something was looking for a keyword called "NAXIS" that
does not exist.  If an error like this occurs in some other software that uses
PyFITS, it may indicate a bug in that software, in that it expected to find a
keyword that didn't exist in a file.

Most "expected" exceptions will output a message at the end of the traceback
giving some idea of why the exception occurred and what to do about it.  The
more vague and mysterious the error message in an exception appears, the more
likely that it was caused by a bug in PyFITS.  So if you're getting an
exception and you really don't know why or what to do about it, feel free to
report it as a bug.

.. _stack trace: http://en.wikipedia.org/wiki/Stack_trace

Why does opening a file work in CFITSIO, ds9, etc. but not in PyFITS?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

As mentioned elsewhere in this FAQ, there are many unusual corner cases when
dealing with FITS files.  It's possible that a file should work, but isn't
support due to a bug.  Sometimes it's even possible for a file to work in an
older version of PyFITS, but not a newer version due to a regression that
isn't tested for yet.

Another problem with the FITS format is that, as old as it is, there are many
conventions that appear in files from certain sources that do not meet the
FITS standard.  And yet they are so common-place that it is necessary to
support them in any FITS readers.  CONTINUE cards are one such example.  There
are non-standard conventions supported by PyFITS that are not supported by
CFITSIO and vice-versa.  You may have hit one of those cases.

If PyFITS is having trouble opening a file, a good way to rule out whether not
the problem is with PyFITS is to run the file through the `fitsverify`_.  For
smaller files you can even use the `online FITS verifier`_.  These use CFITSIO
under the hood, and should give a good indication of whether or not there is
something erroneous about the file.  If the file is malformatted, fitsverify
will output errors and warnings.

If fitsverify confirms no problems with a file, and PyFITS is still having
trouble opening it (especially if it produces a traceback) then it's likely
there is a bug in PyFITS.

.. _fitsverify: http://heasarc.gsfc.nasa.gov/docs/software/ftools/fitsverify/
.. _online FITS verifier: http://fits.gsfc.nasa.gov/fits_verify.html

How do I turn off the warning messages PyFITS keeps outputting to my console?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

PyFITS uses Python's built-in `warnings`_ subsystem for informing about
exceptional conditions in the code that are recoverable, but that the user may
want to be informed of.  One of the most common warnings in PyFITS occurs when
updating a header value in such a way that the comment must be truncated to
preserve space::

    Card is too long, comment is truncated.

Any console output generated by PyFITS can be assumed to be from the warnings
subsystem.  Fortunately there are two easy ways to quiet these warnings:

 1. Using the `-W option`_ to the ``python`` executable.  Just start Python
    like::

        $ python -Wignore <scriptname>

    or for short::

        $ python -Wi <scriptname>

    and all warning output will be silenced.

 2. Warnings can be silenced programatically from anywhere within a script.
    For example, to disable all warnings in a script, add something like::

        import warnings
        warnings.filterwarnings('ignore')

 Another option, instead of ``ignore`` is ``once``, which causes any warning
 to be output only once within the session, rather than repeatedly (such as in
 a loop).  There are many more ways to filter warnings with ``-W`` and the
 warnings module.  For example, it is possible to silence only specific
 warning messages.  Please refer to the Python documentation for more details,
 or ask at help@stsci.edu.

.. _warnings: http://docs.python.org/library/warnings.html
.. _-W option: http://docs.python.org/using/cmdline.html#cmdoption-W

How can I check if my code is using deprecated PyFITS features?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

PyFITS 3.0 included a major reworking of the code and some of the APIs.  Most
of the differences are just renaming functions to use a more consistent naming
scheme.  For example the ``createCard()`` function was renamed to
``create_card()`` for consistency with a ``lower_case_underscore`` naming
scheme for functions.

There are a few other functions and attributes that were deprecated either
because they were renamed to something simpler or more consistent, or because
they were redundant or replaced.

Eventually all deprecated features will be removed in future PyFITS versions
(though there will be significant warnings in advance).  It is best to check
whether your code is using deprecated features sooner rather than later.

On Python 2.5, all deprecation warnings are displayed by default, so you may
have already discovered them.  However, on Python 2.6 and up, deprecation
warnings are *not* displayed by default.  To show all deprecation warnings,
start Python like::

    $ python -Wd <scriptname>

Most deprecation issues can be fixed with a simple find/replace.  The warnings
displayed will let you know how to replace the old interface.

If you have a lot of old code that was written for older versions of PyFITS it
would be worth doing this.  PyFITS 3.1 introduces a significant rewrite of the
Header interface, and contains even more deprecations.

What convention does PyFITS use for indexing, such as of image coordinates?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

All arrays and sequences in PyFITS use a zero-based indexing scheme.  For
example, the first keyword in a header is ``header[0]``, not ``header[1]``.
This is in accordance with Python itself, as well as C, on which PyFITS is
based.

This may come as a surprise to veteran FITS users coming from IRAF, where
1-based indexing is typically used, due to its origins in FORTRAN.

Likewise, the top-left pixel in an N x N array is ``data[0,0]``.  The indices
for 2-dimensional arrays are row-major order, in that the first index is the
row number, and the second index is the column number.  Or put in terms of
axes, the first axis is the y-axis, and the second axis is the x-axis.  This is
the opposite of column-major order, which is used by FORTRAN and hence FITS.
For example, the second index refers to the axis specified by NAXIS1 in the
FITS header.

In general, for N-dimensional arrays, row-major orders means that the
right-most axis is the one that varies the fastest while moving over the
array data linearly.  For example, the 3-dimensional array::

    [[[1, 2],
      [3, 4]],
     [[5, 6],
      [7, 8]]]

is represented linearly in row-major order as::

    [1, 2, 3, 4, 5, 6, 7, 8]

Since 2 immediately follows 1, you can see that the right-most (or inner-most)
axis is the one that varies the fastest.

The discrepancy in axis-ordering may take some getting used to, but it is a
necessary evil.  Since most other Python and C software assumes row-major
ordering, trying to enforce column-major ordering in arrays returned by PyFITS
is likely to cause more difficulties than it's worth.

How do I open a very large image that won't fit in memory?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Prior to PyFITS 3.1, when the data portion of an HDU is accessed, the data is
read into memory in its entirety.  For example::

    >>> hdul = pyfits.open('myimage.fits')
    >>> hdul[0].data
    ...

reads the entire image array from disk into memory.  For very large images or
tables this is clearly undesirable, if not impossible given the available
resources.

However, ``pyfits.open()`` has an option to access the data portion of an HDU
by memory mapping using `mmap`_.  What this means is that accessing the data
as in the example above only reads portions of the data into memory on demand.
For example, if I request just a slice of the image, such as
``hdul[0].data[100:200]``, then just rows 100-200 will be read into memory.
This happens transparently, as though the entire image were already in memory.
This works the same way for tables.  For most cases this is your best bet for
working with large files.

To use memory mapping, just add the ``memmap=True`` argument to
``pyfits.open()``.

In PyFITS 3.1, the mmap support is improved enough that ``memmap=True`` is the
default for all ``pyfits.open()`` calls.  The default can also be controlled
through an environment variable called ``PYFITS_USE_MEMMAP``.  Setting this to
``0`` will disable mmap by default.

Unfortunately, memory mapping does not currently work as well with scaled
image data, where BSCALE and BZERO factors need to be applied to the data to
yield physical values.  Currently this requires enough memory to hold the
entire array, though this is an area that will see improvement in the future.

An alternative, which currently only works for image data (that is,
non-tables) is the sections interface.  It is largely replaced by the better
support for memmap, but may still be useful on systems with more limited
virtual-memory space, such as on 32-bit systems.  Support for scaled image
data is flakey with sections too, though that will be fixed.  See `the PyFITS
documentation
<http://packages.python.org/pyfits/users_guide/users_image.html#data-section>`_
for more details on working with sections.

.. _mmap: http://en.wikipedia.org/wiki/Mmap

How can I create a very large FITS file from scratch?
"""""""""""""""""""""""""""""""""""""""""""""""""""""

This is a very common issue, but unfortunately PyFITS does not come with any
built-in facilities for creating large files (larger than will fit in memory)
from scratch (though it may in the future).

Normally to create a single image FITS file one would do something like::

    >> data = numpy.zeros((40000, 40000), dtype=numpy.float64)
    >> hdu = pyfits.PrimaryHDU(data=data)
    >> hdu.writeto('large.fits')

However, a 40000 x 40000 array of doubles is nearly twelve gigabytes!  Most
systems won't be able to create that in memory just to write out to disk.  In
order to create such a large file efficiently requires a little extra work,
and a few assumptions.

First, it is helpful to anticipate about how large (as in, how many keywords)
the header will have in it.  FITS headers must be written in 2880 byte
blocks--large enough for 36 keywords per block (including the END keyword in
the final block).  Typical headers have somewhere between 1 and 4 blocks,
though sometimes more.

Since the first thing we write to a FITS file is the header, we want to write
enough header blocks so that there is plenty of padding in which to add new
keywords without having to resize the whole file.  Say you want the header to
use 4 blocks by default.  Then, excluding the END card which PyFITS will add
automatically, create the header and pad it out to 36 * 4 cards like so::

    >>> data = numpy.zeros((100, 100), dtype=numpy.float64)
    # This is a stub array that we'll be using the initialize the HDU; its
    # exact size is irrelevant, as long as it has the desired number of
    # dimensions
    >>> hdu = pyfits.PrimaryHDU(data=data)
    >>> header = hdu.header
    >>> while len(header) < (36 * 4 - 1):
    ...     header.append()  # Adds a blank card to the end

Now adjust the NAXISn keywords to the desired size of the array, and write
*only* the header out to a file.  Using the ``hdu.writeto()`` method will
cause PyFITS to "helpfully" reset the NAXISn keywords to match the size of the
dummy array::

    >>> header['NAXIS1'] = 40000
    >>> header['NAXIS2'] = 40000
    >>> header.tofile('large.fits')

Finally, we need to grow out the end of the file to match the length of the
data (plus the length of the header).  This can be done very efficiently on
most systems by seeking past the end of the file and writing a single byte,
like so::

    >>> with open('large.fits', 'rb+') as fobj:
    ...     fobj.seek(len(header.tostring()) + (40000 * 40000 * 8) - 1)
    ...     # The -1 is to account for the final byte that we are about to
    ...     # write
    ...     fobj.write('\0')

On modern operating systems this will cause the file (past the header) to be
filled with zeros out to the ~12GB needed to hold a 40000 x 40000 image.  On
filesystems that support sparse file creation (most Linux filesystems, but not
HFS+) this is a very fast, efficient operation.  On other systems your mileage
may vary.

This isn't the only way to build up a large file, but probably one of the
safest.  This method can also be used to create large multi-extension FITS
files, with a little care.

For creating very large tables, this method may also be used.  Though it can
be difficult to determine ahead of time how many rows a table will need.  In
general, use of PyFITS is discouraged for the creation and manipulation of
large tables.  The FITS format itself is not designed for efficient on-disk or
in-memory manipulation of table structures.  For large, heavy-duty table data
it might be better too look into using `HDF5`_ through the `PyTables`_
library.

PyTables makes use of Numpy under the hood, and can be used to write binary
table data to disk in the same format required by FITS.  It is then possible
to serialize your table to the FITS format for distribution.  At some point
this FAQ might provide an example of how to do this.

.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _PyTables: http://www.pytables.org/moin

How do I create a multi-extension FITS file from scratch?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When you open a FITS file with ``pyfits.open()``, a ``pyfits.HDUList`` object
is returned, which holds all the HDUs in the file.  This ``HDUList`` class is
a subclass of Python's builtin ``list``, and can be created from scratch and
used as such::

    >>> new_hdul = pyfits.HDUList()
    >>> new_hdul.append(pyfits.ImageHDU())
    >>> new_hdul.append(pyfits.ImageHDU())
    >>> new_hdul.writeto('test.fits')

That will create a new multi-extension FITS file with two empty IMAGE
extensions (a default PRIMARY HDU is prepended automatically if one was not
provided manually).

Why is an image containing integer data being converted unexpectedly to floats?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If the header for your image contains non-trivial values for the optional
BSCALE and/or BZERO keywords (that is, BSCALE != 1 and/or BZERO != 0), then
the raw data in the file must be rescaled to its physical values according to
the formula::

    physical_value = BZERO + BSCALE * array_value

As BZERO and BSCALE are floating point values, the resulting value must be a
float as well.  If the original values were 16-bit integers, the resulting
values are single-precision (32-bit) floats.  If the original values were
32-bit integers the resulting values are double-precision (64-bit floats).

This automatic scaling can easily catch you of guard if you're not expecting
it, because it doesn't happen until the data portion of the HDU is accessed
(to allow things like updating the header without rescaling the data).  For
example::

    >>> hdul = pyfits.open('scaled.fits')
    >>> image = hdul['SCI', 1]
    >>> image.header['BITPIX']
    32
    >>> image.header['BSCALE']
    2.0
    >>> data = image.data  # Read the data into memory
    >>> data.dtype
    dtype('float64')  # Got float64 despite BITPIX = 32 (32-bit int)
    >>> image.header['BITPIX']  # The BITPIX will automatically update too
    -64
    >>> 'BSCALE' in image.header  # And the BSCALE keyword removed
    False

The reason for this is that once a user accesses the data they may also
manipulate it and perform calculations on it.  If the data were forced to
remain as integers, a great deal of precision is lost.  So it is best to err
on the side of not losing data, at the cost of causing some confusion at
first.

If the data must be returned to integers before saving, use the ``scale()``
method::

    >>> image.scale('int32')
    >>> image.header['BITPIX']
    32

Alternatively, if a file is opened with ``mode='update'`` along with the
``scale_back=True`` argument, the original BSCALE and BZERO scaling will
be automatically re-applied to the data before saving.  Usually this is
not desireable, especially when converting from floating point back to
unsigned integer values.  But this may be useful in cases where the raw
data needs to be modified corresponding to changes in the physical values.

To prevent rescaling from occurring at all (good for updating headers--even if
you don't intend for the code to access the data, it's good to err on the side
of caution here), use the ``do_not_scale_image_data`` argument when opening
the file::

    >>> hdul = pyfits.open('scaled.fits', do_not_scale_image_data=True)
    >>> image = hdul['SCI', 1]
    >>> image.data.dtype
    dtype('int32')

Why am I losing precision when I assign floating point values in the header?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The FITS standard allows two formats for storing floating-point numbers in a
header value.  The "fixed" format requires the ASCII representation of the
number to be in bytes 11 through 30 of the header card, and to be
right-justified.  This leaves a standard number of characters for any comment
string.

The fixed format is not wide enough to represent the full range of values that
can be stored in a 64-bit float with full precision.  So FITS also supports a
"free" format in which the ASCII representation can be stored anywhere, using
the full 70 bytes of the card (after the keyword).

Currently PyFITS only supports writing fixed format (it can read both
formats), so all floating point values assigned to a header are stored in the
fixed format.  There are plans to add support for more flexible formatting.

In the meantime it is possible to add or update cards by manually formatting
the card image::

    >>> c = pyfits.Card.fromstring('FOO     = 1234567890.123456789')
    >>> h = pyfits.Header()
    >>> h.append(c)
    >>> h
    FOO     = 1234567890.123456789

As long as you don't assign new values to 'FOO' via ``h['FOO'] = 123``, PyFITS
will maintain the header value exactly as you formatted it.
