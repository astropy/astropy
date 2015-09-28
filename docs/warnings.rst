.. _python-warnings:

**********************
Python warnings system
**********************

.. doctest-skip-all

Astropy uses the Python :mod:`warnings` module to issue warning messages.  The
details of using the warnings module are general to Python, and apply to any
Python software that uses this system.  The user can suppress the warnings
using the python command line argument ``-W"ignore"`` when starting an
interactive python session.  For example::

     $ python -W"ignore"

The user may also use the command line argument when running a python script as
follows::

     $ python -W"ignore" myscript.py

It is also possible to suppress warnings from within a python script.  For
instance, the warnings issued from a single call to the
`astropy.io.fits.writeto` function may be suppressed from within a Python
script using the `warnings.filterwarnings` function as follows::

     >>> import warnings
     >>> from astropy.io import fits
     >>> warnings.filterwarnings('ignore', category=UserWarning, append=True)
     >>> fits.writeto(filename, data, clobber=True)

An equivalent way to insert an entry into the list of warning filter specifications
for simple call `warnings.simplefilter`::

    >>> warnings.simplefilter('ignore', UserWarning)

Astropy includes its own warning classes,
`~astropy.utils.exceptions.AstropyWarning` and
`~astropy.utils.exceptions.AstropyUserWarning`.  All warnings from Astropy are
based on these warning classes (see below for the distinction between them). One
can thus ignore all warnings from Astropy (while still allowing through
warnings from other libraries like Numpy) by using something like::

    >>> from astropy.utils.exceptions import AstropyWarning
    >>> warnings.simplefilter('ignore', category=AstropyWarning)

Warning filters may also be modified just within a certain context using the
`warnings.catch_warnings` context manager::

    >>> with warnings.catch_warnings():
    ...     warnings.simplefilter('ignore', AstropyWarning)
    ...     fits.writeto(filename, data, clobber=True)

As mentioned above, there are actually *two* base classes for Astropy warnings.
The main distinction is that `~astropy.utils.exceptions.AstropyUserWarning` is
for warnings that are *intended* for typical users (e.g. "Warning: Ambiguous
unit", something that might be because of improper input).  In contrast,
`~astropy.utils.exceptions.AstropyWarning` warnings that are *not*
`~astropy.utils.exceptions.AstropyUserWarning` may be for lower-level warnings
more useful for developers writing code that *uses* Astropy (e.g., the
deprecation warnings discussed below).  So if you're a user that just wants to
silence everything, the code above will suffice, but if you are a developer and
want to hide development-related warnings from your users, you may wish to still
allow through `~astropy.utils.exceptions.AstropyUserWarning`.

Astropy also issues warnings when deprecated API features are used.  If you
wish to *squelch* deprecation warnings, you can start Python with
``-Wi::Deprecation``.  This sets all deprecation warnings to ignored.  There is
also an Astropy-specific `~astropy.utils.exceptions.AstropyDeprecationWarning`
which can be used to disable deprecation warnings from Astropy only.

The base class of `~astropy.utils.exceptions.AstropyUserWarning` and
`~astropy.utils.exceptions.AstropyDeprecationWarning` is called simply
`~astropy.utils.exceptions.AstropyWarning`.  So *all* warnings coming from
Astropy can be be disabled with::

    >>> from astropy.utils.exceptions import AstropyWarning
    >>> warnings.simplefilter('ignore', category=AstropyWarning)

However, it is recommendable to use
`~astropy.utils.exceptions.AstropyUserWarning` in most cases, and to only
squelch deprecation warnings once you know they exist and have determined that
they can be ignored for the time being.

See `the CPython documentation
<http://docs.python.org/2/using/cmdline.html#cmdoption-W>`__ for more
information on the -W argument.
