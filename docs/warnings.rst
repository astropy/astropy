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

Astropy includes its own warning class,
`~astropy.utils.exceptions.AstropyUserWarning`, on which all warnings from
Astropy are based.  So one can also ignore warnings from Astropy (while still
allowing through warnings from other libraries like Numpy) by using something
like::

    >>> from astropy.utils.exceptions import AstropyUserWarning
    >>> warnings.simplefilter('ignore', category=AstropyUserWarning)

However, warning filters may also be modified just within a certain context
using the `warnings.catch_warnings` context manager::

    >>> with warnings.catch_warnings():
    ...     warnings.simplefilter('ignore', AstropyUserWarning)
    ...     fits.writeto(filename, data, clobber=True)

Astropy also issues warnings when deprecated API features are used.  If you
wish to *squelch* deprecation warnings, you can start Python with
``-Wi::Deprecation``.  This sets all deprecation warnings to ignored.  There is
also an Astropy-specific `~astropy.utils.exceptions.AstropyDeprecationWarning`
which can be used to disable deprecation warnings from Astropy only.

See `the CPython documentation
<http://docs.python.org/2/using/cmdline.html#cmdoption-W>`__ for more
information on the -W argument.
