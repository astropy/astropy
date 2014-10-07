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
script as follows::

     >>> import warnings
     >>> from astropy.io import fits
     >>> warnings.filterwarnings('ignore', category=UserWarning, append=True)
     >>> fits.writeto(filename, data, clobber=True)

Astropy includes its own warning class,
`~astropy.utils.exceptions.AstropyUserWarning`, on which all warnings from
Astropy are based.  So one can also ignore warnings from Astropy (while still
allowing through warnings from other libraries like Numpy) by using something
like::

    >>> warnings.filterwarnings('ignore', category=AstropyUserWarning)

However, warning filters may also be modified just within a certain context
using the `~warnings.catch_warnings` context manager::

    >>> from warnings import catch_warnings
    >>> with catch_warnings():
    ...     warnings.filterwarnings('ignore', AstropyUserWarning)
    ...     fits.writeto(filename, data, clobber=True)

Astropy also issues warnings when deprecated API features are used.  If you
wish to *squelch* deprecation warnings, you can start Python with
``-Wi::Deprecation``.  This sets all deprecation warnings to ignored.  There is
also an Astropy-specific `~astropy.utils.exceptions.AstropyDeprecationWarning`
which can be used to disable deprecation warnings from Astropy only.

See http://docs.python.org/using/cmdline.html#cmdoption-unittest-discover-W for
more information on the -W argument.
