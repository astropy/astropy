**********************
Python warnings system
**********************

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
instance, the warnings issued from a single call to the writeto convenience
function may be suppressed from within a python script as follows::

     >>> import warnings
     >>> from astropy.io import fits
     >>> warnings.filterwarnings('ignore', category=UserWarning, append=True)
     >>> fits.writeto(file, im, clobber=True)  # doctest: +SKIP

However, warning filters may also be modified just within a certain context
using `warnings.catch_warnings`.

Astropy also issues warnings when deprecated API features are used.  In Python
2.7 and up deprecation warnings are ignored by default.  To run Python with
deprecation warnings enabled, either start Python with the ``-Wall`` argument,
or you can enable deprecation warnings specifically with ``-Wd``.

In Python versions below 2.7, if you wish to *squelch* deprecation warnings,
you can start Python with ``-Wi::Deprecation``.  This sets all deprecation
warnings to ignored.  See
http://docs.python.org/using/cmdline.html#cmdoption-unittest-discover-W
for more information on the -W argument.
