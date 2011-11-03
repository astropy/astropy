**********************
Miscellaneous Features
**********************

In this chapter, we'll describe some of the miscellaneous features of PyFITS.


Warning Messages
================

PyFITS uses the Python warnings module to issue warning messages.  The user can
suppress the warnings using the python command line argument ``-W"ignore"``
when starting an interactive python session.  For example:

.. parsed-literal::

     python -W"ignore"

The user may also use the command line argument when running a python script as
follows:

.. parsed-literal::

     python -W"ignore" myscript.py

It is also possible to suppress warnings from within a python script.  For
instance, the warnings issued from a single call to the writeto convenience
function may be suppressed from within a python script as follows:

.. parsed-literal::

     import warnings
     import pyfits

     # ...

     warnings.resetwarnings()
     warnings.filterwarnings('ignore', category=UserWarning, append=True)
     pyfits.writeto(file, im, clobber=True)
     warnings.resetwarnings()
     warnings.filterwarnings('always', category=UserWarning, append=True)

     # ...

PyFITS also issues warnings when deprecated API features are used.  In Python
2.7 and up deprecation warnings are ignored by default.  To run Python with
deprecation warnings enabled, either start Python with the ``-Wall`` argument,
or you can enable deprecation warnings specifically with ``-Wd``.

In Python versions below 2.7, if you wish to *squelch* deprecation warnings,
you can start Python with ``-Wi::Deprecation``.  This sets all deprecation
warnings to ignored.  See
http://docs.python.org/using/cmdline.html#cmdoption-unittest-discover-W
for more information on the -W argument.
