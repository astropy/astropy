============
Known Issues
============

While most bugs and issues are managed using the `astropy issue
tracker <https://github.com/astropy/astropy/issues>`_, this document
lists issues that are too difficult to fix, may require some
intervention from the user to workaround, or require upstream fixes in
other projects.

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

The IPython issue: `https://github.com/ipython/ipython/pull/2738`_

Floating point precision issues on Python 2.6 on Microsoft Windows
------------------------------------------------------------------

When converting floating point numbers to strings on Python 2.6 on a
Microsoft Windows platform, some of the requested precision may be
lost.

The easiest workaround is to install Python 2.7.

The Python issue: `http://bugs.python.org/issue7117`_
