# STDLIB
import io
import os
import re

# LOCAL
from .. import exceptions
from astropy import setup_helpers

warnings_header = u"""
.. _warnings:

Warnings
========

.. note::
    Most of the following warnings indicate violations of the VOTable
    specification.  They should be reported to the authors of the
    tools that produced the VOTable file.

    To control the warnings emitted, use the standard Python
    :mod:`warnings` module.  Most of these are of the type
    `VOTableSpecWarning`.

.. contents::

"""


exceptions_header = u"""
.. _exceptions:

Exceptions
==========

.. note::

    This is a list of many of the fatal exceptions emitted by vo.table
    when the file does not conform to spec.  Other exceptions may be
    raised due to unforeseen cases or bugs in vo.table itself.

.. contents::

"""


_find_dedent_regex = re.compile("(?:(?:\n\r?)|^)( *)\S")
_dedent_regex = {}
def dedent(s):
    if not s:      # includes case of s is None
        return u''

    if not isinstance(s, unicode):
        s = s.decode('utf-8')

    match = _find_dedent_regex.match(s)
    if match is None:
        return s

    nshift = match.end(1) - match.start(1)
    if nshift == 0:
        return s

    unindent = _dedent_regex.get(nshift, None)
    if unindent is None:
        unindent = re.compile("\n\r? {0,%d}" % nshift)
        _dedent_regex[nshift] = unindent

    result = unindent.sub("\n", s).strip()
    return result


def generate_set(filename, prefix, header):
    classes = []
    for key, val in exceptions.__dict__.items():
        if re.match(prefix + "[0-9]{2}", key):
            classes.append((key, val))
    classes.sort()

    out = io.StringIO()

    out.write(header)
    out.write(u'\n')
    for name, cls in classes:
        out.write(u".. _%s:\n\n" % name)
        msg = "%s: %s" % (cls.__name__, cls.get_short_name())
        if not isinstance(msg, unicode):
            msg = msg.decode('utf-8')
        out.write(msg)
        out.write(u'\n')
        out.write(u'-' * len(msg))
        out.write(u'\n\n')
        out.write(dedent(cls.__doc__))
        out.write(u'\n\n')

    setup_helpers.write_if_different(
        filename, out.getvalue().encode('utf-8'))


def setup(app):
    root = os.path.join(app.srcdir, 'vo')

    generate_set(os.path.join(root, 'warnings.rst'), 'W', warnings_header)
    generate_set(os.path.join(root, 'exceptions.rst'), 'E', exceptions_header)
