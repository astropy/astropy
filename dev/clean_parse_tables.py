#!/usr/bin/env python

"""
A utility to fix PLY-generated lex and yacc tables to be Python 2 and
3 compatible.
"""

import os

for root, dirs, files in os.walk(os.path.join(os.path.dirname(__file__), '..')):
    for fname in files:
        if not (fname.endswith('lextab.py') or fname.endswith('parsetab.py')):
            continue

        path = os.path.join(root, fname)
        with open(path, 'rb') as fd:
            lines = fd.readlines()

        with open(path, 'wb') as fd:
            fd.write("# Licensed under a 3-clause BSD style license - see LICENSE.rst\n")
            fd.write("from __future__ import (absolute_import, division, print_function, unicode_literals)\n")
            fd.write('\n')

            lines = [x.replace("u'", "'").replace('u"', '"') for x in lines]
            lines = [x for x in lines if not (fname in x and x[0] == '#')]

            fd.write(''.join(lines))
