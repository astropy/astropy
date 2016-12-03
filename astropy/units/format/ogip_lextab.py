# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function, unicode_literals)

_tabversion   = '3.8'
_lextokens    = set(('UINT', 'WHITESPACE', 'STAR', 'OPEN_PAREN', 'UNIT', 'SIGN', 'UFLOAT', 'STARSTAR', 'DIVISION', 'LIT10', 'CLOSE_PAREN', 'UNKNOWN'))
_lexreflags   = 0
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>(((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+))|(((\\d+\\.\\d*)|(\\.\\d+))([eE][+-]?\\d+)?))|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_X>[x×])|(?P<t_LIT10>10)|(?P<t_UNKNOWN>[Uu][Nn][Kk][Nn][Oo][Ww][Nn])|(?P<t_UNIT>[a-zA-Z][a-zA-Z_]*)|(?P<t_WHITESPACE>[ \t]+)|(?P<t_STARSTAR>\\*\\*)|(?P<t_CLOSE_PAREN>\\))|(?P<t_OPEN_PAREN>\\()|(?P<t_STAR>\\*)|(?P<t_DIVISION>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, None, None, None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_X', 'X'), ('t_LIT10', 'LIT10'), ('t_UNKNOWN', 'UNKNOWN'), ('t_UNIT', 'UNIT'), (None, 'WHITESPACE'), (None, 'STARSTAR'), (None, 'CLOSE_PAREN'), (None, 'OPEN_PAREN'), (None, 'STAR'), (None, 'DIVISION')])]}
_lexstateignore = {'INITIAL': ''}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
