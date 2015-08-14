# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function, unicode_literals)

_tabversion   = '3.5'
_lextokens    = set(['DIVISION', 'STAR', 'WHITESPACE', 'STARSTAR', 'UNKNOWN', 'LIT10', 'SIGN', 'OPEN_PAREN', 'UINT', 'UNIT', 'CLOSE_PAREN', 'UFLOAT'])
_lexreflags   = 0
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>(((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+))|(((\\d+\\.\\d*)|(\\.\\d+))([eE][+-]?\\d+)?))|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_X>[x\xd7])|(?P<t_LIT10>10)|(?P<t_UNKNOWN>[Uu][Nn][Kk][Nn][Oo][Ww][Nn])|(?P<t_UNIT>[a-zA-Z][a-zA-Z_]*)|(?P<t_WHITESPACE>[ \t]+)|(?P<t_STARSTAR>\\*\\*)|(?P<t_CLOSE_PAREN>\\))|(?P<t_OPEN_PAREN>\\()|(?P<t_STAR>\\*)|(?P<t_DIVISION>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, None, None, None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_X', 'X'), ('t_LIT10', 'LIT10'), ('t_UNKNOWN', 'UNKNOWN'), ('t_UNIT', 'UNIT'), (None, 'WHITESPACE'), (None, 'STARSTAR'), (None, 'CLOSE_PAREN'), (None, 'OPEN_PAREN'), (None, 'STAR'), (None, 'DIVISION')])]}
_lexstateignore = {'INITIAL': ''}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
