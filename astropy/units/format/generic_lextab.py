from __future__ import absolute_import, division, print_function, unicode_literals

# generic_lextab.py. This file automatically created by PLY (version 3.4). Don't edit!
_tabversion   = '3.4'
_lextokens    = {'CARET': 1, 'SOLIDUS': 1, 'STAR': 1, 'DOUBLE_STAR': 1, 'PERIOD': 1, 'SQRT': 1, 'SIGN': 1, 'OPEN_PAREN': 1, 'UINT': 1, 'UNIT': 1, 'CLOSE_PAREN': 1, 'UFLOAT': 1}
_lexreflags   = 0
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_SQRT>sqrt)|(?P<t_UNIT>%|[a-zA-Z][a-zA-Z_]*)|(?P<t_DOUBLE_STAR>\\*\\*)|(?P<t_CLOSE_PAREN>\\))|(?P<t_CARET>\\^)|(?P<t_PERIOD>\\.)|(?P<t_STAR>\\*)|(?P<t_OPEN_PAREN>\\()|(?P<t_SOLIDUS>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_SQRT', 'SQRT'), ('t_UNIT', 'UNIT'), (None, 'DOUBLE_STAR'), (None, 'CLOSE_PAREN'), (None, 'CARET'), (None, 'PERIOD'), (None, 'STAR'), (None, 'OPEN_PAREN'), (None, 'SOLIDUS')])]}
_lexstateignore = {'INITIAL': ' '}
_lexstateerrorf = {'INITIAL': 't_error'}
