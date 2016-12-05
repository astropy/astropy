# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function, unicode_literals)

_tabversion   = '3.8'
_lextokens    = set(('CARET', 'CLOSE_PAREN', 'SOLIDUS', 'FUNCNAME', 'UFLOAT', 'STAR', 'SIGN', 'UINT', 'DOUBLE_STAR', 'OPEN_PAREN', 'PERIOD', 'UNIT'))
_lexreflags   = 32
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [("(?P<t_UFLOAT>((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_FUNCNAME>((sqrt)|(ln)|(exp)|(log)|(mag)|(dB)|(dex))(?=\\ *\\())|(?P<t_UNIT>%|([YZEPTGMkhdcmunpfazy]?'((?!\\d)\\w)+')|((?!\\d)\\w)+)|(?P<t_DOUBLE_STAR>\\*\\*)|(?P<t_CLOSE_PAREN>\\))|(?P<t_STAR>\\*)|(?P<t_PERIOD>\\.)|(?P<t_OPEN_PAREN>\\()|(?P<t_CARET>\\^)|(?P<t_SOLIDUS>/)", [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_FUNCNAME', 'FUNCNAME'), None, None, None, None, None, None, None, None, ('t_UNIT', 'UNIT'), None, None, None, (None, 'DOUBLE_STAR'), (None, 'CLOSE_PAREN'), (None, 'STAR'), (None, 'PERIOD'), (None, 'OPEN_PAREN'), (None, 'CARET'), (None, 'SOLIDUS')])]}
_lexstateignore = {'INITIAL': ' '}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
