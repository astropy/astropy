# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function, unicode_literals)

_tabversion   = '3.5'
_lextokens    = set(['DIVISION', 'PRODUCT', 'SIGN', 'OPEN_PAREN', 'UINT', 'UNIT', 'X', 'CLOSE_PAREN', 'UFLOAT'])
_lexreflags   = 32
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>((\\d+\\.?\\d+)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_X>[x\xd7])|(?P<t_UNIT>\\%|\xb0|\\\\h|((?!\\d)\\w)+)|(?P<t_CLOSE_PAREN>\\))|(?P<t_OPEN_PAREN>\\()|(?P<t_PRODUCT>\\.)|(?P<t_DIVISION>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_X', 'X'), ('t_UNIT', 'UNIT'), None, (None, 'CLOSE_PAREN'), (None, 'OPEN_PAREN'), (None, 'PRODUCT'), (None, 'DIVISION')])]}
_lexstateignore = {'INITIAL': ''}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
