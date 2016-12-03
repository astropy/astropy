# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function, unicode_literals)

_tabversion   = '3.8'
_lextokens    = set(('UINT', 'PRODUCT', 'OPEN_PAREN', 'UNIT', 'SIGN', 'UFLOAT', 'DIVISION', 'CLOSE_PAREN', 'X'))
_lexreflags   = 32
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>((\\d+\\.?\\d+)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_X>[x×])|(?P<t_UNIT>\\%|°|\\\\h|((?!\\d)\\w)+)|(?P<t_CLOSE_PAREN>\\))|(?P<t_OPEN_PAREN>\\()|(?P<t_PRODUCT>\\.)|(?P<t_DIVISION>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_X', 'X'), ('t_UNIT', 'UNIT'), None, (None, 'CLOSE_PAREN'), (None, 'OPEN_PAREN'), (None, 'PRODUCT'), (None, 'DIVISION')])]}
_lexstateignore = {'INITIAL': ''}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
