# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/units
#
# You can then commit the changes to this file.

# cds_lextab.py. This file automatically created by PLY (version 3.11). Don't edit!
_tabversion   = '3.10'
_lextokens    = set(('CLOSE_BRACKET', 'CLOSE_PAREN', 'DIVISION', 'OPEN_BRACKET', 'OPEN_PAREN', 'PRODUCT', 'SIGN', 'UFLOAT', 'UINT', 'UNIT', 'X'))
_lexreflags   = 32
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>((\\d+\\.?\\d+)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_X>[x×])|(?P<t_UNIT>\\%|°|\\\\h|((?!\\d)\\w)+)|(?P<t_PRODUCT>\\.)|(?P<t_OPEN_PAREN>\\()|(?P<t_CLOSE_PAREN>\\))|(?P<t_OPEN_BRACKET>\\[)|(?P<t_CLOSE_BRACKET>\\])|(?P<t_DIVISION>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_X', 'X'), ('t_UNIT', 'UNIT'), None, (None, 'PRODUCT'), (None, 'OPEN_PAREN'), (None, 'CLOSE_PAREN'), (None, 'OPEN_BRACKET'), (None, 'CLOSE_BRACKET'), (None, 'DIVISION')])]}
_lexstateignore = {'INITIAL': ''}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
