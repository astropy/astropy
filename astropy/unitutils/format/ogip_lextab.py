# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/units
#
# You can then commit the changes to this file.

# ogip_lextab.py. This file automatically created by PLY (version 3.10). Don't edit!
_tabversion   = '3.10'
_lextokens    = set(('CLOSE_PAREN', 'DIVISION', 'LIT10', 'OPEN_PAREN', 'SIGN', 'STAR', 'STARSTAR', 'UFLOAT', 'UINT', 'UNIT', 'UNKNOWN', 'WHITESPACE'))
_lexreflags   = 64
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>(((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+))|(((\\d+\\.\\d*)|(\\.\\d+))([eE][+-]?\\d+)?))|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_X>[x×])|(?P<t_LIT10>10)|(?P<t_UNKNOWN>[Uu][Nn][Kk][Nn][Oo][Ww][Nn])|(?P<t_UNIT>[a-zA-Z][a-zA-Z_]*)|(?P<t_WHITESPACE>[ \t]+)|(?P<t_STARSTAR>\\*\\*)|(?P<t_STAR>\\*)|(?P<t_CLOSE_PAREN>\\))|(?P<t_OPEN_PAREN>\\()|(?P<t_DIVISION>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, None, None, None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_X', 'X'), ('t_LIT10', 'LIT10'), ('t_UNKNOWN', 'UNKNOWN'), ('t_UNIT', 'UNIT'), (None, 'WHITESPACE'), (None, 'STARSTAR'), (None, 'STAR'), (None, 'CLOSE_PAREN'), (None, 'OPEN_PAREN'), (None, 'DIVISION')])]}
_lexstateignore = {'INITIAL': ''}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
