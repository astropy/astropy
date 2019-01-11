# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/units
#
# You can then commit the changes to this file.

# generic_lextab.py. This file automatically created by PLY (version 3.10). Don't edit!
_tabversion   = '3.10'
_lextokens    = set(('CARET', 'CLOSE_PAREN', 'DOUBLE_STAR', 'FUNCNAME', 'OPEN_PAREN', 'PERIOD', 'SIGN', 'SOLIDUS', 'STAR', 'UFLOAT', 'UINT', 'UNIT'))
_lexreflags   = 32
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [("(?P<t_UFLOAT>((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_FUNCNAME>((sqrt)|(ln)|(exp)|(log)|(mag)|(dB)|(dex))(?=\\ *\\())|(?P<t_UNIT>%|([YZEPTGMkhdcmunpfazy]?'((?!\\d)\\w)+')|((?!\\d)\\w)+)|(?P<t_DOUBLE_STAR>\\*\\*)|(?P<t_CLOSE_PAREN>\\))|(?P<t_OPEN_PAREN>\\()|(?P<t_CARET>\\^)|(?P<t_PERIOD>\\.)|(?P<t_STAR>\\*)|(?P<t_SOLIDUS>/)", [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_FUNCNAME', 'FUNCNAME'), None, None, None, None, None, None, None, None, ('t_UNIT', 'UNIT'), None, None, None, (None, 'DOUBLE_STAR'), (None, 'CLOSE_PAREN'), (None, 'OPEN_PAREN'), (None, 'CARET'), (None, 'PERIOD'), (None, 'STAR'), (None, 'SOLIDUS')])]}
_lexstateignore = {'INITIAL': ' '}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
