# generic_lextab.py. This file automatically created by PLY (version 3.4). Don't edit!
_tabversion   = '3.4'
_lextokens    = {'STAR': 1, 'UINT': 1, 'DOUBLE_STAR': 1, 'CARET': 1, 'UNIT': 1, 'OPEN_PAREN': 1, 'UFLOAT': 1, 'SOLIDUS': 1, 'PERIOD': 1, 'SIGN': 1, 'CLOSE_PAREN': 1, 'FUNCNAME': 1}
_lexreflags   = 32
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [("(?P<t_UFLOAT>((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_FUNCNAME>((sqrt)|(ln)|(exp)|(log)|(mag)|(dB)|(dex))(?=\\ *\\())|(?P<t_UNIT>%|([YZEPTGMkhdcmunpfazy]?'((?!\\d)\\w)+')|((?!\\d)\\w)+)|(?P<t_DOUBLE_STAR>\\*\\*)|(?P<t_OPEN_PAREN>\\()|(?P<t_CLOSE_PAREN>\\))|(?P<t_PERIOD>\\.)|(?P<t_STAR>\\*)|(?P<t_CARET>\\^)|(?P<t_SOLIDUS>/)", [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_FUNCNAME', 'FUNCNAME'), None, None, None, None, None, None, None, None, ('t_UNIT', 'UNIT'), None, None, None, (None, 'DOUBLE_STAR'), (None, 'OPEN_PAREN'), (None, 'CLOSE_PAREN'), (None, 'PERIOD'), (None, 'STAR'), (None, 'CARET'), (None, 'SOLIDUS')])]}
_lexstateignore = {'INITIAL': ' '}
_lexstateerrorf = {'INITIAL': 't_error'}
