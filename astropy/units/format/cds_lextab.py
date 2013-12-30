# cds_lextab.py. This file automatically created by PLY (version 3.4). Don't edit!
_tabversion   = '3.4'
_lextokens    = {'DIVISION': 1, 'PRODUCT': 1, 'SIGN': 1, 'OPEN_PAREN': 1, 'UINT': 1, 'UNIT': 1, 'X': 1, 'CLOSE_PAREN': 1, 'UFLOAT': 1}
_lexreflags   = 32
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>((\\d+\\.?\\d+)|(\\.\\d+))([eE][+-]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+-](?=\\d))|(?P<t_X>[x\xd7])|(?P<t_UNIT>\\%|\xb0|\\\\h|((?!\\d)\\w)+)|(?P<t_CLOSE_PAREN>\\))|(?P<t_PRODUCT>\\.)|(?P<t_OPEN_PAREN>\\()|(?P<t_DIVISION>/)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_X', 'X'), ('t_UNIT', 'UNIT'), None, (None, 'CLOSE_PAREN'), (None, 'PRODUCT'), (None, 'OPEN_PAREN'), (None, 'DIVISION')])]}
_lexstateignore = {'INITIAL': ''}
_lexstateerrorf = {'INITIAL': 't_error'}
