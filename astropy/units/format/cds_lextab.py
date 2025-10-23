# cds_lextab.py. This file automatically created by PLY (version 3.11). Don't edit!
_lexmap = {
    't_UFLOAT': 0,
    't_UINT': 1,
    't_SIGN': 2,
    't_X': 3,
    't_UNIT': 4,
    't_DIMENSIONLESS': 5,
}

states = {}

lexer = {
    'literals': '',
    'reflags': 32,
    'rules': [
        (("((\\d+\\.?\\d+)|(\\.\\d+))([eE][+-]?\\d+)?"), 't_UFLOAT'),
        (("\\d+"), 't_UINT'),
        (("[+-](?=\\d)"), 't_SIGN'),
        (("[x×]"), 't_X'),
        (("\\%|°|\\\\h|((?!\\d)\\w)+"), 't_UNIT'),
        (("---|-"), 't_DIMENSIONLESS'),
    ],
    'error': 't_error'
}

