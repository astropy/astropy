# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICNSE.rst

"""
Handles units in `Office of Guest Investigator Programs (OGIP)
FITS files
<http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six

import keyword
import math
import os
import warnings

from . import core, generic, utils
from ...utils.compat.fractions import Fraction


class OGIP(generic.Generic):
    """
    Support the units in `Office of Guest Investigator Programs (OGIP)
    FITS files
    <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
    """

    _tokens = (
        'DIVISION',
        'OPEN_PAREN',
        'CLOSE_PAREN',
        'WHITESPACE',
        'STARSTAR',
        'STAR',
        'SIGN',
        'UFLOAT',
        'LIT10',
        'UINT',
        'UNKNOWN',
        'UNIT'
    )

    @staticmethod
    def _generate_unit_names():

        from ... import units as u
        names = {}
        deprecated_names = set()

        bases = [
            'A', 'C', 'cd', 'eV', 'F', 'g', 'H', 'Hz', 'J',
            'Jy', 'K', 'lm', 'lx', 'm', 'mol', 'N', 'ohm', 'Pa',
            'pc', 'rad', 's', 'S', 'sr', 'T', 'V', 'W', 'Wb'
        ]
        deprecated_bases = []
        prefixes = [
            'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', 'c', 'd',
            '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'
        ]

        for base in bases + deprecated_bases:
            for prefix in prefixes:
                key = prefix + base
                if keyword.iskeyword(key):
                    continue
                names[key] = getattr(u, key)
        for base in deprecated_bases:
            for prefix in prefixes:
                deprecated_names.add(prefix + base)

        simple_units = [
            'angstrom', 'arcmin', 'arcsec', 'AU', 'barn', 'bin',
            'byte', 'chan', 'count', 'day', 'deg', 'erg', 'G',
            'h', 'lyr', 'mag', 'min', 'photon', 'pixel',
            'voxel', 'yr'
        ]
        for unit in simple_units:
            names[unit] = getattr(u, unit)

        # Create a separate, disconnected unit for the special case of
        # Crab and mCrab, since OGIP doesn't define their quantities.
        Crab = u.def_unit(['Crab'], prefixes=False, doc='Crab (X-ray flux)')
        mCrab = u.Unit(10 ** -3 * Crab)
        names['Crab'] = Crab
        names['mCrab'] = mCrab

        deprecated_units = ['Crab', 'mCrab']
        for unit in deprecated_units:
            deprecated_names.add(unit)

        # Define the function names, so we can parse them, even though
        # we can't use any of them (other than sqrt) meaningfully for
        # now.
        functions = [
            'log', 'ln', 'exp', 'sqrt', 'sin', 'cos', 'tan', 'asin',
            'acos', 'atan', 'sinh', 'cosh', 'tanh'
        ]
        for name in functions:
            names[name] = name

        return names, deprecated_names, functions

    @classmethod
    def _make_lexer(cls):
        from ...extern.ply import lex

        tokens = cls._tokens

        t_DIVISION = r'/'
        t_OPEN_PAREN = r'\('
        t_CLOSE_PAREN = r'\)'
        t_WHITESPACE = '[ \t]+'
        t_STARSTAR = r'\*\*'
        t_STAR = r'\*'

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens
        def t_UFLOAT(t):
            r'(((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+))|(((\d+\.\d*)|(\.\d+))([eE][+-]?\d+)?)'
            t.value = float(t.value)
            return t

        def t_UINT(t):
            r'\d+'
            t.value = int(t.value)
            return t

        def t_SIGN(t):
            r'[+-](?=\d)'
            t.value = float(t.value + '1')
            return t

        def t_X(t):  # multiplication for factor in front of unit
            r'[x×]'
            return t

        def t_LIT10(t):
            r'10'
            return 10

        def t_UNKNOWN(t):
            r'[Uu][Nn][Kk][Nn][Oo][Ww][Nn]'
            return None

        def t_UNIT(t):
            r'[a-zA-Z][a-zA-Z_]*'
            t.value = cls._get_unit(t)
            return t

        # Don't ignore whitespace
        t_ignore = ''

        # Error handling rule
        def t_error(t):
            raise ValueError(
                "Invalid character at col {0}".format(t.lexpos))

        lexer = lex.lex(optimize=True, lextab='ogip_lextab',
                        outputdir=os.path.dirname(__file__))

        return lexer

    @classmethod
    def _make_parser(cls):
        """
        The grammar here is based on the description in the
        `Specification of Physical Units within OGIP FITS files
        <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__,
        which is not terribly precise.  The exact grammar is here is
        based on the YACC grammar in the `unity library
        <https://bitbucket.org/nxg/unity/>`_.
        """

        from ...extern.ply import yacc

        tokens = cls._tokens

        def p_main(p):
            '''
            main : UNKNOWN
                 | complete_expression
                 | scale_factor complete_expression
                 | scale_factor WHITESPACE complete_expression
            '''
            if len(p) == 4:
                p[0] = p[1] * p[3]
            elif len(p) == 3:
                p[0] = p[1] * p[2]
            else:
                p[0] = p[1]

        def p_complete_expression(p):
            '''
            complete_expression : product_of_units
            '''
            p[0] = p[1]

        def p_product_of_units(p):
            '''
            product_of_units : unit_expression
                             | division unit_expression
                             | product_of_units product unit_expression
                             | product_of_units division unit_expression
            '''
            if len(p) == 4:
                if p[2] == 'DIVISION':
                    p[0] = p[1] / p[3]
                else:
                    p[0] = p[1] * p[3]
            elif len(p) == 3:
                p[0] = p[2] ** -1
            else:
                p[0] = p[1]

        def p_unit_expression(p):
            '''
            unit_expression : unit
                            | UNIT OPEN_PAREN complete_expression CLOSE_PAREN
                            | OPEN_PAREN complete_expression CLOSE_PAREN
                            | UNIT OPEN_PAREN complete_expression CLOSE_PAREN power numeric_power
                            | OPEN_PAREN complete_expression CLOSE_PAREN power numeric_power
            '''
            if p[1] in cls._functions and p[1] != 'sqrt':
                raise ValueError(
                    "The function '{0}' is valid in OGIP, but not understood "
                    "by astropy.units.".format(
                        p[1]))

            if len(p) == 7:
                if p[1] == 'sqrt':
                    p[0] = p[1] * p[3] ** (0.5 * p[6])
                else:
                    p[0] = p[1] * p[3] ** p[6]
            elif len(p) == 6:
                p[0] = p[2] ** p[5]
            elif len(p) == 5:
                if p[1] == 'sqrt':
                    p[0] = p[3] ** 0.5
                else:
                    p[0] = p[1] * p[3]
            elif len(p) == 4:
                p[0] = p[2]
            else:
                p[0] = p[1]

        def p_scale_factor(p):
            '''
            scale_factor : LIT10 power numeric_power
                         | LIT10
                         | signed_float
                         | signed_float power numeric_power
                         | signed_int power numeric_power
            '''
            if len(p) == 4:
                p[0] = 10 ** p[3]
            else:
                p[0] = p[1]
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(p[0]) % 1.0 != 0.0:
                from ..core import UnitsWarning
                warnings.warn(
                    "'{0}' scale should be a power of 10 in "
                    "OGIP format".format(p[0]), UnitsWarning)

        def p_division(p):
            '''
            division : DIVISION
                     | WHITESPACE DIVISION
                     | WHITESPACE DIVISION WHITESPACE
                     | DIVISION WHITESPACE
            '''
            p[0] = 'DIVISION'

        def p_product(p):
            '''
            product : WHITESPACE
                    | STAR
                    | WHITESPACE STAR
                    | WHITESPACE STAR WHITESPACE
                    | STAR WHITESPACE
            '''
            p[0] = 'PRODUCT'

        def p_power(p):
            '''
            power : STARSTAR
            '''
            p[0] = 'POWER'

        def p_unit(p):
            '''
            unit : UNIT
                 | UNIT power numeric_power
            '''
            if len(p) == 4:
                p[0] = p[1] ** p[3]
            else:
                p[0] = p[1]

        def p_numeric_power(p):
            '''
            numeric_power : UINT
                          | signed_float
                          | OPEN_PAREN signed_int CLOSE_PAREN
                          | OPEN_PAREN signed_float CLOSE_PAREN
                          | OPEN_PAREN signed_float division UINT CLOSE_PAREN
            '''
            if len(p) == 6:
                p[0] = Fraction(int(p[2]), int(p[4]))
            elif len(p) == 4:
                p[0] = p[2]
            else:
                p[0] = p[1]

        def p_sign(p):
            '''
            sign : SIGN
                 |
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_signed_int(p):
            '''
            signed_int : SIGN UINT
            '''
            p[0] = p[1] * p[2]

        def p_signed_float(p):
            '''
            signed_float : sign UINT
                         | sign UFLOAT
            '''
            p[0] = p[1] * p[2]

        def p_error(p):
            raise ValueError()

        parser = yacc.yacc(debug=False, tabmodule='ogip_parsetab',
                           outputdir=os.path.dirname(__file__),
                           write_tables=True)

        return parser

    @classmethod
    def _get_unit(cls, t):
        try:
            return cls._parse_unit(t.value)
        except ValueError as e:
            raise ValueError(
                "At col {0}, '{1}': {2}".format(
                    t.lexpos, t.value, six.text_type(e)))

    @classmethod
    def _validate_unit(cls, unit, detailed_exception=True):
        if unit not in cls._units:
            if detailed_exception:
                raise ValueError(
                    "Unit '{0}' not supported by the OGIP "
                    "standard. {1}".format(
                        unit, utils.did_you_mean_units(
                            unit, cls._units, cls._deprecated_units,
                            cls._to_decomposed_alternative)))
            else:
                raise ValueError()

        if unit in cls._deprecated_units:
            utils.unit_deprecation_warning(
                unit, cls._units[unit], 'OGIP',
                cls._to_decomposed_alternative)

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        cls._validate_unit(unit, detailed_exception=detailed_exception)
        return cls._units[unit]

    @classmethod
    def parse(cls, s, debug=False):
        s = s.strip()
        try:
            # This is a short circuit for the case where the string is
            # just a single unit name
            return cls._parse_unit(s, detailed_exception=False)
        except ValueError:
            try:
                return core.Unit(
                    cls._parser.parse(s, lexer=cls._lexer, debug=debug))
            except ValueError as e:
                if six.text_type(e):
                    raise
                else:
                    raise ValueError(
                        "Syntax error parsing unit '{0}'".format(s))

    @classmethod
    def _get_unit_name(cls, unit):
        name = unit.get_format_name('ogip')
        cls._validate_unit(name)
        return name

    @classmethod
    def _format_unit_list(cls, units):
        out = []
        units.sort(key=lambda x: cls._get_unit_name(x[0]).lower())

        for base, power in units:
            if power == 1:
                out.append(cls._get_unit_name(base))
            else:
                power = utils.format_power(power)
                if '/' in power:
                    out.append('{0}**({1})'.format(
                        cls._get_unit_name(base), power))
                else:
                    out.append('{0}**{1}'.format(
                        cls._get_unit_name(base), power))
        return ' '.join(out)

    @classmethod
    def to_string(cls, unit):
        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(unit.scale) % 1.0 != 0.0:
                warnings.warn(
                    "'{0}' scale should be a power of 10 in "
                    "OGIP format".format(
                        unit.scale),
                    core.UnitsWarning)

        return generic._to_string(cls, unit)

    @classmethod
    def _to_decomposed_alternative(cls, unit):
        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(unit.scale) % 1.0 != 0.0:
                scale = unit.scale
                unit = copy.copy(unit)
                unit._scale = 1.0
                return '{0} (with data multiplied by {1})'.format(
                    generic._to_string(cls, unit), scale)

        return generic._to_string(unit)
