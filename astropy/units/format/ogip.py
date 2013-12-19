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

from .generic import Generic
from . import utils
from ...utils.compat.fractions import Fraction


class OGIP(Generic):
    """
    Support the units in `Office of Guest Investigator Programs (OGIP)
    FITS files
    <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__.
    """
    def __init__(self):
        # Build this on the class, so it only gets generated once.
        if not '_units' in OGIP.__dict__:
            OGIP._units, OGIP._deprecated_units, OGIP._functions = self._generate_unit_names()

        if '_parser' not in OGIP.__dict__:
            OGIP._parser, OGIP._lexer = self._make_parser()

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
        mCrab = 10 ** -3 * Crab
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
    def _make_parser(cls):
        """
        The grammar here is based on the description in the
        `Specification of Physical Units within OGIP FITS files
        <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__,
        which is not terribly precise.  The exact grammar is here is
        based on the YACC grammar in the `unity library
        <https://bitbucket.org/nxg/unity/>`_.
        """
        from ...extern.ply import lex, yacc

        tokens = (
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

        try:
            from . import ogip_lextab
            lexer = lex.lex(optimize=True, lextab=ogip_lextab)
        except ImportError:
            lexer = lex.lex(optimize=True, lextab='ogip_lextab',
                            outputdir=os.path.dirname(__file__))

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

        try:
            from . import ogip_parsetab
            parser = yacc.yacc(debug=False, tabmodule=ogip_parsetab,
                               write_tables=False)
        except ImportError:
            parser = yacc.yacc(debug=False, tabmodule='ogip_parsetab',
                               outputdir=os.path.dirname(__file__))

        return parser, lexer

    @classmethod
    def _get_unit(cls, t):
        try:
            return cls._parse_unit(t.value)
        except ValueError:
            raise ValueError(
                "At col {0}, {1!r} is not a valid unit according to the "
                "OGIP standard".format(
                    t.lexpos, t.value))

    @classmethod
    def _parse_unit(cls, unit):
        if unit not in cls._units:
            raise ValueError(
                "Unit {0!r} not supported by the OGIP "
                "standard".format(unit))

        if unit in cls._deprecated_units:
            from ..core import UnitsWarning
            warnings.warn(
                "The unit {0!r} is discouraged in the OGIP "
                "standard".format(unit),
                UnitsWarning)

        return cls._units[unit]

    def parse(self, s, debug=False):
        s = s.strip()
        try:
            # This is a short circuit for the case where the string is
            # just a single unit name
            return self._parse_unit(s)
        except ValueError:
            from ..core import Unit
            try:
                return Unit(
                    self._parser.parse(s, lexer=self._lexer, debug=debug))
            except ValueError as e:
                if six.text_type(e):
                    raise ValueError("{0} in unit {1!r}".format(
                        six.text_type(e), s))
                else:
                    raise ValueError(
                        "Syntax error parsing unit {0!r}".format(s))

    def _get_unit_name(self, unit):
        name = unit.get_format_name('ogip')

        if name not in self._units:
            raise ValueError(
                "Unit {0!r} is not part of the OGIP standard".format(name))

        if unit in self._deprecated_units:
            from ..core import UnitsWarning
            warnings.warn(
                "The unit {0!r} is discouraged in the OGIP "
                "standard".format(unit),
                UnitsWarning)

        return name

    def _format_unit_list(self, units):
        out = []
        units.sort(key=lambda x: self._get_unit_name(x[0]).lower())

        for base, power in units:
            if power == 1:
                out.append(self._get_unit_name(base))
            else:
                power = utils.format_power(power)
                if '/' in power:
                    out.append('{0}**({1})'.format(
                        self._get_unit_name(base), power))
                else:
                    out.append('{0}**{1}'.format(
                        self._get_unit_name(base), power))
        return ' '.join(out)

    def to_string(self, unit):
        from .. import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, self._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            # Can't use np.log10 here, because p[0] may be a Python long.
            if math.log10(unit.scale) % 1.0 != 0.0:
                warnings.warn(
                    "'{0}' scale should be a power of 10 in "
                    "OGIP format".format(
                        unit.scale),
                    core.UnitsWarning)

        return Generic.to_string(self, unit)
