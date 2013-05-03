# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICNSE.rst

"""
Handles a "generic" string format for units
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import re

from .base import Base
from . import utils


# TODO: Support logarithmic units using bracketed syntax

class CDS(Base):
    """
    Support the `Centre de Données astronomiques de Strasbourg
    <http://cds.u-strasbg.fr/>`_ `Standards for Astronomical
    Catalogues 2.0 <http://cds.u-strasbg.fr/doc/catstd-3.2.htx>`_
    format.  This format is used by VOTable up to version 1.2.
    """
    def __init__(self):
        # Build this on the class, so it only gets generated once.
        if not '_units' in CDS.__dict__:
            CDS._units = self._generate_unit_names()

        if '_parser' not in CDS.__dict__:
            CDS._parser, CDS._lexer = self._make_parser()

    @staticmethod
    def _generate_unit_names():
        import keyword
        from ... import units as u
        names = {}

        names['%'] = u.Unit(0.01)

        bases = [
            'A', 'C', 'cd', 'eV', 'F', 'g', 'H', 'Hz', 'J', 'K',
            'lm', 'lx', 'm', 'mol', 'N', 'Ohm', 'Pa', 'rad', 's', 'S',
            'sr', 'T', 'V', 'W', 'Wb']

        prefixes = [
            'y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', 'c', 'd',
            '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']

        for base in bases:
            for prefix in prefixes:
                key = prefix + base
                if keyword.iskeyword(key):
                    continue
                names[key] = getattr(u, key)

        simple_bases = [
            'a', 'AU', 'arcmin', 'arcsec', 'barn', 'bit',
            'byte', 'ct', 'D', 'd', 'deg', 'h', 'Jy', 'mag', 'mas',
            'min', 'pc', 'pix', 'Ry', 'solLum', 'solMass', 'solRad',
            'Sun', 'yr']

        for base in simple_bases:
            names[base] = getattr(u, base)

        return names

    @classmethod
    def _make_parser(cls):
        """
        The grammar here is based on the description in the `Standards
        for Astronomical Catalogues 2.0
        <http://cds.u-strasbg.fr/doc/catstd-3.2.htx>`_, which is not
        terribly precise.  The exact grammar is here is based on the
        YACC grammar in the `unity library
        <https://bitbucket.org/nxg/unity/>`_.
        """
        from ...extern.ply import lex, yacc

        tokens = (
            'PRODUCT',
            'DIVISION',
            'OPEN_PAREN',
            'CLOSE_PAREN',
            'X',
            'SIGN',
            'UINT',
            'UFLOAT',
            'UNIT'
            )

        t_PRODUCT = r'\.'
        t_DIVISION = r'/'
        t_OPEN_PAREN = r'\('
        t_CLOSE_PAREN = r'\)'

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens
        def t_UFLOAT(t):
            r'((\d+\.?\d+)|(\.\d+))([eE][+-]?\d+)?'
            if not re.search(r'[eE\.]', t.value):
                t.type = 'UINT'
                t.value = int(t.value)
            else:
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

        def t_X(t):
            r'x'
            return t

        def t_UNIT(t):
            r'\%|([a-zA-Z][a-zA-Z_]*)'
            t.value = cls._get_unit(t)
            return t

        t_ignore = ''

        # Error handling rule
        def t_error(t):
            raise ValueError(
                "Invalid character at col {0}".format(t.lexpos))

        lexer = lex.lex()

        def p_main(p):
            '''
            main : factor combined_units
                 | combined_units
            '''
            if len(p) == 3:
                from ..core import Unit
                p[0] = Unit(p[1] * p[2])
            else:
                p[0] = p[1]

        def p_combined_units(p):
            '''
            combined_units : product_of_units
                           | division_of_units
            '''
            p[0] = p[1]

        def p_product_of_units(p):
            '''
            product_of_units : unit_expression PRODUCT combined_units
                             | unit_expression
            '''
            if len(p) == 4:
                p[0] = p[1] * p[3]
            else:
                p[0] = p[1]

        def p_division_of_units(p):
            '''
            division_of_units : DIVISION unit_expression
                              | unit_expression DIVISION combined_units
            '''
            if len(p) == 3:
                p[0] = p[1] ** -1
            else:
                p[0] = p[1] / p[3]

        def p_unit_expression(p):
            '''
            unit_expression : unit_with_power
                            | OPEN_PAREN combined_units CLOSE_PAREN
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = p[2]

        def p_factor(p):
            '''
            factor : signed_float X UINT signed_int
                   | UINT X UINT signed_int
                   | UINT signed_int
                   | UINT
                   | signed_float
            '''
            if len(p) == 5:
                if p[3] != 10:
                    raise ValueError(
                        "Only base ten exponents are allowed in CDS")
                p[0] = p[1] * 10.0 ** p[4]
            elif len(p) == 3:
                if p[1] != 10:
                    raise ValueError(
                        "Only base ten exponents are allowed in CDS")
                p[0] = 10.0 ** p[2]
            elif len(p) == 2:
                p[0] = p[1]

        def p_unit_with_power(p):
            '''
            unit_with_power : UNIT numeric_power
                            | UNIT
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = p[1] ** p[2]

        def p_numeric_power(p):
            '''
            numeric_power : sign UINT
            '''
            p[0] = p[1] * p[2]

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

        parser = yacc.yacc()

        return parser, lexer

    @classmethod
    def _get_unit(cls, t):
        try:
            return cls._parse_unit(t.value)
        except ValueError:
            raise ValueError(
                "At col {0}, {1!r} is not a valid unit".format(
                    t.lexpos, t.value))

    @classmethod
    def _parse_unit(cls, unit):
        if unit not in cls._units:
            raise ValueError(
                "Unit {0!r} not supported by the CDS SAC "
                "standard.".format(unit))

        return cls._units[unit]

    def parse(self, s, debug=False):
        if ' ' in s:
            raise ValueError('CDS unit must not contain whitespace')

        # This is a short circuit for the case where the string
        # is just a single unit name
        try:
            return self._parse_unit(s)
        except ValueError:
            try:
                return self._parser.parse(s, lexer=self._lexer, debug=debug)
            except ValueError as e:
                if str(e):
                    raise ValueError("{0} in unit {1!r}".format(
                        str(e), s))
                else:
                    raise ValueError(
                        "Syntax error parsing unit {0!r}".format(s))

    def _get_unit_name(self, unit):
        return unit.get_format_name('cds')

    def _format_unit_list(self, units):
        out = []
        for base, power in units:
            if power == 1:
                out.append(self._get_unit_name(base))
            else:
                out.append('{0}{1}'.format(
                    self._get_unit_name(base), int(power)))
        return '.'.join(out)

    def to_string(self, unit):
        from .. import core

        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, self._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            s = ''
            if unit.scale != 1:
                m, e = utils.split_mantissa_exponent(unit.scale)
                parts = []
                if m:
                    parts.append(m)
                if e:
                    if not e.startswith('-'):
                        e = "+" + e
                    parts.append('10{0}'.format(e))
                s = 'x'.join(parts)
            else:
                s = ''

            pairs = zip(unit.bases, unit.powers)
            pairs.sort(key=lambda x: x[1], reverse=True)

            s += self._format_unit_list(pairs)
        elif isinstance(unit, core.NamedUnit):
            s = self._get_unit_name(unit)

        return s
