# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICNSE.rst

# This module includes files automatically generated from ply (these end in
# _lextab.py and _parsetab.py). To generate these files, remove them from this
# folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/units
#
# You can then commit the changes to the re-generated _lextab.py and
# _parsetab.py files.

"""
Handles the CDS string format for units
"""


import operator
import os
import re


from .base import Base
from . import core, utils
from astropy.units.utils import is_effectively_unity
from astropy.utils import classproperty
from astropy.utils.misc import did_you_mean


# TODO: Support logarithmic units using bracketed syntax

class CDS(Base):
    """
    Support the `Centre de Données astronomiques de Strasbourg
    <http://cds.u-strasbg.fr/>`_ `Standards for Astronomical
    Catalogues 2.0 <http://vizier.u-strasbg.fr/vizier/doc/catstd-3.2.htx>`_
    format, and the `complete set of supported units
    <http://vizier.u-strasbg.fr/cgi-bin/Unit>`_.  This format is used
    by VOTable up to version 1.2.
    """

    _tokens = (
        'PRODUCT',
        'DIVISION',
        'OPEN_PAREN',
        'CLOSE_PAREN',
        'OPEN_BRACKET',
        'CLOSE_BRACKET',
        'X',
        'SIGN',
        'UINT',
        'UFLOAT',
        'UNIT'
    )

    @classproperty(lazy=True)
    def _units(cls):
        return cls._generate_unit_names()

    @classproperty(lazy=True)
    def _parser(cls):
        return cls._make_parser()

    @classproperty(lazy=True)
    def _lexer(cls):
        return cls._make_lexer()

    @staticmethod
    def _generate_unit_names():
        from astropy.units import cds
        from astropy import units as u

        names = {}

        for key, val in cds.__dict__.items():
            if isinstance(val, u.UnitBase):
                names[key] = val

        return names

    @classmethod
    def _make_lexer(cls):

        from astropy.extern.ply import lex

        tokens = cls._tokens

        t_PRODUCT = r'\.'
        t_DIVISION = r'/'
        t_OPEN_PAREN = r'\('
        t_CLOSE_PAREN = r'\)'
        t_OPEN_BRACKET = r'\['
        t_CLOSE_BRACKET = r'\]'

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

        def t_X(t):  # multiplication for factor in front of unit
            r'[x×]'
            return t

        def t_UNIT(t):
            r'\%|°|\\h|((?!\d)\w)+'
            t.value = cls._get_unit(t)
            return t

        t_ignore = ''

        # Error handling rule
        def t_error(t):
            raise ValueError(
                f"Invalid character at col {t.lexpos}")

        lexer_exists = os.path.exists(os.path.join(os.path.dirname(__file__),
                                      'cds_lextab.py'))

        lexer = lex.lex(optimize=True, lextab='cds_lextab',
                        outputdir=os.path.dirname(__file__),
                        reflags=int(re.UNICODE))

        if not lexer_exists:
            cls._add_tab_header('cds_lextab')

        return lexer

    @classmethod
    def _make_parser(cls):
        """
        The grammar here is based on the description in the `Standards
        for Astronomical Catalogues 2.0
        <http://vizier.u-strasbg.fr/vizier/doc/catstd-3.2.htx>`_, which is not
        terribly precise.  The exact grammar is here is based on the
        YACC grammar in the `unity library
        <https://bitbucket.org/nxg/unity/>`_.
        """

        from astropy.extern.ply import yacc

        tokens = cls._tokens

        def p_main(p):
            '''
            main : factor combined_units
                 | combined_units
                 | OPEN_BRACKET combined_units CLOSE_BRACKET
                 | factor
            '''
            from astropy.units.core import Unit
            from astropy.units import dex
            if len(p) == 3:
                p[0] = Unit(p[1] * p[2])
            elif len(p) == 4:
                p[0] = dex(p[2])
            else:
                p[0] = Unit(p[1])

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
                p[0] = p[2] ** -1
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

        parser_exists = os.path.exists(os.path.join(os.path.dirname(__file__),
                                       'cds_parsetab.py'))

        parser = yacc.yacc(debug=False, tabmodule='cds_parsetab',
                           outputdir=os.path.dirname(__file__),
                           optimize=True, write_tables=True)

        if not parser_exists:
            cls._add_tab_header('cds_parsetab')

        return parser

    @classmethod
    def _get_unit(cls, t):
        try:
            return cls._parse_unit(t.value)
        except ValueError as e:
            raise ValueError(
                "At col {}, {}".format(
                    t.lexpos, str(e)))

    @classmethod
    def _parse_unit(cls, unit, detailed_exception=True):
        if unit not in cls._units:
            if detailed_exception:
                raise ValueError(
                    "Unit '{}' not supported by the CDS SAC "
                    "standard. {}".format(
                        unit, did_you_mean(
                            unit, cls._units)))
            else:
                raise ValueError()

        return cls._units[unit]

    @classmethod
    def parse(cls, s, debug=False):
        if ' ' in s:
            raise ValueError('CDS unit must not contain whitespace')

        if not isinstance(s, str):
            s = s.decode('ascii')

        # This is a short circuit for the case where the string
        # is just a single unit name
        try:
            return cls._parse_unit(s, detailed_exception=False)
        except ValueError:
            try:
                return cls._parser.parse(s, lexer=cls._lexer, debug=debug)
            except ValueError as e:
                if str(e):
                    raise ValueError(str(e))
                else:
                    raise ValueError("Syntax error")

    @staticmethod
    def _get_unit_name(unit):
        return unit.get_format_name('cds')

    @classmethod
    def _format_unit_list(cls, units):
        out = []
        for base, power in units:
            if power == 1:
                out.append(cls._get_unit_name(base))
            else:
                out.append('{}{}'.format(
                    cls._get_unit_name(base), int(power)))
        return '.'.join(out)

    @classmethod
    def to_string(cls, unit):
        # Remove units that aren't known to the format
        unit = utils.decompose_to_known_units(unit, cls._get_unit_name)

        if isinstance(unit, core.CompositeUnit):
            if(unit.physical_type == 'dimensionless' and
               is_effectively_unity(unit.scale*100.)):
                return '%'

            if unit.scale == 1:
                s = ''
            else:
                m, e = utils.split_mantissa_exponent(unit.scale)
                parts = []
                if m not in ('', '1'):
                    parts.append(m)
                if e:
                    if not e.startswith('-'):
                        e = "+" + e
                    parts.append(f'10{e}')
                s = 'x'.join(parts)

            pairs = list(zip(unit.bases, unit.powers))
            if len(pairs) > 0:
                pairs.sort(key=operator.itemgetter(1), reverse=True)

                s += cls._format_unit_list(pairs)

        elif isinstance(unit, core.NamedUnit):
            s = cls._get_unit_name(unit)

        return s
