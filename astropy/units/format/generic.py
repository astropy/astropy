# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles a "generic" string format for units
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six

import os
import re
import warnings

from . import utils
from .base import Base
from ...utils.misc import did_you_mean


def _to_string(cls, unit):
    from .. import core

    if isinstance(unit, core.CompositeUnit):
        parts = []

        if cls._show_scale and unit.scale != 1:
            parts.append('{0:g}'.format(unit.scale))

        if len(unit.bases):
            positives, negatives = utils.get_grouped_by_powers(
                unit.bases, unit.powers)
            if len(positives):
                parts.append(cls._format_unit_list(positives))
            elif len(parts) == 0:
                parts.append('1')

            if len(negatives):
                parts.append('/')
                unit_list = cls._format_unit_list(negatives)
                if len(negatives) == 1:
                    parts.append('{0}'.format(unit_list))
                else:
                    parts.append('({0})'.format(unit_list))

        return ' '.join(parts)
    elif isinstance(unit, core.NamedUnit):
        return cls._get_unit_name(unit)


class Generic(Base):
    """
    A "generic" format.

    The syntax of the format is based directly on the FITS standard,
    but instead of only supporting the units that FITS knows about, it
    supports any unit available in the `astropy.units` namespace.
    """

    _show_scale = True

    def __init__(self):
        # Build this on the class, so it only gets generated once.
        if '_parser' not in Generic.__dict__:
            Generic._parser, Generic._lexer = self._make_parser()

    @classmethod
    def _make_parser(cls):
        """
        The grammar here is based on the description in the `FITS
        standard
        <http://fits.gsfc.nasa.gov/standard30/fits_standard30aa.pdf>`_,
        Section 4.3, which is not terribly precise.  The exact grammar
        is here is based on the YACC grammar in the `unity library
        <https://bitbucket.org/nxg/unity/>`_.

        This same grammar is used by the `"fits"` and `"vounit"`
        formats, the only difference being the set of available unit
        strings.
        """
        from ...extern.ply import lex, yacc

        tokens = (
            'DOUBLE_STAR',
            'STAR',
            'PERIOD',
            'SOLIDUS',
            'CARET',
            'OPEN_PAREN',
            'CLOSE_PAREN',
            'FUNCNAME',
            'UNIT',
            'SIGN',
            'UINT',
            'UFLOAT'
            )

        t_STAR = r'\*'
        t_PERIOD = r'\.'
        t_SOLIDUS = r'/'
        t_DOUBLE_STAR = r'\*\*'
        t_CARET = r'\^'
        t_OPEN_PAREN = r'\('
        t_CLOSE_PAREN = r'\)'

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens
        def t_UFLOAT(t):
            r'((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+)?'
            if not re.search(r'[eE\.]', t.value):
                t.type = 'UINT'
                t.value = int(t.value)
            elif t.value.endswith('.'):
                t.type = 'UINT'
                t.value = int(t.value[:-1])
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

        # This needs to be a function so we can force it to happen
        # before t_UNIT
        def t_FUNCNAME(t):
            r'(sqrt)|(ln)|(exp)|(log)'
            return t

        def t_UNIT(t):
            r"%|([YZEPTGMkhdcmunpfazy]?'((?!\d)\w)+')|((?!\d)\w)+"
            t.value = cls._get_unit(t)
            return t

        t_ignore = ' '

        # Error handling rule
        def t_error(t):
            raise ValueError(
                "Invalid character at col {0}".format(t.lexpos))

        try:
            from . import generic_lextab
            lexer = lex.lex(optimize=True, lextab=generic_lextab,
                            reflags=re.UNICODE)
        except ImportError:
            lexer = lex.lex(optimize=True, lextab='generic_lextab',
                            outputdir=os.path.dirname(__file__),
                            reflags=re.UNICODE)

        def p_main(p):
            '''
            main : product_of_units
                 | factor product_of_units
                 | division_product_of_units
                 | factor division_product_of_units
                 | inverse_unit
                 | factor inverse_unit
                 | factor
            '''
            from ..core import Unit
            if len(p) == 2:
                p[0] = Unit(p[1])
            else:
                p[0] = Unit(p[1] * p[2])

        def p_division_product_of_units(p):
            '''
            division_product_of_units : division_product_of_units division product_of_units
                                      | product_of_units
            '''
            from ..core import Unit
            if len(p) == 4:
                p[0] = Unit(p[1] / p[3])
            else:
                p[0] = p[1]

        def p_inverse_unit(p):
            '''
            inverse_unit : division unit_expression
            '''
            p[0] = p[2] ** -1

        def p_factor(p):
            '''
            factor : factor_float
                   | factor_int
            '''
            p[0] = p[1]

        def p_factor_float(p):
            '''
            factor_float : signed_float
                         | signed_float UINT signed_int
                         | signed_float UINT power numeric_power
            '''
            if len(p) == 4:
                p[0] = p[1] * p[2] ** float(p[3])
            elif len(p) == 5:
                p[0] = p[1] * p[2] ** float(p[4])
            elif len(p) == 2:
                p[0] = p[1]

        def p_factor_int(p):
            '''
            factor_int : UINT
                       | UINT signed_int
                       | UINT power numeric_power
                       | UINT UINT signed_int
                       | UINT UINT power numeric_power
            '''
            if len(p) == 2:
                p[0] = p[1]
            elif len(p) == 3:
                p[0] = p[1] ** float(p[2])
            elif len(p) == 4:
                if isinstance(p[2], int):
                    p[0] = p[1] * p[2] ** float(p[3])
                else:
                    p[0] = p[1] ** float(p[3])
            elif len(p) == 5:
                p[0] = p[1] * p[2] ** p[4]

        def p_product_of_units(p):
            '''
            product_of_units : unit_expression product product_of_units
                             | unit_expression product_of_units
                             | unit_expression
            '''
            if len(p) == 2:
                p[0] = p[1]
            elif len(p) == 3:
                p[0] = p[1] * p[2]
            else:
                p[0] = p[1] * p[3]

        def p_unit_expression(p):
            '''
            unit_expression : function
                            | unit_with_power
                            | OPEN_PAREN product_of_units CLOSE_PAREN
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = p[2]

        def p_unit_with_power(p):
            '''
            unit_with_power : UNIT power numeric_power
                            | UNIT numeric_power
                            | UNIT
            '''
            if len(p) == 2:
                p[0] = p[1]
            elif len(p) == 3:
                p[0] = p[1] ** p[2]
            else:
                p[0] = p[1] ** p[3]

        def p_numeric_power(p):
            '''
            numeric_power : sign UINT
                          | OPEN_PAREN paren_expr CLOSE_PAREN
            '''
            if len(p) == 3:
                p[0] = p[1] * p[2]
            elif len(p) == 4:
                p[0] = p[2]

        def p_paren_expr(p):
            '''
            paren_expr : sign UINT
                       | signed_float
                       | frac
            '''
            if len(p) == 3:
                p[0] = p[1] * p[2]
            else:
                p[0] = p[1]

        def p_frac(p):
            '''
            frac : sign UINT division sign UINT
            '''
            p[0] = (p[1] * p[2]) / (p[4] * p[5])

        def p_sign(p):
            '''
            sign : SIGN
                 |
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_product(p):
            '''
            product : STAR
                    | PERIOD
            '''
            pass

        def p_division(p):
            '''
            division : SOLIDUS
            '''
            pass

        def p_power(p):
            '''
            power : DOUBLE_STAR
                  | CARET
            '''
            pass

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

        def p_function_name(p):
            '''
            function_name : FUNCNAME
            '''
            p[0] = p[1]

        def p_function(p):
            '''
            function : function_name OPEN_PAREN unit_expression CLOSE_PAREN
            '''
            if p[1] == 'sqrt':
                p[0] = p[3] ** 0.5
            else:
                raise ValueError(
                   "'{0}' is not a recognized function".format(p[1]))

        def p_error(p):
            raise ValueError()

        try:
            from . import generic_parsetab
            parser = yacc.yacc(debug=False, tabmodule=generic_parsetab,
                               write_tables=False)
        except ImportError:
            parser = yacc.yacc(debug=False, tabmodule='generic_parsetab',
                               outputdir=os.path.dirname(__file__))

        return parser, lexer

    @classmethod
    def _get_unit(cls, t):
        try:
            return cls._parse_unit(t.value)
        except ValueError as e:
            raise ValueError(
                "At col {0}, {1}".format(
                    t.lexpos, six.text_type(e)))

    @classmethod
    def _parse_unit(cls, s, detailed_exception=True):
        from ..core import get_current_unit_registry
        registry = get_current_unit_registry().registry
        if s == '%':
            return registry['percent']
        elif s in registry:
            return registry[s]

        if detailed_exception:
            raise ValueError(
                '{0} is not a valid unit. {1}'.format(
                    s, did_you_mean(s, registry)))
        else:
            raise ValueError()

    def parse(self, s, debug=False):
        if not isinstance(s, six.text_type):
            s = s.decode('ascii')

        result = self._do_parse(s, debug=debug)
        if s.count('/') > 1:
            from ..core import UnitsWarning
            warnings.warn(
                "'{0}' contains multiple slashes, which is "
                "discouraged by the FITS standard".format(s),
                UnitsWarning)
        return result

    def _do_parse(self, s, debug=False):
        try:
            # This is a short circuit for the case where the string
            # is just a single unit name
            return self._parse_unit(s, detailed_exception=False)
        except ValueError as e:
            try:
                return self._parser.parse(s, lexer=self._lexer, debug=debug)
            except ValueError as e:
                if six.text_type(e):
                    raise
                else:
                    raise ValueError(
                        "Syntax error parsing unit '{0}'".format(s))

    @classmethod
    def _get_unit_name(cls, unit):
        return unit.get_format_name('generic')

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
                    out.append('{0}({1})'.format(
                        cls._get_unit_name(base), power))
                else:
                    out.append('{0}{1}'.format(
                        cls._get_unit_name(base), power))
        return ' '.join(out)

    @classmethod
    def to_string(cls, unit):
        return _to_string(cls, unit)


class Unscaled(Generic):
    """
    A format that doesn't display the scale part of the unit, other
    than that, it is identical to the `Generic` format.

    This is used in some error messages where the scale is irrelevant.
    """
    _show_scale = False
