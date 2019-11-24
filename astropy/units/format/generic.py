# Licensed under a 3-clause BSD style license - see LICENSE.rst

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
Handles a "generic" string format for units
"""

import os
import re
import warnings
import sys
from fractions import Fraction
import unicodedata

from . import core, utils
from .base import Base
from astropy.utils import classproperty
from astropy.utils.misc import did_you_mean


def _is_ascii(s):
    if sys.version_info >= (3, 7, 0):
        return s.isascii()
    else:
        try:
            s.encode('ascii')
            return True
        except UnicodeEncodeError:
            return False


def _to_string(cls, unit):
    if isinstance(unit, core.CompositeUnit):
        parts = []

        if cls._show_scale and unit.scale != 1:
            parts.append(f'{unit.scale:g}')

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
                    parts.append(f'{unit_list}')
                else:
                    parts.append(f'({unit_list})')

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

    _tokens = (
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

    @classproperty(lazy=True)
    def _all_units(cls):
        return cls._generate_unit_names()

    @classproperty(lazy=True)
    def _units(cls):
        return cls._all_units[0]

    @classproperty(lazy=True)
    def _deprecated_units(cls):
        return cls._all_units[1]

    @classproperty(lazy=True)
    def _functions(cls):
        return cls._all_units[2]

    @classproperty(lazy=True)
    def _parser(cls):
        return cls._make_parser()

    @classproperty(lazy=True)
    def _lexer(cls):
        return cls._make_lexer()

    @classmethod
    def _make_lexer(cls):
        from astropy.extern.ply import lex

        tokens = cls._tokens

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
            t.value = int(t.value + '1')
            return t

        # This needs to be a function so we can force it to happen
        # before t_UNIT
        def t_FUNCNAME(t):
            r'((sqrt)|(ln)|(exp)|(log)|(mag)|(dB)|(dex))(?=\ *\()'
            return t

        def t_UNIT(t):
            "%|([YZEPTGMkhdcmu\N{MICRO SIGN}npfazy]?'((?!\\d)\\w)+')|((?!\\d)\\w)+"
            t.value = cls._get_unit(t)
            return t

        t_ignore = ' '

        # Error handling rule
        def t_error(t):
            raise ValueError(
                f"Invalid character at col {t.lexpos}")

        lexer_exists = os.path.exists(os.path.join(os.path.dirname(__file__),
                                      'generic_lextab.py'))

        lexer = lex.lex(optimize=True, lextab='generic_lextab',
                        outputdir=os.path.dirname(__file__),
                        reflags=int(re.UNICODE))

        if not lexer_exists:
            cls._add_tab_header('generic_lextab')

        return lexer

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
        from astropy.extern.ply import yacc

        tokens = cls._tokens

        def p_main(p):
            '''
            main : product_of_units
                 | factor product_of_units
                 | factor product product_of_units
                 | division_product_of_units
                 | factor division_product_of_units
                 | factor product division_product_of_units
                 | inverse_unit
                 | factor inverse_unit
                 | factor product inverse_unit
                 | factor
            '''
            from astropy.units.core import Unit
            if len(p) == 2:
                p[0] = Unit(p[1])
            elif len(p) == 3:
                p[0] = Unit(p[1] * p[2])
            elif len(p) == 4:
                p[0] = Unit(p[1] * p[3])

        def p_division_product_of_units(p):
            '''
            division_product_of_units : division_product_of_units division product_of_units
                                      | product_of_units
            '''
            from astropy.units.core import Unit
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
            factor : factor_fits
                   | factor_float
                   | factor_int
            '''
            p[0] = p[1]

        def p_factor_float(p):
            '''
            factor_float : signed_float
                         | signed_float UINT signed_int
                         | signed_float UINT power numeric_power
            '''
            if cls.name == 'fits':
                raise ValueError("Numeric factor not supported by FITS")
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
            if cls.name == 'fits':
                raise ValueError("Numeric factor not supported by FITS")
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

        def p_factor_fits(p):
            '''
            factor_fits : UINT power OPEN_PAREN signed_int CLOSE_PAREN
                        | UINT power OPEN_PAREN UINT CLOSE_PAREN
                        | UINT power signed_int
                        | UINT power UINT
                        | UINT SIGN UINT
                        | UINT OPEN_PAREN signed_int CLOSE_PAREN
            '''
            if p[1] != 10:
                if cls.name == 'fits':
                    raise ValueError("Base must be 10")
                else:
                    return
            if len(p) == 4:
                if p[2] in ('**', '^'):
                    p[0] = 10 ** p[3]
                else:
                    p[0] = 10 ** (p[2] * p[3])
            elif len(p) == 5:
                p[0] = 10 ** p[3]
            elif len(p) == 6:
                p[0] = 10 ** p[4]

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
            p[0] = Fraction(p[1] * p[2], p[4] * p[5])

        def p_sign(p):
            '''
            sign : SIGN
                 |
            '''
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1

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
            p[0] = p[1]

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
            function : function_name OPEN_PAREN main CLOSE_PAREN
            '''
            if p[1] == 'sqrt':
                p[0] = p[3] ** 0.5
                return
            elif p[1] in ('mag', 'dB', 'dex'):
                function_unit = cls._parse_unit(p[1])
                # In Generic, this is callable, but that does not have to
                # be the case in subclasses (e.g., in VOUnit it is not).
                if callable(function_unit):
                    p[0] = function_unit(p[3])
                    return

            raise ValueError("'{}' is not a recognized function".format(p[1]))

        def p_error(p):
            raise ValueError()

        parser_exists = os.path.exists(os.path.join(os.path.dirname(__file__),
                                       'generic_parsetab.py'))

        parser = yacc.yacc(debug=False, tabmodule='generic_parsetab',
                           outputdir=os.path.dirname(__file__))

        if not parser_exists:
            cls._add_tab_header('generic_parsetab')

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
    def _parse_unit(cls, s, detailed_exception=True):
        registry = core.get_current_unit_registry().registry
        if s == '%':
            return registry['percent']

        if not _is_ascii(s):
            if s[0] == '\N{MICRO SIGN}':
                s = 'u' + s[1:]
            if s[-1] == '\N{GREEK CAPITAL LETTER OMEGA}':
                s = s[:-1] + 'Ohm'
            elif s[-1] == '\N{LATIN CAPITAL LETTER A WITH RING ABOVE}':
                s = s[:-1] + 'Angstrom'

        if s in registry:
            return registry[s]

        if detailed_exception:
            raise ValueError(
                '{} is not a valid unit. {}'.format(
                    s, did_you_mean(s, registry)))
        else:
            raise ValueError()

    _translations = str.maketrans({
        '\N{GREEK SMALL LETTER MU}': '\N{MICRO SIGN}',
        '\N{MINUS SIGN}': '-',
    })
    """Character translations that should be applied before parsing a string.

    Note that this does explicitly *not* generally translate MICRO SIGN to u,
    since then a string like 'µ' would be interpreted as unit mass.
    """

    _superscripts = (
        '\N{SUPERSCRIPT MINUS}'
        '\N{SUPERSCRIPT PLUS SIGN}'
        '\N{SUPERSCRIPT ZERO}'
        '\N{SUPERSCRIPT ONE}'
        '\N{SUPERSCRIPT TWO}'
        '\N{SUPERSCRIPT THREE}'
        '\N{SUPERSCRIPT FOUR}'
        '\N{SUPERSCRIPT FIVE}'
        '\N{SUPERSCRIPT SIX}'
        '\N{SUPERSCRIPT SEVEN}'
        '\N{SUPERSCRIPT EIGHT}'
        '\N{SUPERSCRIPT NINE}'
    )

    _superscript_translations = str.maketrans(_superscripts, '-+0123456789')
    _regex_superscript = re.compile(f'[{_superscripts}]+')
    _regex_deg = re.compile('°([CF])?')

    @classmethod
    def _convert_superscript(cls, m):
        return '({})'.format(
            m.group().translate(cls._superscript_translations)
        )

    @classmethod
    def _convert_deg(cls, m):
        if len(m.string) == 1:
            return 'deg'
        return m.string.replace('°', 'deg_')

    @classmethod
    def parse(cls, s, debug=False):
        if not isinstance(s, str):
            s = s.decode('ascii')
        elif not _is_ascii(s):
            # common normalization of unicode strings to avoid
            # having to deal with multiple representations of
            # the same character. This normalizes to "composed" form
            # and will e.g. convert OHM SIGN to GREEK CAPITAL LETTER OMEGA
            s = unicodedata.normalize('NFC', s)
            # Translate some basic unicode items that we'd like to support on
            # input but are not standard.
            s = s.translate(cls._translations)

            # TODO: might the below be better done in the parser/lexer?
            # Translate superscripts to parenthesized numbers; this ensures
            # that mixes of superscripts and regular numbers fail.
            s = cls._regex_superscript.sub(cls._convert_superscript, s)
            # Translate possible degrees.
            s = cls._regex_deg.sub(cls._convert_deg, s)

        result = cls._do_parse(s, debug=debug)
        # Check for excess solidi, but exclude fractional exponents (accepted)
        n_slashes = s.count('/')
        if n_slashes > 1 and (n_slashes - len(re.findall(r'\(\d+/\d+\)', s))) > 1:
            warnings.warn(
                "'{}' contains multiple slashes, which is "
                "discouraged by the FITS standard".format(s),
                core.UnitsWarning)
        return result

    @classmethod
    def _do_parse(cls, s, debug=False):
        try:
            # This is a short circuit for the case where the string
            # is just a single unit name
            return cls._parse_unit(s, detailed_exception=False)
        except ValueError as e:
            try:
                return cls._parser.parse(s, lexer=cls._lexer, debug=debug)
            except ValueError as e:
                if str(e):
                    raise
                else:
                    raise ValueError(f"Syntax error parsing unit '{s}'")

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
                if '/' in power or '.' in power:
                    out.append('{}({})'.format(
                        cls._get_unit_name(base), power))
                else:
                    out.append('{}{}'.format(
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
