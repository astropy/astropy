# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handles a "generic" string format for units
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from fractions import Fraction
import functools

from .base import Base
from . import utils


DEBUG = False


def _trace(func):
    """
    A utility decorator to help debug the parser.
    """
    def run(self, s, loc, toks):
        if DEBUG:
            print(func.__name__, toks, end=' ')
            try:
                result = func(self, s, loc, toks)
            except Exception as e:
                print("Exception: ", e.message)
                raise
            print(result)
        else:
            return func(self, s, loc, toks)
        return result

    return functools.update_wrapper(run, func)


class Generic(Base):
    """
    A "generic" format.

    The syntax of the format is based directly on the FITS standard,
    but instead of only supporting the units that FITS knows about, it
    supports any unit available in the `astropy.units` namespace.
    """
    def __init__(self):
        # Build this on the class, so it only gets generated once.
        if '_parser' not in Generic.__dict__:
            Generic._parser = self._make_parser()

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
        from astropy.extern import pyparsing as p

        product = p.Literal("*") | p.Literal(".") | p.White()
        division = p.Literal("/")
        power = p.Literal("**") | p.Literal("^") | p.Empty()
        open_p = p.Literal("(")
        close_p = p.Literal(")")
        # TODO: We only support the sqrt function for now because it's
        # obvious how to handle it.
        function_name = p.Literal("sqrt")

        unsigned_integer = p.Regex(r'\d+')
        signed_integer = p.Regex(r'[+-]\d+')
        integer = p.Regex(r'[+-]?\d+')
        floating_point = p.Regex(r'[+-]?((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+)?')

        division_product_of_units = p.Forward()
        factor = p.Forward()
        factor_product_of_units = p.Forward()
        frac = p.Forward()
        function = p.Forward()
        main = p.Forward()
        numeric_power = p.Forward()
        product_of_units = p.Forward()
        unit_expression = p.Forward()
        unit = p.Forward()
        unit_with_power = p.Forward()

        main << (
            (factor_product_of_units) |
            (division_product_of_units))

        factor_product_of_units << (
            (p.Optional(factor +
                        p.Suppress(p.Optional(p.White())), default=1.0) +
             product_of_units +
             p.StringEnd()))

        division_product_of_units << (
            (p.Optional(factor, default=1.0) +
             p.Optional(product_of_units, default=1.0) +
             p.Suppress(division) +
             unit_expression +
             p.StringEnd()))

        product_of_units << (
            (unit_expression + p.Suppress(product) + product_of_units) |
            (unit_expression))

        function << (
            function_name +
            p.Suppress(open_p) + unit_expression + p.Suppress(close_p))

        unit_expression << (
            (function) |
            (unit_with_power) |
            (p.Suppress(open_p) + product_of_units + p.Suppress(close_p))
            )

        factor << (
            (unsigned_integer + signed_integer) |
            (unsigned_integer + p.Suppress(power) + numeric_power) |
            (floating_point + p.Suppress(p.White()) +
             unsigned_integer + signed_integer) |
            (floating_point + p.Suppress(p.White()) +
             unsigned_integer + p.Suppress(power) + numeric_power) |
            (floating_point)
            )

        unit << p.Word(p.alphas, p.alphas + '_')

        unit_with_power << (
            (unit + p.Suppress(power) + numeric_power) |
            (unit))

        numeric_power << (
            integer |
            (p.Suppress(open_p) + integer + p.Suppress(close_p)) |
            (p.Suppress(open_p) + floating_point + p.Suppress(close_p)) |
            (p.Suppress(open_p) + frac + p.Suppress(close_p)))

        frac << (
            integer + p.Suppress(division) + integer)

        # Set actions
        for key, val in locals().items():
            if isinstance(val, p.ParserElement):
                val.setName(key)
                val.leaveWhitespace()
            method_name = "_parse_{0}".format(key)
            if hasattr(cls, method_name):
                val.setParseAction(getattr(cls, method_name))

        return main

    @classmethod
    @_trace
    def _parse_unsigned_integer(cls, s, loc, toks):
        return int(toks[0])

    @classmethod
    @_trace
    def _parse_signed_integer(cls, s, loc, toks):
        return int(toks[0])

    @classmethod
    @_trace
    def _parse_integer(cls, s, loc, toks):
        return int(toks[0])

    @classmethod
    @_trace
    def _parse_floating_point(cls, s, loc, toks):
        return float(toks[0])

    @classmethod
    @_trace
    def _parse_factor(cls, s, loc, toks):
        if len(toks) == 1:
            return toks[0]
        elif len(toks) == 2:
            return toks[0] ** float(toks[1])
        elif len(toks) == 3:
            return float(toks[0]) * toks[1] ** float(toks[2])

    @classmethod
    @_trace
    def _parse_frac(cls, s, loc, toks):
        return Fraction(toks[0], toks[1])

    @classmethod
    @_trace
    def _parse_unit(cls, s, loc, toks):
        if toks[0] in cls._unit_namespace:
            return cls._unit_namespace[toks[0]]
        raise ValueError(
            s, loc, "{0!r} is not a recognized unit".format(toks[0]))

    @classmethod
    @_trace
    def _parse_product_of_units(cls, s, loc, toks):
        if len(toks) == 1:
            return toks[0]
        else:
            return toks[0] * toks[1]

    @classmethod
    @_trace
    def _parse_division_product_of_units(cls, s, loc, toks):
        return (toks[0] * toks[1]) / toks[2]

    @classmethod
    @_trace
    def _parse_factor_product_of_units(cls, s, loc, toks):
        if toks[0] != 1.0:
            return toks[0] * toks[1]
        else:
            return toks[1]

    @classmethod
    @_trace
    def _parse_unit_with_power(cls, s, loc, toks):
        if len(toks) == 1:
            return toks[0]
        else:
            return toks[0] ** toks[1]

    @classmethod
    @_trace
    def _parse_function(cls, s, loc, toks):
        # TODO: Add support for more functions here
        if toks[0] == 'sqrt':
            return toks[1] ** -2.0
        else:
            raise ValueError(
                s, loc, "{0!r} is not a recognized function".format(
                    toks[0]))

    def parse(self, s):
        if DEBUG:
            print("parse", s)

        if '_unit_namespace' not in Generic.__dict__:
            from ... import units as u
            ns = {}
            for key, val in u.__dict__.items():
                if isinstance(val, u.UnitBase):
                   ns[key] = val
            Generic._unit_namespace = ns

        return self._parser.parseString(s, parseAll=True)[0]

    def _get_unit_name(self, unit):
        return unit.get_format_name('generic')

    def _format_unit_list(self, units):
        out = []
        for base, power in units:
            if power == 1:
                out.append(self._get_unit_name(base))
            else:
                out.append('{0}{1}'.format(
                    self._get_unit_name(base), power))
        return ' '.join(out)

    def to_string(self, unit):
        from .. import core

        if isinstance(unit, core.CompositeUnit):
            if unit.scale != 1:
                s = '{0:e} '.format(unit.scale)
            else:
                s = ''

            if len(unit.bases):
                positives, negatives = utils.get_grouped_by_powers(
                    unit.bases, unit.powers)
                if len(positives):
                    s += self._format_unit_list(positives)
                elif s == '':
                    s = '1'

                if len(negatives):
                    s += ' / ({0})'.format(self._format_unit_list(negatives))
        elif isinstance(unit, core.NamedUnit):
            s = self._get_unit_name(unit)

        return s
