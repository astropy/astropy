# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICNSE.rst
"""
Handles a "generic" string format for units
"""
from __future__ import absolute_import, division, print_function, unicode_literals

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
            CDS._parser = self._make_parser()

    @staticmethod
    def _generate_unit_names():
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
        from ...extern import pyparsing as p

        product = p.Literal(".")
        division = p.Literal("/")
        open_p = p.Literal("(")
        close_p = p.Literal(")")
        literal10 = p.Literal("10")
        literalx = p.Literal("x")

        unsigned_integer = p.Regex(r'\d+')
        signed_integer = p.Regex(r'[+-]\d+')
        integer = p.Regex(r'[+-]?\d+')
        floating_point = p.Regex(r'[+-]?((\d+\.?\d*)|(\.\d+))([eE][+-]?\d+)?')

        factor = p.Forward()
        main = p.Forward()
        numeric_power = p.Forward()
        product_of_units = p.Forward()
        unit_expression = p.Forward()
        unit = p.Forward()
        unit_with_power = p.Forward()

        main << (
            (p.Optional(factor, default=1.0) +
             product_of_units +
             p.StringEnd()))

        product_of_units << (
            (unit_expression + p.StringEnd()) ^
            (division + unit_expression) ^
            (unit_expression + product + product_of_units) ^
            (unit_expression + division + product_of_units))

        unit_expression << (
            (unit_with_power) ^
            (p.Suppress(open_p) + product_of_units + p.Suppress(close_p)))

        factor << (
            (floating_point + literalx + literal10 + signed_integer) ^
            (literal10 + signed_integer) ^
            (unsigned_integer) ^
            (literal10) ^
            (floating_point))

        unit_with_power << (
            (unit + p.Optional(numeric_power)))

        numeric_power << (
            integer)

        unit << (
            p.Literal('%') |
            p.Word(p.alphas, p.alphas + '_'))

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
    @utils._trace
    def _parse_unsigned_integer(cls, s, loc, toks):
        return int(toks[0])

    @classmethod
    @utils._trace
    def _parse_signed_integer(cls, s, loc, toks):
        return int(toks[0])

    @classmethod
    @utils._trace
    def _parse_integer(cls, s, loc, toks):
        return int(toks[0])

    @classmethod
    @utils._trace
    def _parse_floating_point(cls, s, loc, toks):
        return float(toks[0])

    @classmethod
    @utils._trace
    def _parse_factor(cls, s, loc, toks):
        if len(toks) == 1:
            return toks[0]
        elif len(toks) == 2:
            return 10.0 ** float(toks[1])
        elif len(toks) == 4:
            return toks[0] * 10.0 ** toks[3]

    @classmethod
    @utils._trace
    def _parse_unit_with_power(cls, s, loc, toks):
        if len(toks) == 1:
            return toks[0]
        elif len(toks) == 2:
            return toks[0] ** toks[1]

    @classmethod
    @utils._trace
    def _parse_unit(cls, s, loc, toks):
        from ...extern import pyparsing as p

        unit = toks[0]
        if unit not in cls._units:
            raise p.ParseException(
                "Unit {0!r} not supported by the CDS SAC "
                "standard.".format(unit))

        return cls._units[unit]

    @classmethod
    @utils._trace
    def _parse_product_of_units(cls, s, loc, toks):
        if len(toks) == 1:
            return toks[0]
        elif len(toks) == 2:
            from ..core import Unit
            return Unit(1.0 / toks[1])
        elif len(toks) == 3:
            if toks[1] == '/':
                return toks[0] / toks[2]
            else:
                return toks[0] * toks[2]

    @classmethod
    @utils._trace
    def _parse_main(cls, s, loc, toks):
        from ..core import Unit
        return Unit(toks[0] * toks[1])

    def parse(self, s):
        from ...extern import pyparsing as p

        if utils.DEBUG:
            print("parse", s)

        if ' ' in s:
            raise ValueError('CDS unit must not contain whitespace')

        # This is a short circuit for the case where the string
        # is just a single unit name
        try:
            return self._parse_unit(s, 0, [s])
        except p.ParseException as e:
            try:
                return self._parser.parseString(s, parseAll=True)[0]
            except p.ParseException as e:
                raise ValueError("{0} in {1!r}".format(
                    utils.cleanup_pyparsing_error(e), s))

    def _get_unit_name(self, unit):
        return unit.get_format_name('cds')

    def _format_unit_list(self, units):
        out = []
        for base, power in units:
            if power == 1:
                out.append(self._get_unit_name(base))
            else:
                out.append('{0}{1}'.format(
                    self._get_unit_name(base), power))
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
