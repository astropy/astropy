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

"""Handles the CDS string format for units."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from astropy.units.core import CompositeUnit, Unit
from astropy.units.utils import is_effectively_unity
from astropy.utils import classproperty, parsing

from .base import Base, _ParsingFormatMixin

if TYPE_CHECKING:
    from typing import ClassVar, Literal

    from astropy.extern.ply.lex import Lexer
    from astropy.units import UnitBase
    from astropy.utils.parsing import ThreadSafeParser


class CDS(Base, _ParsingFormatMixin):
    """
    Support the `Centre de Données astronomiques de Strasbourg
    <https://cds.unistra.fr/>`_ `Standards for Astronomical
    Catalogues 2.0 <https://vizier.unistra.fr/vizier/doc/catstd-3.2.htx>`_
    format, and the `complete set of supported units
    <https://vizier.unistra.fr/viz-bin/Unit>`_.  This format is used
    by VOTable up to version 1.2.
    """

    _space: ClassVar[str] = "."
    _times: ClassVar[str] = "x"
    _scale_unit_separator: ClassVar[str] = ""

    _tokens: ClassVar[tuple[str, ...]] = (
        "PRODUCT",
        "DIVISION",
        "OPEN_PAREN",
        "CLOSE_PAREN",
        "OPEN_BRACKET",
        "CLOSE_BRACKET",
        "X",
        "SIGN",
        "UINT",
        "UFLOAT",
        "UNIT",
        "DIMENSIONLESS",
    )

    @classproperty(lazy=True)
    def _units(cls) -> dict[str, UnitBase]:
        from astropy import units as u
        from astropy.units import cds

        return {k: v for k, v in cds.__dict__.items() if isinstance(v, u.UnitBase)}

    @classproperty(lazy=True)
    def _lexer(cls) -> Lexer:
        tokens = cls._tokens

        t_PRODUCT = r"\."
        t_DIVISION = r"/"
        t_OPEN_PAREN = r"\("
        t_CLOSE_PAREN = r"\)"
        t_OPEN_BRACKET = r"\["
        t_CLOSE_BRACKET = r"\]"

        # NOTE THE ORDERING OF THESE RULES IS IMPORTANT!!
        # Regular expression rules for simple tokens

        def t_UFLOAT(t):
            r"((\d+\.?\d+)|(\.\d+))([eE][+-]?\d+)?"
            if not re.search(r"[eE\.]", t.value):
                t.type = "UINT"
                t.value = int(t.value)
            else:
                t.value = float(t.value)
            return t

        def t_UINT(t):
            r"\d+"
            t.value = int(t.value)
            return t

        def t_SIGN(t):
            r"[+-](?=\d)"
            t.value = float(t.value + "1")
            return t

        def t_X(t):  # multiplication for factor in front of unit
            r"[x×]"
            return t

        # Most units are just combinations of letters with no numbers, but there
        # are a few special ones (\h is Planch constant) and three that end in 0.
        def t_UNIT(t):
            r"%|°|\\h|(a|eps|mu)0|((?!\d)\w)+"
            t.value = cls._get_unit(t)
            return t

        def t_DIMENSIONLESS(t):
            r"---|-"
            # These are separate from t_UNIT since they cannot have a prefactor.
            t.value = cls._get_unit(t)
            return t

        t_ignore = ""

        # Error handling rule
        def t_error(t):
            raise ValueError(f"Invalid character at col {t.lexpos}")

        return parsing.lex(
            lextab="cds_lextab", package="astropy/units", reflags=int(re.UNICODE)
        )

    @classproperty(lazy=True)
    def _parser(cls) -> ThreadSafeParser:
        """
        The grammar here is based on the description in the `Standards
        for Astronomical Catalogues 2.0
        <https://vizier.unistra.fr/vizier/doc/catstd-3.2.htx>`_, which is not
        terribly precise.  The exact grammar is here is based on the
        YACC grammar in the `unity library <https://purl.org/nxg/dist/unity/>`_.
        """
        tokens = cls._tokens

        def p_main(p):
            """
            main : factor combined_units
                 | combined_units
                 | DIMENSIONLESS
                 | OPEN_BRACKET combined_units CLOSE_BRACKET
                 | OPEN_BRACKET DIMENSIONLESS CLOSE_BRACKET
                 | factor
            """
            from astropy.units import dex

            if len(p) == 3:
                p[0] = CompositeUnit(p[1] * p[2].scale, p[2].bases, p[2].powers)
            elif len(p) == 4:
                p[0] = dex(p[2])
            else:
                p[0] = Unit(p[1])

        def p_combined_units(p):
            """
            combined_units : product_of_units
                           | division_of_units
            """
            p[0] = p[1]

        def p_product_of_units(p):
            """
            product_of_units : unit_expression PRODUCT combined_units
                             | unit_expression
            """
            if len(p) == 4:
                p[0] = p[1] * p[3]
            else:
                p[0] = p[1]

        def p_division_of_units(p):
            """
            division_of_units : DIVISION unit_expression
                              | combined_units DIVISION unit_expression
            """
            if len(p) == 3:
                p[0] = p[2] ** -1
            else:
                p[0] = p[1] / p[3]

        def p_unit_expression(p):
            """
            unit_expression : unit_with_power
                            | OPEN_PAREN combined_units CLOSE_PAREN
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = p[2]

        def p_factor(p):
            """
            factor : signed_float X UINT signed_int
                   | UINT X UINT signed_int
                   | UINT signed_int
                   | UINT
                   | signed_float
            """
            if len(p) == 5:
                if p[3] != 10:
                    raise ValueError("Only base ten exponents are allowed in CDS")
                p[0] = p[1] * 10.0 ** p[4]
            elif len(p) == 3:
                if p[1] != 10:
                    raise ValueError("Only base ten exponents are allowed in CDS")
                p[0] = 10.0 ** p[2]
            elif len(p) == 2:
                p[0] = p[1]

        def p_unit_with_power(p):
            """
            unit_with_power : UNIT numeric_power
                            | UNIT
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = p[1] ** p[2]

        def p_numeric_power(p):
            """
            numeric_power : sign UINT
            """
            p[0] = p[1] * p[2]

        def p_sign(p):
            """
            sign : SIGN
                 |
            """
            if len(p) == 2:
                p[0] = p[1]
            else:
                p[0] = 1.0

        def p_signed_int(p):
            """
            signed_int : SIGN UINT
            """
            p[0] = p[1] * p[2]

        def p_signed_float(p):
            """
            signed_float : sign UINT
                         | sign UFLOAT
            """
            p[0] = p[1] * p[2]

        def p_error(p):
            raise ValueError()

        return parsing.yacc(tabmodule="cds_parsetab", package="astropy/units")

    @classmethod
    def parse(cls, s: str, debug: bool = False) -> UnitBase:
        if " " in s:
            raise ValueError("CDS unit must not contain whitespace")
        if not isinstance(s, str):
            s = s.decode("ascii")

        return cls._do_parse(s, debug)

    @classmethod
    def _format_mantissa(cls, m: str) -> str:
        return "" if m == "1" else m

    @classmethod
    def _format_superscript(cls, number: str) -> str:
        return number if number.startswith("-") else "+" + number

    @classmethod
    def to_string(
        cls, unit: UnitBase, fraction: bool | Literal["inline", "multiline"] = False
    ) -> str:
        # Remove units that aren't known to the format
        unit = cls._decompose_to_known_units(unit)

        if not unit.bases:
            if unit.scale == 1:
                return "---"
            elif is_effectively_unity(unit.scale * 100.0):
                return "%"

        return super().to_string(unit, fraction=fraction)
