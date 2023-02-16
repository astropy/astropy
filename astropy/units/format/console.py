# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Handles the "Console" unit format.
"""


from . import base, core, utils


class Console(base.Base):
    """
    Output-only format for to display pretty formatting at the
    console.

    For example::

      >>> import astropy.units as u
      >>> print(u.Ry.decompose().to_string('console'))  # doctest: +FLOAT_CMP
      2.1798721*10^-18m^2 kg s^-2
      >>> print(u.Ry.decompose().to_string('console', inline=False))  # doctest: +FLOAT_CMP
                       m^2 kg
      2.1798721*10^-18 ------
                        s^2
    """

    _times = "*"
    _line = "-"

    @classmethod
    def _get_unit_name(cls, unit):
        return unit.get_format_name("console")

    @classmethod
    def _format_superscript(cls, number):
        return f"^{number}"

    @classmethod
    def _format_unit_list(cls, units):
        out = []
        for base_, power in units:
            if power == 1:
                out.append(cls._get_unit_name(base_))
            else:
                out.append(
                    cls._get_unit_name(base_)
                    + cls._format_superscript(utils.format_power(power))
                )
        return " ".join(out)

    @classmethod
    def format_exponential_notation(cls, val):
        m, ex = utils.split_mantissa_exponent(val)

        parts = []
        if m:
            parts.append(m)

        if ex:
            parts.append(f"10{cls._format_superscript(ex)}")

        return cls._times.join(parts)

    @classmethod
    def to_string(cls, unit, inline=True):
        if isinstance(unit, core.CompositeUnit):
            if unit.scale == 1:
                s = ""
            else:
                s = cls.format_exponential_notation(unit.scale)

            if len(unit.bases):
                if inline:
                    nominator = zip(unit.bases, unit.powers)
                    denominator = []
                else:
                    nominator, denominator = utils.get_grouped_by_powers(
                        unit.bases, unit.powers
                    )
                if len(denominator):
                    if len(nominator):
                        nominator = cls._format_unit_list(nominator)
                    else:
                        nominator = "1"
                    denominator = cls._format_unit_list(denominator)
                    fraclength = max(len(nominator), len(denominator))
                    f = f"{{0:^{len(s)}s}} {{1:^{fraclength}s}}"

                    lines = [
                        f.format("", nominator),
                        f.format(s, cls._line * fraclength),
                        f.format("", denominator),
                    ]

                    s = "\n".join(lines)
                else:
                    nominator = cls._format_unit_list(nominator)
                    s += nominator
        elif isinstance(unit, core.NamedUnit):
            s = cls._get_unit_name(unit)

        return s
