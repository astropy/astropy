# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package defines units that can also be used as functions of other units.
If called, their arguments are used to initialize the corresponding function
unit (e.g., ``u.mag(u.ct/u.s)``).  Note that the prefixed versions cannot be
called, as it would be unclear what, e.g., ``u.mmag(u.ct/u.s)`` would mean.

It also defines a few commonly used magnitude unit instances, like STmag,
for which the physical units are defined in `astropy.units.photometric`.

All units are also available in (and should be used through) the
`astropy.units` namespace.
"""

from astropy.units import photometric
from astropy.units.core import CompositeUnit, UnitBase, _add_prefixes

from .logarithmic import DecibelUnit, DexUnit, MagUnit
from .mixin import IrreducibleFunctionUnit, RegularFunctionUnit

__all__: list[str] = []  #  Units are added at the end

_ns = globals()

###########################################################################
# Logarithmic units

# These calls are what core.def_unit would do, but we need to use the callable
# unit versions.  The actual function unit classes get added in logarithmic.

dex = IrreducibleFunctionUnit(
    ["dex"], namespace=_ns, doc="Dex: Base 10 logarithmic unit"
)
dex._function_unit_class = DexUnit

dB = RegularFunctionUnit(
    ["dB", "decibel"],
    0.1 * dex,
    namespace=_ns,
    doc="Decibel: ten per base 10 logarithmic unit",
)
dB._function_unit_class = DecibelUnit

mag = RegularFunctionUnit(
    ["mag"],
    -0.4 * dex,
    namespace=_ns,
    doc="Astronomical magnitude: -2.5 per base 10 logarithmic unit",
)
_add_prefixes(mag, namespace=_ns, prefixes=True)
mag._function_unit_class = MagUnit

STmag = mag(photometric.STflux)
STmag.__doc__ = "ST magnitude: STmag=-21.1 corresponds to 1 erg/s/cm2/A"

ABmag = mag(photometric.ABflux)
ABmag.__doc__ = "AB magnitude: ABmag=-48.6 corresponds to 1 erg/s/cm2/Hz"

M_bol = mag(photometric.Bol)
M_bol.__doc__ = (
    f"Absolute bolometric magnitude: M_bol=0 corresponds to {photometric.Bol}"
)

m_bol = mag(photometric.bol)
m_bol.__doc__ = (
    f"Apparent bolometric magnitude: m_bol=0 corresponds to {photometric.bol}"
)


###########################################################################
# DOCSTRING

__all__ += [n for n, v in _ns.items() if isinstance(v, (UnitBase, MagUnit))]

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    from astropy.units.utils import generate_unit_summary as _generate_unit_summary

    def _description(unit):
        pu = unit.physical_unit.represents
        if unit.__doc__[:2] in {"AB", "ST"}:
            pu = 1.0 * CompositeUnit(1.0, pu.bases, pu.powers)
        return "".join(
            unit.__doc__.partition("corresponds to ")[:-1]
            + (f":math:`{pu.to_string(format='latex')[1:-1]}`",)
        )

    template = """
   * - ``{}``
     - {}
     - {}
"""

    __doc__ += (
        _generate_unit_summary(globals())
        + """
.. list-table:: Available Magnitude Units
   :header-rows: 1
   :widths: 10 50 10

   * - Unit
     - Description
     - Represents
"""
        + "".join(
            [
                template.format(key, _description(val), val)
                for key, val in globals().items()
                if isinstance(val, MagUnit)
            ]
        )
    )
