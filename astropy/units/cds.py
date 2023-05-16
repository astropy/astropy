# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This package defines units used in the CDS format, both the units
defined in `Centre de Données astronomiques de Strasbourg
<https://cds.unistra.fr/>`_ `Standards for Astronomical Catalogues 2.0
<https://vizier.unistra.fr/vizier/doc/catstd-3.2.htx>`_ format and the `complete
set of supported units <https://vizier.unistra.fr/viz-bin/Unit>`_.
This format is used by VOTable up to version 1.2.

These units are not available in the top-level `astropy.units`
namespace.  To use these units, you must import the `astropy.units.cds`
module::

    >>> from astropy.units import cds
    >>> q = 10. * cds.lyr  # doctest: +SKIP

To include them in `~astropy.units.UnitBase.compose` and the results of
`~astropy.units.UnitBase.find_equivalent_units`, do::

    >>> from astropy.units import cds
    >>> cds.enable()  # doctest: +SKIP
"""

_ns = globals()


def _initialize_module():
    """Initialize CDS units module."""
    # Local imports to avoid polluting top-level namespace
    import numpy as np

    from astropy import units as u
    from astropy.constants import si as _si

    from . import core

    # The CDS format also supports power-of-2 prefixes as defined here:
    # http://physics.nist.gov/cuu/Units/binary.html
    prefixes = core.si_prefixes + core.binary_prefixes

    # CDS only uses the short prefixes
    prefixes = [(short, short, factor) for (short, long, factor) in prefixes]

    # The following units are defined in alphabetical order, directly from
    # here: https://vizier.unistra.fr/viz-bin/Unit

    mapping = [
        (["A"], u.A, "Ampere"),
        (["a"], u.a, "year", ["P"]),
        (["a0"], _si.a0, "Bohr radius"),
        (["al"], u.lyr, "Light year", ["c", "d"]),
        (["lyr"], u.lyr, "Light year"),
        (["alpha"], _si.alpha, "Fine structure constant"),
        ((["AA", "Å"], ["Angstrom", "Angstroem"]), u.AA, "Angstrom"),
        (["arcmin", "arcm"], u.arcminute, "minute of arc"),
        (["arcsec", "arcs"], u.arcsecond, "second of arc"),
        (["atm"], _si.atm, "atmosphere"),
        (["AU", "au"], u.au, "astronomical unit"),
        (["bar"], u.bar, "bar"),
        (["barn"], u.barn, "barn"),
        (["bit"], u.bit, "bit"),
        (["byte"], u.byte, "byte"),
        (["C"], u.C, "Coulomb"),
        (["c"], _si.c, "speed of light", ["p"]),
        (["cal"], 4.1854 * u.J, "calorie"),
        (["cd"], u.cd, "candela"),
        (["ct"], u.ct, "count"),
        (["D"], u.D, "Debye (dipole)"),
        (["d"], u.d, "Julian day", ["c"]),
        ((["deg", "°"], ["degree"]), u.degree, "degree"),
        (["dyn"], u.dyn, "dyne"),
        (["e"], _si.e, "electron charge", ["m"]),
        (["eps0"], _si.eps0, "electric constant"),
        (["erg"], u.erg, "erg"),
        (["eV"], u.eV, "electron volt"),
        (["F"], u.F, "Farad"),
        (["G"], _si.G, "Gravitation constant"),
        (["g"], u.g, "gram"),
        (["gauss"], u.G, "Gauss"),
        (["geoMass", "Mgeo"], u.M_earth, "Earth mass"),
        (["H"], u.H, "Henry"),
        (["h"], u.h, "hour", ["p"]),
        (["hr"], u.h, "hour"),
        (["\\h"], _si.h, "Planck constant"),
        (["Hz"], u.Hz, "Hertz"),
        (["inch"], 0.0254 * u.m, "inch"),
        (["J"], u.J, "Joule"),
        (["JD"], u.d, "Julian day", ["M"]),
        (["jovMass", "Mjup"], u.M_jup, "Jupiter mass"),
        (["Jy"], u.Jy, "Jansky"),
        (["K"], u.K, "Kelvin"),
        (["k"], _si.k_B, "Boltzmann"),
        (["l"], u.l, "litre", ["a"]),
        (["lm"], u.lm, "lumen"),
        (["Lsun", "solLum"], u.solLum, "solar luminosity"),
        (["lx"], u.lx, "lux"),
        (["m"], u.m, "meter"),
        (["mag"], u.mag, "magnitude"),
        (["me"], _si.m_e, "electron mass"),
        (["min"], u.minute, "minute"),
        (["MJD"], u.d, "Julian day"),
        (["mmHg"], 133.322387415 * u.Pa, "millimeter of mercury"),
        (["mol"], u.mol, "mole"),
        (["mp"], _si.m_p, "proton mass"),
        (["Msun", "solMass"], u.solMass, "solar mass"),
        ((["mu0", "µ0"], []), _si.mu0, "magnetic constant"),
        (["muB"], _si.muB, "Bohr magneton"),
        (["N"], u.N, "Newton"),
        (["Ohm"], u.Ohm, "Ohm"),
        (["Pa"], u.Pa, "Pascal"),
        (["pc"], u.pc, "parsec"),
        (["ph"], u.ph, "photon"),
        (["pi"], u.Unit(np.pi), "π"),
        (["pix"], u.pix, "pixel"),
        (["ppm"], u.Unit(1e-6), "parts per million"),
        (["R"], _si.R, "gas constant"),
        (["rad"], u.radian, "radian"),
        (["Rgeo"], _si.R_earth, "Earth equatorial radius"),
        (["Rjup"], _si.R_jup, "Jupiter equatorial radius"),
        (["Rsun", "solRad"], u.solRad, "solar radius"),
        (["Ry"], u.Ry, "Rydberg"),
        (["S"], u.S, "Siemens"),
        (["s", "sec"], u.s, "second"),
        (["sr"], u.sr, "steradian"),
        (["Sun"], u.Sun, "solar unit"),
        (["T"], u.T, "Tesla"),
        (["t"], 1e3 * u.kg, "metric tonne", ["c"]),
        (["u"], _si.u, "atomic mass", ["da", "a"]),
        (["V"], u.V, "Volt"),
        (["W"], u.W, "Watt"),
        (["Wb"], u.Wb, "Weber"),
        (["yr"], u.a, "year"),
    ]

    for entry in mapping:
        if len(entry) == 3:
            names, unit, doc = entry
            excludes = []
        else:
            names, unit, doc, excludes = entry
        core.def_unit(
            names,
            unit,
            prefixes=prefixes,
            namespace=_ns,
            doc=doc,
            exclude_prefixes=excludes,
        )

    core.def_unit(["µas"], u.microarcsecond, doc="microsecond of arc", namespace=_ns)
    core.def_unit(["mas"], u.milliarcsecond, doc="millisecond of arc", namespace=_ns)
    core.def_unit(
        ["---", "-"],
        u.dimensionless_unscaled,
        doc="dimensionless and unscaled",
        namespace=_ns,
    )
    core.def_unit(["%"], u.percent, doc="percent", namespace=_ns)
    # The Vizier "standard" defines this in units of "kg s-3", but
    # that may not make a whole lot of sense, so here we just define
    # it as its own new disconnected unit.
    core.def_unit(["Crab"], prefixes=prefixes, namespace=_ns, doc="Crab (X-ray) flux")


_initialize_module()


###########################################################################
# DOCSTRING

if __doc__ is not None:
    # This generates a docstring for this module that describes all of the
    # standard units defined here.
    from .utils import generate_unit_summary as _generate_unit_summary

    __doc__ += _generate_unit_summary(globals())


def enable():
    """
    Enable CDS units so they appear in results of
    `~astropy.units.UnitBase.find_equivalent_units` and
    `~astropy.units.UnitBase.compose`.  This will disable
    all of the "default" `astropy.units` units, since there
    are some namespace clashes between the two.

    This may be used with the ``with`` statement to enable CDS
    units only temporarily.
    """
    # Local imports to avoid cyclical import and polluting namespace
    import inspect

    from .core import set_enabled_units

    return set_enabled_units(inspect.getmodule(enable))
