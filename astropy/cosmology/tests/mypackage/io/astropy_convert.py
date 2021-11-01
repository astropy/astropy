# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see Astropy LICENSE.rst

"""
Register conversion methods for cosmology objects with Astropy Cosmology.

With this registered, we can start with a Cosmology from
``mypackage`` and convert it to an astropy Cosmology instance.

    >>> from mypackage.cosmology import myplanck
    >>> from astropy.cosmology import Cosmology
    >>> cosmo = Cosmology.from_format(myplanck, format="mypackage")
    >>> cosmo

We can also do the reverse: start with an astropy Cosmology and convert it
to a ``mypackage`` object.

    >>> from astropy.cosmology import Planck18
    >>> myplanck = Planck18.to_format("mypackage")
    >>> myplanck

"""

# THIRD PARTY
import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import FLRW, Cosmology, FlatLambdaCDM
from astropy.cosmology.connect import convert_registry

# LOCAL
from mypackage.cosmology import MyCosmology

__doctest_skip__ = ['*']


def from_mypackage(mycosmo):
    """Load `~astropy.cosmology.Cosmology` from ``mypackage`` object.

    Parameters
    ----------
    mycosmo : `~mypackage.cosmology.MyCosmology`

    Returns
    -------
    `~astropy.cosmology.Cosmology`
    """
    m = dict(mycosmo)
    m["name"] = mycosmo.name

    # ----------------
    # remap Parameters
    m["H0"] = m.pop("hubble_parameter") * (u.km / u.s / u.Mpc)
    m["Om0"] = m.pop("initial_matter_density")
    m["Tcmb0"] = m.pop("initial_temperature") * u.K
    # m["Neff"] = m.pop("Neff")  # skip b/c unchanged
    m["m_nu"] = m.pop("neutrino_masses") * u.eV
    m["Ob0"] = m.pop("initial_baryon_density")

    # ----------------
    # remap metadata
    m["t0"] = m.pop("current_age") * u.Gyr

    # optional
    if "reionization_redshift" in m:
        m["z_reion"] = m.pop("reionization_redshift")

    # ...  # keep building `m`

    # ----------------
    # Detect which type of Astropy cosmology to build.
    # TODO! CUSTOMIZE FOR DETECTION
    # Here we just force FlatLambdaCDM, but if your package allows for
    # non-flat cosmologies...
    m["cosmology"] = FlatLambdaCDM

    # build cosmology
    return Cosmology.from_format(m, format="mapping", move_to_meta=True)


def to_mypackage(cosmology, *args):
    """Return the cosmology as a ``mycosmo``.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology`

    Returns
    -------
    `~mypackage.cosmology.MyCosmology`
    """
    if not isinstance(cosmology, FLRW):
        raise TypeError("format 'mypackage' only supports FLRW cosmologies.")

    # ----------------
    # Cosmology provides a nice method "mapping", so all that needs to
    # be done here is initialize from the dictionary
    m = cosmology.to_format("mapping")

    # Detect which type of MyCosmology to build.
    # Here we have forced FlatLambdaCDM, but if your package allows for
    # non-flat cosmologies...
    m.pop("cosmology")

    # MyCosmology doesn't support metadata. If your cosmology class does...
    meta = m.pop("meta")
    m = {**meta, **m}  # merge, preferring current values

    # ----------------
    # remap values
    # MyCosmology doesn't support units, so take values.
    m["hubble_parameter"] = m.pop("H0").to_value(u.km/u.s/u.Mpc)
    m["initial_matter_density"] = m.pop("Om0")
    m["initial_temperature"] = m.pop("Tcmb0").to_value(u.K)
    # m["Neff"] = m.pop("Neff")  # skip b/c unchanged
    m["neutrino_masses"] = m.pop("m_nu").to_value(u.eV)
    m["initial_baryon_density"] = m.pop("Ob0")
    m["current_age"] = m.pop("t0", cosmology.age(0 * cu.redshift)).to_value(u.Gyr)

    # optional
    if "z_reion" in m:
        m["reionization_redshift"] = (m.pop("z_reion") << cu.redshift).value

    # ...  # keep remapping

    return MyCosmology(**m)


def mypackage_identify(origin, format, *args, **kwargs):
    """Identify if object uses format "mypackage"."""
    itis = False
    if origin == "read":
        itis = isinstance(args[1], MyCosmology) and (format in (None, "mypackage"))
    return itis


# -------------------------------------------------------------------
# Register to/from_format & identify methods with Astropy Unified I/O

convert_registry.register_reader("mypackage", Cosmology, from_mypackage, force=True)
convert_registry.register_writer("mypackage", Cosmology, to_mypackage, force=True)
convert_registry.register_identifier("mypackage", Cosmology, mypackage_identify, force=True)
