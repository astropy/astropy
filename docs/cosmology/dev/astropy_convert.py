# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see Astropy LICENSE.rst

"""
Register conversion methods for cosmology objects with Astropy Cosmology.

With this registered, we can start with a Cosmology from
``mypackage`` and convert it an astropy Cosmology instance.

    >>> from mypackage.cosmology import somecosmologyobject
    >>> from astropy.cosmology import Cosmology
    >>> cosmo = Cosmology.from_format(somecosmologyobject, format="mypackage")
    >>> cosmo

We can also do the reverse: start with an astropy Cosmology and convert it
to a ``mypackage`` object.

    >>> from astropy.cosmology import Planck18
    >>> myplanck = Planck18.to_format("mypackage")
    >>> myplanck

"""

# THIRD PARTY
from astropy.cosmology import Cosmology
from astropy.io import registry as io_registry

# LOCAL
from mypackage.cosmology import MyCosmology


def from_mypackage(mycosmo):
    """Load `~astropy.cosmology.Cosmology` from ``mypackage`` object."""
    # Cosmology provides a nice method "mapping", so all that needs to
    # be done here is create a dictionary of the parameters
    mapping = {}
    mapping["H0"] = mycosmo.hubble_parameter
    mapping["Om0"] = mycosmo.Omega_matter_initial
    ...  # keep building mapping

    return Cosmology.from_format(
        mapping, format="mapping", move_to_meta=True
    )  # extra info -> meta


def to_myformat(cosmology, file, *, overwrite=False, **kwargs):
    """Return the cosmology as a ``mycosmo``."""
    # Cosmology provides a nice method "mapping", so all that needs to
    # be done here is initialize from the dictionary
    mapping = cosmology.to_format("mapping")
    # correct entries
    mapping["hubble_parameter"] = mapping.pop("H0")
    mapping["Omega_matter_initial"] = mapping.pop("Om0")
    ...  # keep remapping

    return MyCosmology(**mapping)


def mypackage_identify(origin, format, *args, **kwargs):
    itis = False
    if origin == "read":
        itis = isinstance(args[1], MyCosmology) and (format in (None, "mypackage"))

    return itis


# -------------------------------------------------------------------
# Register to/from_format & identify methods with Astropy Unified I/O

io_registry.register_reader("mypackage", Cosmology, read_mypackage)
io_registry.register_writer("mypackage", Cosmology, write_mypackage)
io_registry.register_identifier("mypackage", Cosmology, mypackage_identify)
