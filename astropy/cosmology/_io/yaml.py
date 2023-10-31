# Licensed under a 3-clause BSD style license - see LICENSE.rst
r"""|Cosmology| <-> YAML I/O, using |Cosmology.to_format| and |Cosmology.from_format|.

This module provides functions to transform a |Cosmology| object to and from a `yaml
<https://yaml.org>`_ representation. The functions are registered with
``convert_registry`` under the format name "yaml". This format is primarily intended for
use by other I/O functions, e.g. |Table|'s metadata serialization, which themselves
require YAML serialization.

    >>> from astropy.cosmology import Planck18
    >>> yml = Planck18.to_format("yaml")
    >>> yml  # doctest: +NORMALIZE_WHITESPACE
    "!astropy.cosmology...FlatLambdaCDM\nH0: !astropy.units.Quantity...

    >>> print(Cosmology.from_format(yml, format="yaml"))
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                  Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
"""  # this is shown in the docs.

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology.connect import convert_registry
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Cosmology
from astropy.io.misc.yaml import AstropyDumper, AstropyLoader, dump, load

from .mapping import from_mapping
from .utils import FULLQUALNAME_SUBSTITUTIONS as QNS

__all__ = []  # nothing is publicly scoped


##############################################################################
# Serializer Functions
# these do Cosmology <-> YAML through a modified dictionary representation of
# the Cosmology object. The Unified-I/O functions are just wrappers to the YAML
# that calls these functions.


def yaml_representer(tag):
    """`yaml <https://yaml.org>`_ representation of |Cosmology| object.

    Parameters
    ----------
    tag : str
        The class tag, e.g. '!astropy.cosmology.LambdaCDM'

    Returns
    -------
    representer : callable[[`~astropy.io.misc.yaml.AstropyDumper`, |Cosmology|], str]
        Function to construct :mod:`yaml` representation of |Cosmology| object.
    """

    def representer(dumper, obj):
        """Cosmology yaml representer function for {}.

        Parameters
        ----------
        dumper : `~astropy.io.misc.yaml.AstropyDumper`
        obj : `~astropy.cosmology.Cosmology`

        Returns
        -------
        str
            :mod:`yaml` representation of |Cosmology| object.
        """
        # convert to mapping
        map = obj.to_format("mapping")
        # remove the cosmology class info. It's already recorded in `tag`
        map.pop("cosmology")
        # make the metadata serializable in an order-preserving way.
        map["meta"] = tuple(map["meta"].items())

        return dumper.represent_mapping(tag, map)

    representer.__doc__ = representer.__doc__.format(tag)

    return representer


def yaml_constructor(cls):
    """Cosmology| object from :mod:`yaml` representation.

    Parameters
    ----------
    cls : type
        The class type, e.g. `~astropy.cosmology.LambdaCDM`.

    Returns
    -------
    constructor : callable
        Function to construct |Cosmology| object from :mod:`yaml` representation.
    """

    def constructor(loader, node):
        """Cosmology yaml constructor function.

        Parameters
        ----------
        loader : `~astropy.io.misc.yaml.AstropyLoader`
        node : `yaml.nodes.MappingNode`
            yaml representation of |Cosmology| object.

        Returns
        -------
        `~astropy.cosmology.Cosmology` subclass instance
        """
        # create mapping from YAML node
        map = loader.construct_mapping(node)
        # restore metadata to dict
        map["meta"] = dict(map["meta"])
        # get cosmology class qualified name from node
        cosmology = str(node.tag).split(".")[-1]
        # create Cosmology from mapping
        return from_mapping(map, move_to_meta=False, cosmology=cosmology)

    return constructor


def register_cosmology_yaml(cosmo_cls):
    """Register :mod:`yaml` for Cosmology class.

    Parameters
    ----------
    cosmo_cls : `~astropy.cosmology.Cosmology` class
    """
    fqn = f"{cosmo_cls.__module__}.{cosmo_cls.__qualname__}"
    tag = "!" + QNS.get(
        fqn, fqn
    )  # Possibly sub fully qualified name for a preferred path

    AstropyDumper.add_representer(cosmo_cls, yaml_representer(tag))
    AstropyLoader.add_constructor(tag, yaml_constructor(cosmo_cls))


##############################################################################
# Unified-I/O Functions


def from_yaml(yml, *, cosmology=None):
    """Load `~astropy.cosmology.Cosmology` from :mod:`yaml` object.

    Parameters
    ----------
    yml : str
        :mod:`yaml` representation of |Cosmology| object
    cosmology : str, `~astropy.cosmology.Cosmology` class, or None (optional, keyword-only)
        The expected cosmology class (or string name thereof). This argument is
        is only checked for correctness if not `None`.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Raises
    ------
    TypeError
        If the |Cosmology| object loaded from ``yml`` is not an instance of
        the ``cosmology`` (and ``cosmology`` is not `None`).

    Examples
    --------
    >>> from astropy.cosmology import Cosmology, Planck18
    >>> yml = Planck18.to_format("yaml")
    >>> print(Cosmology.from_format(yml, format="yaml"))
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                  Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
    """
    with u.add_enabled_units(cu):
        cosmo = load(yml)

    # Check argument `cosmology`, if not None
    # This kwarg is required for compatibility with |Cosmology.from_format|
    if isinstance(cosmology, str):
        cosmology = _COSMOLOGY_CLASSES[cosmology]
    if cosmology is not None and not isinstance(cosmo, cosmology):
        raise TypeError(f"cosmology {cosmo} is not an {cosmology} instance.")

    return cosmo


def to_yaml(cosmology, *args):
    r"""Return the cosmology class, parameters, and metadata as a :mod:`yaml` object.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
        The cosmology to serialize.
    *args
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`

    Returns
    -------
    str
        :mod:`yaml` representation of |Cosmology| object

    Examples
    --------
    >>> from astropy.cosmology import Planck18
    >>> Planck18.to_format("yaml")
    "!astropy.cosmology...FlatLambdaCDM\nH0: !astropy.units.Quantity...
    """
    return dump(cosmology)


# ``read`` cannot handle non-path strings.
#  TODO! this says there should be different types of I/O registries.
#        not just hacking object conversion on top of file I/O.
# def yaml_identify(origin, format, *args, **kwargs):
#     """Identify if object uses the yaml format.
#
#     Returns
#     -------
#     bool
#     """
#     itis = False
#     if origin == "read":
#         itis = isinstance(args[1], str) and args[1][0].startswith("!")
#         itis &= format in (None, "yaml")
#
#     return itis


# ===================================================================
# Register

convert_registry.register_reader("yaml", Cosmology, from_yaml)
convert_registry.register_writer("yaml", Cosmology, to_yaml)
# convert_registry.register_identifier("yaml", Cosmology, yaml_identify)
