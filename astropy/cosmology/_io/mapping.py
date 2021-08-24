# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
from collections.abc import Mapping

import numpy as np

from astropy.cosmology.core import _COSMOLOGY_CLASSES, Cosmology
from astropy.io import registry as io_registry
from astropy.table import QTable


def from_mapping(map, *, move_to_meta=False, **kwargs):
    """Load `~astropy.cosmology.Cosmology` from mapping object.

    Parameters
    ----------
    map : mapping
        Must have field "cosmology" which can be either the string name of the
        cosmology class (e.g. "FlatLambdaCDM") or the class itself.
        The remaining fields are arguments into the class --
        like "name" or "meta".

    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class'
        signature to the Cosmology's metadata. This will only be applied if the
        Cosmology does NOT have a keyword-only argument (e.g. ``**kwargs``).
        Arguments moved to the metadata will be merged with existing metadata,
        preferring specified metadata in the case of a merge conflict
        (e.g. for ``Cosmology(meta={'key':10}, key=42)``, the ``Cosmology.meta``
        will be ``{'key': 10}``).
    **kwargs
        Not used.
        If 'format' is a kwarg, it must be 'mapping'.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a dictionary with
    ``from_mapping``, we will first make a mapping using
    :meth:`~astropy.cosmology.Cosmology.to_format`.

        >>> from astropy.cosmology import Cosmology, Planck18
        >>> cm = Planck18.to_format('mapping'); cm
        {'cosmology': <class 'astropy.cosmology.core.FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966,
         'Tcmb0': 2.7255, 'Neff': 3.046, 'm_nu': <Quantity [0. , 0. , 0.06] eV>,
         'Ob0': 0.04897, 'meta': ...

    Now this dict can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = Cosmology.from_format(cm, format="mapping")
        >>> cosmo
        FlatLambdaCDM(name="Planck18", H0=67.7 km / (Mpc s), Om0=0.31,
                      Tcmb0=2.725 K, Neff=3.05, m_nu=[0. 0. 0.06] eV, Ob0=0.049)
    """
    params = copy.deepcopy(map)  # so can pop

    # get cosmology
    # 1st from 'kwargs'. Allows for override of the cosmology, if on file.
    # 2nd from params. This MUST have the cosmology if 'kwargs' did not.
    if "cosmology" in kwargs:
        cosmology = kwargs.pop("cosmology")
        params.pop("cosmology", None)  # pop, but don't use
    else:
        cosmology = params.pop("cosmology")
    # if string, parse to class
    if isinstance(cosmology, str):
        cosmology = _COSMOLOGY_CLASSES[cosmology]

    # select arguments from mapping that are in the cosmo's signature.
    ba = cosmology._init_signature.bind_partial()  # blank set of args
    ba.apply_defaults()  # fill in the defaults
    for k in cosmology._init_signature.parameters.keys():
        if k in params:  # transfer argument, if in params
            ba.arguments[k] = params.pop(k)

    # deal with remaining params. If there is a **kwargs use that, else
    # allow to transfer to metadata. Raise TypeError if can't.
    lastp = tuple(cosmology._init_signature.parameters.values())[-1]
    if lastp.kind == 4:  # variable keyword-only
        ba.arguments[lastp.name] = params
    elif move_to_meta:  # prefers current meta, which was explicitly set
        meta = ba.arguments["meta"] or {}  # (None -> dict)
        ba.arguments["meta"] = {**params, **meta}
    elif bool(params):
        raise TypeError(f"there are unused parameters {params}.")
    # else: pass  # no kwargs, no move-to-meta, and all the params are used

    return cosmology(*ba.args, **ba.kwargs)


def to_mapping(cosmology, *args, **kwargs):
    """Return the cosmology class, inputs, and metadata as a dict.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance

    Returns
    -------
    dict
        with key-values for the cosmology parameters and also:

        - 'cosmology' : the class
        - 'meta' : the contents of the cosmology's metadata attribute
    *args, **kwargs
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`

    Examples
    --------
    A Cosmology as a mapping will have the cosmology's name and
    parameters as items, and the metadata as a nested dictionary.

        >>> from astropy.cosmology import Planck18
        >>> Planck18.to_format('mapping')
        {'cosmology': <class 'astropy.cosmology.core.FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966,
         'Tcmb0': 2.7255, 'Neff': 3.046, 'm_nu': <Quantity [0.  , 0.  , 0.06] eV>,
         'Ob0': 0.04897, 'meta': ...

    """

    m = {}
    # start with the cosmology class & name
    m["cosmology"] = cosmology.__class__
    m["name"] = cosmology.name  # here only for dict ordering
    # get all the immutable inputs
    m.update({k: v for k, v in cosmology._init_arguments.items()
              if k not in ("meta", "name")})
    # add the mutable metadata
    m["meta"] = copy.deepcopy(cosmology.meta)

    return m


def mapping_identify(origin, format, *args, **kwargs):
    """Identify if object uses the mapping format.

    Returns
    -------
    bool
    """
    itis = False
    if origin == "read":
        itis = isinstance(args[1], Mapping) and (format in (None, "mapping"))

    return itis


# ===================================================================
# Register

io_registry.register_reader("mapping", Cosmology, from_mapping)
io_registry.register_writer("mapping", Cosmology, to_mapping)
io_registry.register_identifier("mapping", Cosmology, mapping_identify)
