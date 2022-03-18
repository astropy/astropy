# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
The following are private functions, included here **FOR REFERENCE ONLY** since
the io registry cannot be displayed. These functions are registered into
:meth:`~astropy.cosmology.Cosmology.to_format` and
:meth:`~astropy.cosmology.Cosmology.from_format` and should only be accessed
via these methods.
"""  # this is shown in the docs.

import copy
from collections.abc import Mapping

from astropy.cosmology.connect import convert_registry
from astropy.cosmology.core import _COSMOLOGY_CLASSES, Cosmology

__all__ = []  # nothing is publicly scoped


def from_mapping(map, *, move_to_meta=False, cosmology=None):
    """Load `~astropy.cosmology.Cosmology` from mapping object.

    Parameters
    ----------
    map : mapping
        Arguments into the class -- like "name" or "meta".
        If 'cosmology' is None, must have field "cosmology" which can be either
        the string name of the cosmology class (e.g. "FlatLambdaCDM") or the
        class itself.

    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class'
        signature to the Cosmology's metadata. This will only be applied if the
        Cosmology does NOT have a keyword-only argument (e.g. ``**kwargs``).
        Arguments moved to the metadata will be merged with existing metadata,
        preferring specified metadata in the case of a merge conflict
        (e.g. for ``Cosmology(meta={'key':10}, key=42)``, the ``Cosmology.meta``
        will be ``{'key': 10}``).

    cosmology : str, `~astropy.cosmology.Cosmology` class, or None (optional, keyword-only)
        The cosmology class (or string name thereof) to use when constructing
        the cosmology instance. The class also provides default parameter values,
        filling in any non-mandatory arguments missing in 'map'.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a dictionary with
    ``from_mapping``, we will first make a mapping using
    :meth:`~astropy.cosmology.Cosmology.to_format`.

        >>> from astropy.cosmology import Cosmology, Planck18
        >>> cm = Planck18.to_format('mapping')
        >>> cm
        {'cosmology': <class 'astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966,
         'Tcmb0': <Quantity 2.7255 K>, 'Neff': 3.046,
         'm_nu': <Quantity [0. , 0. , 0.06] eV>, 'Ob0': 0.04897,
         'meta': ...

    Now this dict can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = Cosmology.from_format(cm, format="mapping")
        >>> cosmo
        FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                      Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)

    Specific cosmology classes can be used to parse the data. The class'
    default parameter values are used to fill in any information missing in the
    data.

        >>> from astropy.cosmology import FlatLambdaCDM
        >>> del cm["Tcmb0"]  # show FlatLambdaCDM provides default
        >>> FlatLambdaCDM.from_format(cm)
        FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                      Tcmb0=0.0 K, Neff=3.046, m_nu=None, Ob0=0.04897)
    """
    params = dict(map)  # so we are guaranteed to have a poppable map

    # get cosmology
    # 1st from argument. Allows for override of the cosmology, if on file.
    # 2nd from params. This MUST have the cosmology if 'kwargs' did not.
    if cosmology is None:
        cosmology = params.pop("cosmology")
    else:
        params.pop("cosmology", None)  # pop, but don't use
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
    elif params:
        raise TypeError(f"there are unused parameters {params}.")
    # else: pass  # no kwargs, no move-to-meta, and all the params are used

    return cosmology(*ba.args, **ba.kwargs)


def to_mapping(cosmology, *args, cls=dict, cosmology_as_str=False, move_from_meta=False):
    """Return the cosmology class, parameters, and metadata as a `dict`.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    *args
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`
    cls : type (optional, keyword-only)
        `dict` or `collections.Mapping` subclass.
        The mapping type to return. Default is `dict`.
    cosmology_as_str : bool (optional, keyword-only)
        Whether the cosmology value is the class (if `False`, default) or
        the semi-qualified name (if `True`).
    move_from_meta : bool (optional, keyword-only)
        Whether to add the Cosmology's metadata as an item to the mapping (if
        `False`, default) or to merge with the rest of the mapping, preferring
        the original values (if `True`)

    Returns
    -------
    dict
        with key-values for the cosmology parameters and also:
        - 'cosmology' : the class
        - 'meta' : the contents of the cosmology's metadata attribute.
                   If ``move_from_meta`` is `True`, this key is missing and the
                   contained metadata are added to the main `dict`.

    Examples
    --------
    A Cosmology as a mapping will have the cosmology's name and
    parameters as items, and the metadata as a nested dictionary.

        >>> from astropy.cosmology import Planck18
        >>> Planck18.to_format('mapping')
        {'cosmology': <class 'astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966,
         'Tcmb0': <Quantity 2.7255 K>, 'Neff': 3.046,
         'm_nu': <Quantity [0.  , 0.  , 0.06] eV>, 'Ob0': 0.04897,
         'meta': ...

    The dictionary type may be changed with the ``cls`` keyword argument:

        >>> from collections import OrderedDict
        >>> Planck18.to_format('mapping', cls=OrderedDict)
        OrderedDict([('cosmology', <class 'astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM'>),
          ('name', 'Planck18'), ('H0', <Quantity 67.66 km / (Mpc s)>),
          ('Om0', 0.30966), ('Tcmb0', <Quantity 2.7255 K>), ('Neff', 3.046),
          ('m_nu', <Quantity [0.  , 0.  , 0.06] eV>), ('Ob0', 0.04897),
          ('meta', ...

    Sometimes it is more useful to have the name of the cosmology class, not
    the object itself. The keyword argument ``cosmology_as_str`` may be used:

        >>> Planck18.to_format('mapping', cosmology_as_str=True)
        {'cosmology': 'FlatLambdaCDM', ...

    The metadata is normally included as a nested mapping. To move the metadata
    into the main mapping, use the keyword argument ``move_from_meta``. This
    kwarg inverts ``move_to_meta`` in
    ``Cosmology.to_format("mapping", move_to_meta=...)`` where extra items
    are moved to the metadata (if the cosmology constructor does not have a
    variable keyword-only argument -- ``**kwargs``).

        >>> from astropy.cosmology import Planck18
        >>> Planck18.to_format('mapping', move_from_meta=True)
        {'cosmology': <class 'astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM'>,
         'name': 'Planck18', 'Oc0': 0.2607, 'n': 0.9665, 'sigma8': 0.8102, ...
    """
    if not issubclass(cls, (dict, Mapping)):
        raise TypeError(f"'cls' must be a (sub)class of dict or Mapping, not {cls}")

    m = cls()
    # start with the cosmology class & name
    m["cosmology"] = cosmology.__class__.__qualname__ if cosmology_as_str else cosmology.__class__
    m["name"] = cosmology.name  # here only for dict ordering

    meta = copy.deepcopy(cosmology.meta)  # metadata (mutable)
    if move_from_meta:
        # Merge the mutable metadata. Since params are added later they will
        # be preferred in cases of overlapping keys. Likewise, need to pop
        # cosmology and name from meta.
        meta.pop("cosmology", None)
        meta.pop("name", None)
        m.update(meta)

    # Add all the immutable inputs
    m.update({k: v for k, v in cosmology._init_arguments.items()
              if k not in ("meta", "name")})
    # Lastly, add the metadata, if haven't already (above)
    if not move_from_meta:
        m["meta"] = meta  # TODO? should meta be type(cls)

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

convert_registry.register_reader("mapping", Cosmology, from_mapping)
convert_registry.register_writer("mapping", Cosmology, to_mapping)
convert_registry.register_identifier("mapping", Cosmology, mapping_identify)
