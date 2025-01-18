# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""|Cosmology| <-> Mapping I/O, using |Cosmology.to_format| and |Cosmology.from_format|.

This module provides functions to transform a |Cosmology| instance to a mapping
(`dict`-like) object and vice versa, from a mapping object back to a |Cosmology|
instance. The functions are registered with ``convert_registry`` under the format name
"mapping". The mapping object is a `dict`-like object, with the cosmology's parameters
and metadata as items. `dict` is a fundamental data structure in Python, and this
representation of a |Cosmology| is useful for translating between many serialization and
storage formats, or even passing arguments to functions.

We start with the simple case of outputting a |Cosmology| as a mapping.

    >>> from astropy.cosmology import Cosmology, Planck18
    >>> cm = Planck18.to_format('mapping')
    >>> cm
    {'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>,
     'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966,
     'Tcmb0': <Quantity 2.7255 K>, 'Neff': 3.046,
     'm_nu': <Quantity [0. , 0. , 0.06] eV>, 'Ob0': 0.04897,
     'meta': ...

``cm`` is a `dict`, with the cosmology's parameters and metadata as items.

How might we use this `dict`? One use is to unpack the `dict` into a function:

    >>> def function(H0, Tcmb0, **kwargs): ...
    >>> function(**cm)

Another use is to merge the `dict` with another `dict`:

    >>> cm2 = {'H0': 70, 'Tcmb0': 2.7}
    >>> cm | cm2
    {'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>, ..., 'H0': 70, ...}

Most saliently, the `dict` can also be used to construct a new cosmological instance
identical to the |Planck18| cosmology from which it was generated.

    >>> cosmo = Cosmology.from_format(cm, format="mapping")
    >>> cosmo
    FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

How did |Cosmology.from_format| know to return an instance of the |FlatLambdaCDM| class?
The mapping object has a field ``cosmology`` which can be either the string name of the
cosmology class (e.g. "FlatLambdaCDM") or the class itself.

This field can be omitted under two conditions.

1. If the cosmology class is passed as the ``cosmology`` keyword argument to
   |Cosmology.from_format|,
2. If a specific cosmology class, e.g. |FlatLambdaCDM|, is used to parse the data.

To the first point, we can pass the cosmology class as the ``cosmology`` keyword
argument to |Cosmology.from_format|.

    >>> del cm["cosmology"]  # remove cosmology class

    >>> Cosmology.from_format(cm, cosmology="FlatLambdaCDM")
    FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

To the second point, we can use specific cosmology class to parse the data.

    >>> from astropy.cosmology import FlatLambdaCDM
    >>> FlatLambdaCDM.from_format(cm)
    FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

Also, the class' default parameter values are used to fill in any information missing in
the data. For example, if ``Tcmb0`` is missing, the default value of 0.0 K is used.

    >>> del cm["Tcmb0"]  # show FlatLambdaCDM provides default
    >>> FlatLambdaCDM.from_format(cm)
    FlatLambdaCDM(name='Planck18', H0=..., Tcmb0=<Quantity 0. K>, ...)

If instead of *missing* information, there is *extra* information, there are a few
options. The first is to use the ``move_to_meta`` keyword argument to move fields that
are not in the Cosmology constructor to the Cosmology's metadata.

    >>> cm2 = cm | {"extra": 42, "cosmology": "FlatLambdaCDM"}
    >>> cosmo = Cosmology.from_format(cm2, move_to_meta=True)
    >>> cosmo.meta
    {'extra': 42, ...}

Alternatively, the ``rename`` keyword argument can be used to rename keys in the mapping
to fields of the |Cosmology|. This is crucial when the mapping has keys that are not
valid arguments to the |Cosmology| constructor.

    >>> cm3 = dict(cm)  # copy
    >>> cm3["cosmo_cls"] = "FlatLambdaCDM"
    >>> cm3["cosmo_name"] = cm3.pop("name")

    >>> rename = {'cosmo_cls': 'cosmology', 'cosmo_name': 'name'}
    >>> Cosmology.from_format(cm3, rename=rename)
    FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 0. K>, Neff=3.046, m_nu=None, Ob0=0.04897)

Let's take a closer look at |Cosmology.to_format|, because there a lot of options, to
tailor the output to specific needs.

The dictionary type may be changed with the ``cls`` keyword argument:

    >>> from collections import OrderedDict
    >>> Planck18.to_format('mapping', cls=OrderedDict)
    OrderedDict({'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>, 'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966, 'Tcmb0': <Quantity 2.7255 K>, 'Neff': 3.046, 'm_nu': <Quantity [0.  , 0.  , 0.06] eV>, 'Ob0': 0.04897, 'meta': {...}})

Sometimes it is more useful to have the name of the cosmology class, not the type
itself. The keyword argument ``cosmology_as_str`` may be used:

    >>> Planck18.to_format('mapping', cosmology_as_str=True)
    {'cosmology': 'FlatLambdaCDM', ...

The metadata is normally included as a nested mapping. To move the metadata into the
main mapping, use the keyword argument ``move_from_meta``. This kwarg inverts
``move_to_meta`` in ``Cosmology.to_format("mapping", move_to_meta=...)`` where extra
items are moved to the metadata (if the cosmology constructor does not have a variable
keyword-only argument -- ``**kwargs``).

    >>> from astropy.cosmology import Planck18
    >>> Planck18.to_format('mapping', move_from_meta=True)
    {'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>,
        'name': 'Planck18', 'Oc0': 0.2607, 'n': 0.9665, 'sigma8': 0.8102, ...

Lastly, the keys in the mapping may be renamed with the ``rename`` keyword.

    >>> rename = {'cosmology': 'cosmo_cls', 'name': 'cosmo_name'}
    >>> Planck18.to_format('mapping', rename=rename)
    {'cosmo_cls': <class 'astropy.cosmology...FlatLambdaCDM'>,
     'cosmo_name': 'Planck18', ...
"""

from __future__ import annotations

__all__: list[str] = []  # nothing is publicly scoped

import copy
import inspect
from collections.abc import Mapping, MutableMapping
from typing import TYPE_CHECKING, Any, TypeVar

from astropy.cosmology._src.core import _COSMOLOGY_CLASSES, Cosmology
from astropy.cosmology._src.io.connect import convert_registry

if TYPE_CHECKING:
    from astropy.cosmology._src.typing import _CosmoT

    _MapT = TypeVar("_MapT", MutableMapping[str, Any])


def _rename_map(
    map: Mapping[str, Any], /, renames: Mapping[str, str]
) -> dict[str, Any]:
    """Apply rename to map."""
    if common_names := set(renames.values()).intersection(map):
        raise ValueError(
            "'renames' values must be disjoint from 'map' keys, "
            f"the common keys are: {common_names}"
        )
    return {renames.get(k, k): v for k, v in map.items()}  # dict separate from input


def _get_cosmology_class(
    cosmology: type[_CosmoT] | str | None, params: dict[str, Any], /
) -> type[_CosmoT]:
    # get cosmology
    # 1st from argument. Allows for override of the cosmology, if on file.
    # 2nd from params. This MUST have the cosmology if 'kwargs' did not.
    if cosmology is None:
        cosmology = params.pop("cosmology")
    else:
        params.pop("cosmology", None)  # pop, but don't use
    # if string, parse to class
    return _COSMOLOGY_CLASSES[cosmology] if isinstance(cosmology, str) else cosmology


def from_mapping(
    mapping: Mapping[str, Any],
    /,
    *,
    move_to_meta: bool = False,
    cosmology: str | type[_CosmoT] | None = None,
    rename: Mapping[str, str] | None = None,
) -> _CosmoT:
    """Load `~astropy.cosmology.Cosmology` from mapping object.

    Parameters
    ----------
    mapping : Mapping
        Arguments into the class -- like "name" or "meta". If 'cosmology' is None, must
        have field "cosmology" which can be either the string name of the cosmology
        class (e.g. "FlatLambdaCDM") or the class itself.

    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class' signature
        to the Cosmology's metadata. This will only be applied if the Cosmology does NOT
        have a keyword-only argument (e.g. ``**kwargs``). Arguments moved to the
        metadata will be merged with existing metadata, preferring specified metadata in
        the case of a merge conflict (e.g. for ``Cosmology(meta={'key':10}, key=42)``,
        the ``Cosmology.meta`` will be ``{'key': 10}``).

    cosmology : str, |Cosmology| class, or None (optional, keyword-only)
        The cosmology class (or string name thereof) to use when constructing the
        cosmology instance. The class also provides default parameter values, filling in
        any non-mandatory arguments missing in 'map'.

    rename : Mapping[str, str] or None (optional, keyword-only)
        A mapping of keys in ``map`` to fields of the `~astropy.cosmology.Cosmology`.

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
        {'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966,
         'Tcmb0': <Quantity 2.7255 K>, 'Neff': 3.046,
         'm_nu': <Quantity [0. , 0. , 0.06] eV>, 'Ob0': 0.04897,
         'meta': ...

    Now this dict can be used to load a new cosmological instance identical to the
    |Planck18| cosmology from which it was generated.

        >>> cosmo = Cosmology.from_format(cm, format="mapping")
        >>> cosmo
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

    The ``cosmology`` field can be omitted if the cosmology class (or its string name)
    is passed as the ``cosmology`` keyword argument to |Cosmology.from_format|.

        >>> del cm["cosmology"]  # remove cosmology class
        >>> Cosmology.from_format(cm, cosmology="FlatLambdaCDM")
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

    Alternatively, specific cosmology classes can be used to parse the data.

        >>> from astropy.cosmology import FlatLambdaCDM
        >>> FlatLambdaCDM.from_format(cm)
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

    When using a specific cosmology class, the class' default parameter values are used
    to fill in any missing information.

        >>> del cm["Tcmb0"]  # show FlatLambdaCDM provides default
        >>> FlatLambdaCDM.from_format(cm)
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 0. K>, Neff=3.046, m_nu=None, Ob0=0.04897)

    The ``move_to_meta`` keyword argument can be used to move fields that are not in the
    Cosmology constructor to the Cosmology's metadata. This is useful when the
    dictionary contains extra information that is not part of the Cosmology.

        >>> cm2 = cm | {"extra": 42, "cosmology": "FlatLambdaCDM"}
        >>> cosmo = Cosmology.from_format(cm2, move_to_meta=True)
        >>> cosmo.meta
        {'extra': 42, ...}

    The ``rename`` keyword argument can be used to rename keys in the mapping to fields
    of the |Cosmology|. This is crucial when the mapping has keys that are not valid
    arguments to the |Cosmology| constructor.

        >>> cm3 = dict(cm)  # copy
        >>> cm3["cosmo_cls"] = "FlatLambdaCDM"
        >>> cm3["cosmo_name"] = cm3.pop("name")

        >>> rename = {'cosmo_cls': 'cosmology', 'cosmo_name': 'name'}
        >>> Cosmology.from_format(cm3, rename=rename)
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 0. K>, Neff=3.046, m_nu=None, Ob0=0.04897)
    """
    # Rename keys, if given a ``renames`` dict.
    # Also, make a copy of the mapping, so we can pop from it.
    params = _rename_map(dict(mapping), renames=rename or {})

    # Get cosmology class
    cosmology = _get_cosmology_class(cosmology, params)

    # select arguments from mapping that are in the cosmo's signature.
    sig = inspect.signature(cosmology)
    ba = sig.bind_partial()  # blank set of args
    ba.apply_defaults()  # fill in the defaults
    for k in sig.parameters.keys():
        if k in params:  # transfer argument, if in params
            ba.arguments[k] = params.pop(k)

    # deal with remaining params. If there is a **kwargs use that, else
    # allow to transfer to metadata. Raise TypeError if can't.
    lastp = next(reversed(sig.parameters.values()))
    if lastp.kind == 4:  # variable keyword-only
        ba.arguments[lastp.name] = params
    elif move_to_meta:  # prefers current meta, which was explicitly set
        meta = ba.arguments["meta"] or {}  # (None -> dict)
        ba.arguments["meta"] = {**params, **meta}
    elif params:
        raise TypeError(f"there are unused parameters {params}.")
    # else: pass  # no kwargs, no move-to-meta, and all the params are used

    return cosmology(*ba.args, **ba.kwargs)


def to_mapping(
    cosmology: Cosmology,
    *args: object,
    cls: type[_MapT] = dict,
    cosmology_as_str: bool = False,
    move_from_meta: bool = False,
    rename: Mapping[str, str] | None = None,
) -> _MapT:
    """Return the cosmology class, parameters, and metadata as a `dict`.

    Parameters
    ----------
    cosmology : :class:`~astropy.cosmology.Cosmology`
        The cosmology instance to convert to a mapping.
    *args : object
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
    rename : Mapping[str, str] or None (optional, keyword-only)
        A mapping of field names of the :class:`~astropy.cosmology.Cosmology` to keys in
        the map.

    Returns
    -------
    MutableMapping[str, Any]
        A mapping of type ``cls``, by default a `dict`.
        Has key-values for the cosmology parameters and also:
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
        {'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966,
         'Tcmb0': <Quantity 2.7255 K>, 'Neff': 3.046,
         'm_nu': <Quantity [0.  , 0.  , 0.06] eV>, 'Ob0': 0.04897,
         'meta': ...

    The dictionary type may be changed with the ``cls`` keyword argument:

        >>> from collections import OrderedDict
        >>> Planck18.to_format('mapping', cls=OrderedDict)
        OrderedDict({'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>, 'name': 'Planck18', 'H0': <Quantity 67.66 km / (Mpc s)>, 'Om0': 0.30966, 'Tcmb0': <Quantity 2.7255 K>, 'Neff': 3.046, 'm_nu': <Quantity [0.  , 0.  , 0.06] eV>, 'Ob0': 0.04897, 'meta': {...}})

    Sometimes it is more useful to have the name of the cosmology class, not
    the type itself. The keyword argument ``cosmology_as_str`` may be used:

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
        {'cosmology': <class 'astropy.cosmology...FlatLambdaCDM'>,
         'name': 'Planck18', 'Oc0': 0.2607, 'n': 0.9665, 'sigma8': 0.8102, ...

    Lastly, the keys in the mapping may be renamed with the ``rename`` keyword.

        >>> rename = {'cosmology': 'cosmo_cls', 'name': 'cosmo_name'}
        >>> Planck18.to_format('mapping', rename=rename)
        {'cosmo_cls': <class 'astropy.cosmology...FlatLambdaCDM'>,
         'cosmo_name': 'Planck18', ...
    """
    if not issubclass(cls, (dict, Mapping)):
        raise TypeError(f"'cls' must be a (sub)class of dict or Mapping, not {cls}")

    m = cls()
    # start with the cosmology class & name
    m["cosmology"] = (
        cosmology.__class__.__qualname__ if cosmology_as_str else cosmology.__class__
    )
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
    m.update(cosmology.parameters)
    # Lastly, add the metadata, if haven't already (above)
    if not move_from_meta:
        m["meta"] = meta  # TODO? should meta be type(cls)
    # Rename keys
    return m if rename is None else _rename_map(m, rename)


def mapping_identify(
    origin: str, format: str | None, *args: object, **kwargs: object
) -> bool:
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
