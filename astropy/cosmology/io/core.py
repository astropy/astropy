# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy

from astropy.table import QTable
from astropy.cosmology import Cosmology

__all__ = ["from_mapping", "to_mapping", "from_table", "to_table"]


def from_mapping(mapping, *, move_to_meta=False):
    """Load `~astropy.cosmology.Cosmology` from mapping object.

    Parameters
    ----------
    mapping : mapping
        Must have field "cosmology".

    move_to_meta : bool (optional, keyword-only)
        Whether to move arguments not in the initialization signature to the
        metadata. This will only have an effect if there is not variable
        keyword-only argument.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    The following is an example JSON serialization of a cosmology.

    .. code-block:: json

        {
          "cosmology": "FlatLambdaCDM",
          "H0": 67.66,
          "Om0": 0.30966,
          "name": "Example",
          "meta": {
            "z_reion": 7.82
          }
        }

    """
    params = copy.deepcopy(mapping)  # so can pop

    cosmology = params.pop("cosmology")
    if isinstance(cosmology, str):
        subclasses = {c.__name__: c for c in Cosmology.__subclasses__(deep=True)}
        cosmology = subclasses[cosmology]

    # select arguments from mapping that are in the cosmo's signature.
    ba = cosmology._init_signature.bind_partial()  # blank set of args
    ba.apply_defaults()  # fill in the defaults
    for k in cosmology._init_signature.parameters.keys():  # iter thru sig
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
        raise TypeError(f"There are unused parameters {params}.")

    return cosmology(*ba.args, **ba.kwargs)


def to_mapping(cosmology):
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

    """
    m = {}
    # start with the cosmology class
    m["cosmology"] = cosmology.__class__
    # get all the immutable inputs
    m.update({k: v for k, v in cosmology._init_arguments.items() if k != "meta"})
    # add the mutable metadata
    m["meta"] = copy.deepcopy(cosmology.meta)

    return m


def from_table(table, index=None, *, move_to_meta=False):
    """Instantiate a `~astropy.cosmology.Cosmology` from a |QTable|.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` class
    table : `~astropy.table.QTable`
    index : int or None, optional
        The row from table.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    The following is an example table serialization of a cosmology.
    .. code-block::

        <QTable length=1>
          name      H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
          str8   float64 float64 float64 float64   float64   float64
        -------- ------- ------- ------- ------- ----------- -------
        Planck18   67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    An example when ``index`` is needed.
    .. code-block::

        <QTable length=3>
          name      H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
          str8   float64 float64 float64 float64   float64   float64
        -------- ------- ------- ------- ------- ----------- --------
        Planck13   67.77 0.30712  2.7255   3.046 0.0 .. 0.06 0.048252
        Planck15   67.74  0.3075  2.7255   3.046 0.0 .. 0.06 0.048600
        Planck18   67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.048970
    """
    if isinstance(index, str):  # convert to row index
        if not table.indices:  # no indices, need to make
            table = table.copy()
            table.add_index('name')
        index = table.loc_indices[index]

    if index is None:
        if len(table) != 1:  # N-D table and no index
            raise ValueError(f"Need to specify a row index for N-D table.")
        else:
            index = 0
    row = table[index]

    # special values
    name = row.columns.get("name", [None])[index]  # get name from column
    meta = copy.deepcopy(row.meta)
    # NOTE: there will be a method for row-specific metadata
    # the cosmology class is in the table's metadata
    cosmology = meta.pop("cosmology", None)

    # turn row into mapping (dict of the arguments)
    mapping = {}
    mapping["cosmology"] = cosmology
    mapping["name"] = name  # do before below to ensure order in dict
    mapping.update({k: v for k, v in zip(row.colnames, row.values())})
    mapping["meta"] = meta

    # build cosmology from map
    return from_mapping(mapping, move_to_meta=move_to_meta)


def to_table(cosmology):
    """Serialize the cosmology into a `~astropy.table.QTable`.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance

    Returns
    -------
    `~astropy.table.QTable`
        With columns for the cosmology parameters, and metadata and
        cosmology class name in the Table's ``meta`` attribute
    """
    p = to_mapping(cosmology)

    # create metadata from mapping
    meta = p.pop("meta")
    meta["cosmology"] = p.pop("cosmology").__name__

    # package parameters into lists for Table parsing
    params = {k: [v] for k, v in p.items()}

    return QTable(params, meta=meta)
