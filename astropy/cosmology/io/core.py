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
        Must have field "cosmology" which can be either the string name of the
        cosmology class (e.g. "FlatLambdaCDM") or the class itself.
        The remaining fields are arguments into the class --
        like "name" or "meta".

    move_to_meta : bool (optional, keyword-only)
        Whether to move arguments not in the initialization signature to the
        metadata. This will only have an effect if there is not a variable
        keyword-only argument (e.g. ``**kwargs``). Metadata set in the field
        "meta" has priority and will not be overwritten.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a dictionary with
    ``from_mapping``, we will first make a mapping using
    :func:`~astropy.cosmology.io.to_mapping`.

        >>> from astropy.cosmology import io, Planck18
        >>> cm = io.to_mapping(Planck18); cm
        {'cosmology': <class 'astropy.cosmology.core.FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': 67.66, 'Om0': 0.30966, 'Tcmb0': 2.7255,
         'Neff': 3.046, 'm_nu': [0.0, 0.0, 0.06], 'Ob0': 0.04897,
         'meta': ...

    Now this dict can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = io.from_mapping(cm)
        >>> cosmo
        FlatLambdaCDM(name="Planck18", H0=67.7 km / (Mpc s), Om0=0.31,
                      Tcmb0=2.725 K, Neff=3.05, m_nu=[0. 0. 0.06] eV, Ob0=0.049)

    See Also
    --------
    astropy.cosmology.io.to_mapping : Represent as dictionary.
    astropy.cosmology.io.from_table : Load Cosmology from |QTable|.
    astropy.cosmology.io.to_table : Represent as |QTable|.
    astropy.cosmology.Cosmology.read : Has a ``from_mapping`` convenience method.
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
    # else: pass  # no kwargs, no move-to-meta, and all the params are used

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

    Examples
    --------
    A Cosmology as a mapping will have the cosmology's name and
    parameters as items, and the metadata as nested dictionary.

        >>> from astropy.cosmology import io, Planck18
        >>> io.to_mapping(Planck18)
        {'cosmology': <class 'astropy.cosmology.core.FlatLambdaCDM'>,
         'name': 'Planck18', 'H0': 67.66, 'Om0': 0.30966, 'Tcmb0': 2.7255,
         'Neff': 3.046, 'm_nu': [0.0, 0.0, 0.06], 'Ob0': 0.04897,
         'meta': ...

    See Also
    --------
    astropy.cosmology.io.from_mapping : Load Cosmology from dictionary.
    astropy.cosmology.io.from_table : Load Cosmology from |QTable|.
    astropy.cosmology.io.to_table : Represent as |QTable|.
    astropy.cosmology.Cosmology.write : Has a ``to_mapping`` convenience method.
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


def from_table(table, index=None, *, move_to_meta=False):
    """Instantiate a `~astropy.cosmology.Cosmology` from a |QTable|.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` class
    table : `~astropy.table.QTable`
    index : int or None, optional
        Needed to select the row in tables with multiple rows. ``index`` can be
        an integer for the row number or, if the table is indexed by a column,
        the value of that column. If the table is not indexed and ``indexed``
        is a string, the "name" column is used as the indexing column.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a Table with
    ``from_table``, we will first make a |QTable| using
    :func:`~astropy.cosmology.io.to_table`.

        >>> from astropy.cosmology import io, Planck18
        >>> ct = io.to_table(Planck18); ct
        <QTable length=1>
          name      H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
          str8   float64 float64 float64 float64   float64   float64
        -------- ------- ------- ------- ------- ----------- -------
        Planck18   67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    Now this table can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = io.from_table(ct)
        >>> cosmo
        FlatLambdaCDM(name="Planck18", H0=67.7 km / (Mpc s), Om0=0.31,
                      Tcmb0=2.725 K, Neff=3.05, m_nu=[0. 0. 0.06] eV, Ob0=0.049)

    For tables with multiple rows of cosmological parameters, the ``index``
    argument is needed to select the correct row. The index can be an integer
    for the row number or, if the table is indexed by a column, the value of
    that column. If the table is not indexed and ``indexed``
    is a string, the "name" column is used as the indexing column.

    Here is an example where ``index`` is needed and can be either an integer
    (for the row number) or the name of one of the cosmologies, e.g. 'Planck15'.

    .. code-block::

        <QTable length=3>
          name      H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
          str8   float64 float64 float64 float64   float64   float64
        -------- ------- ------- ------- ------- ----------- --------
        Planck13   67.77 0.30712  2.7255   3.046 0.0 .. 0.06 0.048252
        Planck15   67.74  0.3075  2.7255   3.046 0.0 .. 0.06 0.048600
        Planck18   67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.048970

    See Also
    --------
    astropy.cosmology.io.to_table : Represent as |QTable|.
    astropy.cosmology.io.from_mapping : Load Cosmology from dictionary.
    astropy.cosmology.io.to_mapping : Represent as dictionary.
    astropy.cosmology.Cosmology.read : Has a ``from_table`` convenience method.
    """
    # string index uses the indexed column on the table to find the row index.
    if isinstance(index, str):
        if not table.indices:  # no indexing column, need to make
            table = table.copy()  # unfortunately, need to copy.
            table.add_index('name')
        index = table.loc_indices[index]  # need to convert to row index (int)

    # no index is needed for a 1D table. For an N-D table...
    if index is None:
        if len(table) != 1:  # N-D table and no index
            raise ValueError(f"Need to specify a row index for N-D table.")
        else:
            index = 0
    row = table[index]  # index is now the row index (int)

    # special values
    name = row.columns.get("name", [None])[index]  # get name from column
    meta = copy.deepcopy(row.meta)
    # NOTE: there will be a method for row-specific metadata
    # the cosmology class must be in the table's metadata
    cosmology = meta.pop("cosmology")

    # turn row into mapping (dict of the arguments)
    mapping = {}
    mapping["cosmology"] = cosmology
    mapping["name"] = name  # do before below to ensure order in dict
    mapping.update({k: v for k, v in zip(row.colnames, row.values())
                    if k != "name"})
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

    Examples
    --------
    A Cosmology as a `~astropy.table.QTable` will have the cosmology's name and
    parameters as columns.

        >>> from astropy.cosmology import io, Planck18
        >>> ct = io.to_table(Planck18); ct
        <QTable length=1>
          name      H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
          str8   float64 float64 float64 float64   float64   float64
        -------- ------- ------- ------- ------- ----------- -------
        Planck18   67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    The cosmological class and other metadata, e.g. a paper reference, are in
    the Table's metadata.

        >>> ct.meta
        OrderedDict([..., ('cosmology', 'FlatLambdaCDM')])

    See Also
    --------
    astropy.cosmology.io.from_table : Load Cosmology from |QTable|.
    astropy.cosmology.io.from_mapping : Load Cosmology from dictionary.
    astropy.cosmology.io.to_mapping : Represent as dictionary.
    astropy.cosmology.Cosmology.write : Has a ``to_table`` convenience method.
    """
    # start by getting a map representation. This requires minimal repackaging.
    p = to_mapping(cosmology)

    # create metadata from mapping
    meta = p.pop("meta")
    meta["cosmology"] = p.pop("cosmology").__name__

    # package parameters into lists for Table parsing
    params = {k: [v] for k, v in p.items()}

    return QTable(params, meta=meta)
