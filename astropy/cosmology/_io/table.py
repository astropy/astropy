# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy

import numpy as np

from astropy.io import registry as io_registry
from astropy.table import Table, QTable
from astropy.cosmology import Cosmology, _COSMOLOGY_REGISTRY

from .mapping import from_mapping, to_mapping


def from_table(table, index=None, *, move_to_meta=False, **kwargs):
    """Instantiate a `~astropy.cosmology.Cosmology` from a |QTable|.

    Parameters
    ----------
    table : `~astropy.table.QTable`
    index : int, str, or None, optional
        Needed to select the row in tables with multiple rows. ``index`` can be
        an integer for the row number or, if the table is indexed by a column,
        the value of that column. If the table is not indexed and ``index``
        is a string, the "name" column is used as the indexing column.
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
        If 'format' is a kwarg, it must be 'table'.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a Table with
    ``from_table``, we will first make a |QTable| using
    :func:`~astropy.cosmology._io.to_table`.

        >>> from astropy.cosmology import Planck18
        >>> from astropy.cosmology._io.table import to_table
        >>> ct = to_table(Planck18); ct
        <QTable length=1>
          name      H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
          str8   float64 float64 float64 float64   float64   float64
        -------- ------- ------- ------- ------- ----------- -------
        Planck18   67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    Now this table can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> from astropy.cosmology._io.table import from_table
        >>> cosmo = from_table(ct)
        >>> cosmo
        FlatLambdaCDM(name="Planck18", H0=67.7 km / (Mpc s), Om0=0.31,
                      Tcmb0=2.725 K, Neff=3.05, m_nu=[0. 0. 0.06] eV, Ob0=0.049)

    For tables with multiple rows of cosmological parameters, the ``index``
    argument is needed to select the correct row. The index can be an integer
    for the row number or, if the table is indexed by a column, the value of
    that column. If the table is not indexed and ``index`` is a string, the
    "name" column is used as the indexing column.

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

    For further examples, see :doc:`astropy:cosmology/io`
    """
    # check 'format' correctness
    format = kwargs.pop("format", "table")
    if format != "table":  # check that if specified, it's a table
        raise ValueError(f"'format', if specified, must be 'table' not {format}")

    # ------------------
    # Get row from table

    # string index uses the indexed column on the table to find the row index.
    if isinstance(index, str):
        if not table.indices:  # no indexing column, find by string match
            index = np.where(table['name'] == index)[0][0]
            # TODO! error if no match
        else:
            index = table.loc_indices[index]  # need to convert to row index (int)

    # no index is needed for a 1-row table. For a multi-row table...
    if index is None:
        if len(table) != 1:  # multi-row table and no index
            raise ValueError("need to select a specific row (e.g. index=1) when "
                             "constructing a Cosmology from a multi-row table.")
        else:
            index = 0
    row = table[index]  # index is now the row index (int)

    # ------------------
    # parse row to cosmo

    # special values
    name = row['name'] if 'name' in row.columns else None  # get name from column
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


def to_table(cosmology, *args, **kwargs):
    """Serialize the cosmology into a `~astropy.table.QTable`.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    *args, **kwargs
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`

    Returns
    -------
    `~astropy.table.QTable`
        With columns for the cosmology parameters, and metadata and
        cosmology class name in the Table's ``meta`` attribute

    Examples
    --------
    A Cosmology as a `~astropy.table.QTable` will have the cosmology's name and
    parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> from astropy.cosmology._io.table import to_table
        >>> ct = to_table(Planck18); ct
        <QTable length=1>
          name      H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
          str8   float64 float64 float64 float64   float64   float64
        -------- ------- ------- ------- ------- ----------- -------
        Planck18   67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    The cosmological class and other metadata, e.g. a paper reference, are in
    the Table's metadata.

        >>> ct.meta
        OrderedDict([..., ('cosmology', 'FlatLambdaCDM')])
    """
    # start by getting a map representation. This requires minimal repackaging.
    p = to_mapping(cosmology)

    # create metadata from mapping
    meta = p.pop("meta")
    meta["cosmology"] = p.pop("cosmology").__name__  # move class to Table meta

    # package parameters into lists for Table parsing
    params = {k: [v] for k, v in p.items()}

    return QTable(params, meta=meta)


def table_identify(origin, format, *args, **kwargs):
    """Identify if object uses the Table format.

    Returns
    -------
    bool
    """
    if origin == "write":
        itis = (format == "table")
    elif origin == "read":
        itis = isinstance(args[1], Table) and (format in (None, "table"))

    return itis


# ===================================================================
# Register

io_registry.register_reader("table", Cosmology, from_table)
io_registry.register_writer("table", Cosmology, to_table)
io_registry.register_identifier("table", Cosmology, table_identify)
