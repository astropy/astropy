# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""|Cosmology| <-> |Row| I/O, using |Cosmology.to_format| and |Cosmology.from_format|.

A `~astropy.cosmology.Cosmology` as a `~astropy.table.Row` will have
the cosmology's name and parameters as columns.

    >>> from astropy.cosmology import Planck18
    >>> cr = Planck18.to_format("astropy.row")
    >>> cr
    <Row index=0>
        cosmology     name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                            km / (Mpc s)            K                 eV
        str13       str8     float64    float64 float64 float64  float64[3] float64
    ------------- -------- ------------ ------- ------- ------- ----------- -------
    FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

The cosmological class and other metadata, e.g. a paper reference, are in
the Table's metadata.

    >>> cr.meta
    {'Oc0': 0.2607, 'n': 0.9665, ...}

Now this row can be used to load a new cosmological instance identical
to the ``Planck18`` cosmology from which it was generated.

    >>> cosmo = Cosmology.from_format(cr, format="astropy.row")
    >>> cosmo
    FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

For more information on the argument options, see :ref:`cosmology_io_builtin-table`.
"""

from __future__ import annotations

import copy
from collections import defaultdict
from typing import TYPE_CHECKING

from astropy.table import QTable, Row

# isort: split
from astropy.cosmology._src.core import Cosmology
from astropy.cosmology._src.io.connect import convert_registry

from .mapping import from_mapping

if TYPE_CHECKING:
    from collections.abc import Mapping

    from astropy.cosmology._src.typing import _CosmoT
    from astropy.table import Table


def from_row(
    row: Row,
    *,
    move_to_meta: bool = False,
    cosmology: str | type[_CosmoT] | None = None,
    rename: Mapping[str, str] | None = None,
) -> _CosmoT:
    """Instantiate a `~astropy.cosmology.Cosmology` from a `~astropy.table.Row`.

    Parameters
    ----------
    row : `~astropy.table.Row`
        The object containing the Cosmology information.
    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class'
        signature to the Cosmology's metadata. This will only be applied if the
        Cosmology does NOT have a keyword-only argument (e.g. ``**kwargs``).
        Arguments moved to the metadata will be merged with existing metadata,
        preferring specified metadata in the case of a merge conflict
        (e.g. for ``Cosmology(meta={'key':10}, key=42)``, the ``Cosmology.meta``
        will be ``{'key': 10}``).

    cosmology : str, type, or None (optional, keyword-only)
        The cosmology class (or string name thereof) to use when constructing
        the cosmology instance. The class also provides default parameter values,
        filling in any non-mandatory arguments missing in 'table'.

    rename : Mapping[str, str] or None (optional, keyword-only)
        A mapping of column names in the row to field names of the |Cosmology|.

    Returns
    -------
    `~astropy.cosmology.Cosmology`

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a Row with
    ``from_row``, we will first make a `~astropy.table.Row` using
    :func:`~astropy.cosmology.Cosmology.to_format`.

        >>> from astropy.cosmology import Cosmology, Planck18
        >>> cr = Planck18.to_format("astropy.row")
        >>> cr
        <Row index=0>
          cosmology     name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                               km / (Mpc s)            K                 eV
            str13       str8     float64    float64 float64 float64  float64[3] float64
        ------------- -------- ------------ ------- ------- ------- ----------- -------
        FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    Now this row can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = Cosmology.from_format(cr, format="astropy.row")
        >>> cosmo
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

    The ``cosmology`` information (column or metadata) may be omitted if the cosmology
    class (or its string name) is passed as the ``cosmology`` keyword argument to
    |Cosmology.from_format|.

        >>> del cr.columns["cosmology"]  # remove cosmology from metadata
        >>> Cosmology.from_format(cr, cosmology="FlatLambdaCDM")
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

    Alternatively, specific cosmology classes can be used to parse the data.

        >>> from astropy.cosmology import FlatLambdaCDM
        >>> FlatLambdaCDM.from_format(cr)
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 2.7255 K>, Neff=3.046, m_nu=<Quantity [0.  , 0.  , 0.06] eV>, Ob0=0.04897)

    When using a specific cosmology class, the class' default parameter values are used
    to fill in any missing information.

        >>> del cr.columns["Tcmb0"]  # show FlatLambdaCDM provides default
        >>> FlatLambdaCDM.from_format(cr)
        FlatLambdaCDM(name='Planck18', H0=<Quantity 67.66 km / (Mpc s)>, Om0=0.30966, Tcmb0=<Quantity 0. K>, Neff=3.046, m_nu=None, Ob0=0.04897)

    If a `~astropy.table.Row` object has columns that do not match the fields of the
    `~astropy.cosmology.Cosmology` class, they can be mapped using the ``rename``
    keyword argument.

        >>> renamed = Planck18.to_format("astropy.row", rename={"H0": "Hubble"})
        >>> renamed
        <Row index=0>
          cosmology     name      Hubble      Om0    Tcmb0    Neff      m_nu      Ob0
                               km / (Mpc s)            K                 eV
            str13       str8     float64    float64 float64 float64  float64[3] float64
        ------------- -------- ------------ ------- ------- ------- ----------- -------
        FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

        >>> cosmo = Cosmology.from_format(renamed, format="astropy.row",
        ...                               rename={"Hubble": "H0"})
        >>> cosmo == Planck18
        True
    """
    inv_rename = {v: k for k, v in rename.items()} if rename is not None else {}
    kname = inv_rename.get("name", "name")
    kmeta = inv_rename.get("meta", "meta")
    kcosmo = inv_rename.get("cosmology", "cosmology")

    # special values
    name = row.get(kname)

    meta = defaultdict(dict, copy.deepcopy(row.meta))
    # Now need to add the Columnar metadata. This is only available on the
    # parent table. If Row is ever separated from Table, this should be moved
    # to ``to_table``.
    for col in row._table.itercols():
        if col.info.meta:  # Only add metadata if not empty
            meta[col.name].update(col.info.meta)

    # turn row into mapping, filling cosmo if not in a column
    mapping = dict(row)
    mapping[kname] = name
    mapping.setdefault(kcosmo, meta.pop(kcosmo, None))
    mapping[kmeta] = dict(meta)

    # build cosmology from map
    return from_mapping(
        mapping, move_to_meta=move_to_meta, cosmology=cosmology, rename=rename
    )


def to_row(
    cosmology: Cosmology,
    *args: object,
    cosmology_in_meta: bool = False,
    table_cls: type[Table] = QTable,
    rename: Mapping[str, str] | None = None,
) -> Row:
    """Serialize the cosmology into a `~astropy.table.Row`.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology`
        The cosmology instance to convert to a mapping.
    *args
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`
    table_cls : type (optional, keyword-only)
        Astropy :class:`~astropy.table.Table` class or subclass type to use. Default is
        :class:`~astropy.table.QTable`.
    cosmology_in_meta : bool
        Whether to put the cosmology class in the Table metadata (if `True`) or as the
        first column (if `False`, default).
    rename : Mapping[str, str] or None (optional, keyword-only)
        A mapping of field names of the |Cosmology| to column names in the row.

    Returns
    -------
    `~astropy.table.Row`
        With columns for the cosmology parameters, and metadata in the Table's ``meta``
        attribute. The cosmology class name will either be a column or in ``meta``,
        depending on 'cosmology_in_meta'.

    Examples
    --------
    A `~astropy.cosmology.Cosmology` as a `~astropy.table.Row` will have the cosmology's
    name and parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> cr = Planck18.to_format("astropy.row")
        >>> cr
        <Row index=0>
          cosmology     name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                               km / (Mpc s)            K                 eV
            str13       str8     float64    float64 float64 float64  float64[3] float64
        ------------- -------- ------------ ------- ------- ------- ----------- -------
        FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    The cosmological class and other metadata, e.g. a paper reference, are in the
    Table's metadata.

        >>> cr.meta
        {'Oc0': 0.2607, 'n': 0.9665, ...}

    To move the cosmology class from a column to the Table's metadata, set the
    ``cosmology_in_meta`` argument to `True`:

        >>> Planck18.to_format("astropy.table", cosmology_in_meta=True)
        <QTable length=1>
        name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                km / (Mpc s)            K                 eV
        str8     float64    float64 float64 float64  float64[3] float64
        -------- ------------ ------- ------- ------- ----------- -------
        Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    In Astropy, Row objects are always part of a Table. :class:`~astropy.table.QTable`
    is recommended for tables with |Quantity| columns. However the returned type may be
    overridden using the ``cls`` argument:

        >>> from astropy.table import Table
        >>> Planck18.to_format("astropy.table", cls=Table)
        <Table length=1>
        ...

    The columns can be renamed using the ``rename`` keyword argument.

        >>> renamed = Planck18.to_format("astropy.row", rename={"H0": "Hubble"})
        >>> renamed
        <Row index=0>
          cosmology     name      Hubble      Om0    Tcmb0    Neff      m_nu      Ob0
                               km / (Mpc s)            K                 eV
            str13       str8     float64    float64 float64 float64  float64[3] float64
        ------------- -------- ------------ ------- ------- ------- ----------- -------
        FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897
    """
    from .table import to_table

    table = to_table(
        cosmology, cls=table_cls, cosmology_in_meta=cosmology_in_meta, rename=rename
    )
    return table[0]  # extract row from table


def row_identify(
    origin: str, format: str | None, *args: object, **kwargs: object
) -> bool:
    """Identify if object uses the `~astropy.table.Row` format.

    Returns
    -------
    bool
    """
    itis = False
    if origin == "read":
        itis = isinstance(args[1], Row) and (format in (None, "astropy.row"))
    return itis


# ===================================================================
# Register

convert_registry.register_reader("astropy.row", Cosmology, from_row)
convert_registry.register_writer("astropy.row", Cosmology, to_row)
convert_registry.register_identifier("astropy.row", Cosmology, row_identify)
