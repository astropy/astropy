# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy

import numpy as np

from astropy.cosmology.connect import convert_registry
from astropy.cosmology.core import Cosmology
from astropy.table import Row

from .mapping import from_mapping


def from_row(row, *, move_to_meta=False, cosmology=None):
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

    cosmology : str, `~astropy.cosmology.Cosmology` class, or None (optional, keyword-only)
        The cosmology class (or string name thereof) to use when constructing
        the cosmology instance. The class also provides default parameter values,
        filling in any non-mandatory arguments missing in 'table'.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a Row with
    ``from_row``, we will first make a `~astropy.table.Row` using
    :func:`~astropy.cosmology.Cosmology.to_format`.

        >>> from astropy.cosmology import Cosmology, Planck18
        >>> cr = Planck18.to_format("astropy.row")
        >>> cr
        <Row index=0>
          cosmology     name        H0        Om0    Tcmb0    Neff    m_nu [3]    Ob0
                               km / (Mpc s)            K                 eV
            str13       str8     float64    float64 float64 float64   float64   float64
        ------------- -------- ------------ ------- ------- ------- ----------- -------
        FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    Now this row can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = Cosmology.from_format(cr, format="astropy.row")
        >>> cosmo
        FlatLambdaCDM(name="Planck18", H0=67.7 km / (Mpc s), Om0=0.31,
                      Tcmb0=2.725 K, Neff=3.05, m_nu=[0. 0. 0.06] eV, Ob0=0.049)
    """
    # special values
    name = row['name'] if 'name' in row.columns else None  # get name from column
    meta = copy.deepcopy(row.meta)

    # turn row into mapping, filling cosmo if not in a column
    mapping = dict(row)
    mapping["name"] = name
    mapping.setdefault("cosmology", meta.pop("cosmology", None))
    mapping["meta"] = meta

    # build cosmology from map
    return from_mapping(mapping, move_to_meta=move_to_meta, cosmology=cosmology)


def to_row(cosmology, *args, cosmology_in_meta=False):
    """Serialize the cosmology into a `~astropy.table.Row`.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    *args
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`
    cosmology_in_meta : bool
        Whether to put the cosmology class in the Table metadata (if `True`) or
        as the first column (if `False`, default).

    Returns
    -------
    `~astropy.table.Row`
        With columns for the cosmology parameters, and metadata in the Table's
        ``meta`` attribute. The cosmology class name will either be a column
        or in ``meta``, depending on 'cosmology_in_meta'.

    Examples
    --------
    A Cosmology as a `~astropy.table.Row` will have the cosmology's name and
    parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> cr = Planck18.to_format("astropy.row")
        >>> cr
        <Row index=0>
          cosmology     name        H0        Om0    Tcmb0    Neff    m_nu [3]    Ob0
                               km / (Mpc s)            K                 eV
            str13       str8     float64    float64 float64 float64   float64   float64
        ------------- -------- ------------ ------- ------- ------- ----------- -------
        FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    The cosmological class and other metadata, e.g. a paper reference, are in
    the Table's metadata.
    """
    from .table import to_table

    table = to_table(cosmology, cosmology_in_meta=cosmology_in_meta)
    return table[0]  # extract row from table


def row_identify(origin, format, *args, **kwargs):
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
