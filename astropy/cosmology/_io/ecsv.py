# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology.connect import readwrite_registry
from astropy.cosmology.core import Cosmology
from astropy.table import QTable

from .table import from_table, to_table


def read_ecsv(
    filename, index=None, *, move_to_meta=False, cosmology=None, rename=None, **kwargs
):
    """Read a `~astropy.cosmology.Cosmology` from an ECSV file.

    Parameters
    ----------
    filename : path-like or file-like
        From where to read the Cosmology.
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

    rename : dict or None (optional keyword-only)
        A dictionary mapping column names to fields of the
        `~astropy.cosmology.Cosmology`.

    **kwargs
        Passed to :attr:`astropy.table.QTable.read`

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance
    """
    kwargs["format"] = "ascii.ecsv"
    with u.add_enabled_units(cu):
        table = QTable.read(filename, **kwargs)

    # build cosmology from table
    return from_table(
        table,
        index=index,
        move_to_meta=move_to_meta,
        cosmology=cosmology,
        rename=rename,
    )


def write_ecsv(
    cosmology,
    file,
    *,
    overwrite=False,
    cls=QTable,
    cosmology_in_meta=True,
    rename=None,
    **kwargs
):
    """Serialize the cosmology into a ECSV.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology`
        The cosmology instance to convert to a mapping.
    file : path-like or file-like
        Location to save the serialized cosmology.

    overwrite : bool
        Whether to overwrite the file, if it exists.
    cls : type (optional, keyword-only)
        Astropy :class:`~astropy.table.Table` (sub)class to use when writing.
        Default is :class:`~astropy.table.QTable`.
    cosmology_in_meta : bool (optional, keyword-only)
        Whether to put the cosmology class in the Table metadata (if `True`,
        default) or as the first column (if `False`).
    rename : dict or None (optional keyword-only)
        A dictionary mapping fields of the `~astropy.cosmology.Cosmology` to
        columns of the table.

    **kwargs
        Passed to ``cls.write``

    Raises
    ------
    TypeError
        If kwarg (optional) 'cls' is not a subclass of `astropy.table.Table`
    """
    table = to_table(
        cosmology, cls=cls, cosmology_in_meta=cosmology_in_meta, rename=rename
    )

    kwargs["format"] = "ascii.ecsv"
    table.write(file, overwrite=overwrite, **kwargs)


def ecsv_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if object uses the Table format.

    Returns
    -------
    bool
    """
    return filepath is not None and filepath.endswith(".ecsv")


# ===================================================================
# Register

readwrite_registry.register_reader("ascii.ecsv", Cosmology, read_ecsv)
readwrite_registry.register_writer("ascii.ecsv", Cosmology, write_ecsv)
readwrite_registry.register_identifier("ascii.ecsv", Cosmology, ecsv_identify)
