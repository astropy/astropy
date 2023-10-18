# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""|Cosmology| <-> ECSV I/O, using |Cosmology.read| and |Cosmology.write|.

This module provides functions to write/read a |Cosmology| object to/from an ECSV file.
The functions are registered with ``readwrite_registry`` under the format name
"ascii.ecsv".

We assume the following setup:

    >>> from pathlib import Path
    >>> from tempfile import TemporaryDirectory
    >>> temp_dir = TemporaryDirectory()

To see reading a Cosmology from an ECSV file, we first write a Cosmology to an ECSV
file:

    >>> from astropy.cosmology import Cosmology, Planck18
    >>> file = Path(temp_dir.name) / "file.ecsv"
    >>> Planck18.write(file)

    >>> with open(file) as f: print(f.read())
    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: name, datatype: string}
    ...
    # meta: !!omap
    # - {Oc0: 0.2607}
    ...
    # schema: astropy-2.0
    name H0 Om0 Tcmb0 Neff m_nu Ob0
    Planck18 67.66 0.30966 2.7255 3.046 [0.0,0.0,0.06] 0.04897
    <BLANKLINE>

Now we can read the Cosmology from the ECSV file, constructing a new cosmological
instance identical to the ``Planck18`` cosmology from which it was generated.

    >>> cosmo = Cosmology.read(file)
    >>> print(cosmo)
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
    >>> cosmo == Planck18
    True

If a file already exists, attempting to write will raise an error unless
``overwrite=True``.

    >>> Planck18.write(file, overwrite=True)

By default the cosmology class is written to the Table metadata. This can be changed to
a column of the table using the ``cosmology_in_meta`` keyword argument.

    >>> file = Path(temp_dir.name) / "file2.ecsv"
    >>> Planck18.write(file, cosmology_in_meta=False)
    >>> with open(file) as f: print(f.read())
    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: cosmology, datatype: string}
    # - {name: name, datatype: string}
    ...
    # meta: !!omap
    # - {Oc0: 0.2607}
    ...
    # schema: astropy-2.0
    cosmology name H0 Om0 Tcmb0 Neff m_nu Ob0
    FlatLambdaCDM Planck18 67.66 0.30966 2.7255 3.046 [0.0,0.0,0.06] 0.04897
    <BLANKLINE>

The ``cosmology`` information (column or metadata) may be omitted if the cosmology class
(or its string name) is passed as the ``cosmology`` keyword argument to
|Cosmology.read|. Alternatively, specific cosmology classes can be used to parse the
data.

    >>> from astropy.cosmology import FlatLambdaCDM
    >>> print(FlatLambdaCDM.read(file))
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                    Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)

When using a specific cosmology class, the class' default parameter values are used to
fill in any missing information.

For files with multiple rows of cosmological parameters, the ``index`` argument is
needed to select the correct row. The index can be an integer for the row number or, if
the table is indexed by a column, the value of that column. If the table is not indexed
and ``index`` is a string, the "name" column is used as the indexing column.

Here is an example where ``index`` is needed and can be either an integer (for the row
number) or the name of one of the cosmologies, e.g. 'Planck15'.

    >>> from astropy.cosmology import Planck13, Planck15, Planck18
    >>> from astropy.table import vstack
    >>> cts = vstack([c.to_format("astropy.table")
    ...               for c in (Planck13, Planck15, Planck18)],
    ...              metadata_conflicts='silent')
    >>> file = Path(temp_dir.name) / "file3.ecsv"
    >>> cts.write(file)
    >>> with open(file) as f: print(f.read())
    # %ECSV 1.0
    # ---
    # datatype:
    # - {name: name, datatype: string}
    ...
    # meta: !!omap
    # - {Oc0: 0.2607}
    ...
    # schema: astropy-2.0
    name H0 Om0 Tcmb0 Neff m_nu Ob0
    Planck13 67.77 0.30712 2.7255 3.046 [0.0,0.0,0.06] 0.048252
    Planck15 67.74 0.3075 2.7255 3.046 [0.0,0.0,0.06] 0.0486
    Planck18 67.66 0.30966 2.7255 3.046 [0.0,0.0,0.06] 0.04897

    >>> cosmo = Cosmology.read(file, index="Planck15", format="ascii.ecsv")
    >>> cosmo == Planck15
    True

Fields of the table in the file can be renamed to match the
`~astropy.cosmology.Cosmology` class' signature using the ``rename`` argument. This is
useful when the files's column names do not match the class' parameter names.

    >>> file = Path(temp_dir.name) / "file4.ecsv"
    >>> Planck18.write(file, rename={"H0": "Hubble"})
    >>> with open(file) as f: print(f.read())
     # %ECSV 1.0
    # ---
    # datatype:
    # - {name: name, datatype: string}
    ...
    # meta: !!omap
    # - {Oc0: 0.2607}
    ...
    # schema: astropy-2.0
    name Hubble Om0 Tcmb0 Neff m_nu Ob0
    ...

    >>> cosmo = Cosmology.read(file, rename={"Hubble": "H0"})
    >>> cosmo == Planck18
    True

By default :class:`~astropy.cosmology.Cosmology` instances are written using
`~astropy.table.QTable` as an intermediate representation (for details see
|Cosmology.to_format|, with ``format="astropy.table"``). The `~astropy.table.Table` type
can be changed using the ``cls`` keyword argument.

    >>> from astropy.table import Table
    >>> file = Path(temp_dir.name) / "file5.ecsv"
    >>> Planck18.write(file, cls=Table)

For most use cases, the default ``cls`` of :class:`~astropy.table.QTable` is recommended
and will be largely indistinguishable from other table types, as the ECSV format is
agnostic to the table type. An example of a difference that might necessitate using a
different table type is if a different ECSV schema is desired.

Additional keyword arguments are passed to ``QTable.read`` and ``QTable.write``.

.. testcleanup::

    >>> temp_dir.cleanup()
"""

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology.connect import readwrite_registry
from astropy.cosmology.core import Cosmology
from astropy.table import QTable

from .table import from_table, to_table


def read_ecsv(
    filename, index=None, *, move_to_meta=False, cosmology=None, rename=None, **kwargs
):
    r"""Read a `~astropy.cosmology.Cosmology` from an ECSV file.

    Parameters
    ----------
    filename : path-like or file-like
        From where to read the Cosmology.
    index : int, str, or None, optional
        Needed to select the row in tables with multiple rows. ``index`` can be an
        integer for the row number or, if the table is indexed by a column, the value of
        that column. If the table is not indexed and ``index`` is a string, the "name"
        column is used as the indexing column.

    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class' signature
        to the Cosmology's metadata. This will only be applied if the Cosmology does NOT
        have a keyword-only argument (e.g. ``**kwargs``). Arguments moved to the
        metadata will be merged with existing metadata, preferring specified metadata in
        the case of a merge conflict (e.g. for ``Cosmology(meta={'key':10}, key=42)``,
        the ``Cosmology.meta`` will be ``{'key': 10}``).

    rename : dict or None (optional keyword-only)
        A dictionary mapping column names to fields of the
        `~astropy.cosmology.Cosmology`.

    **kwargs
        Passed to ``QTable.read``

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    We assume the following setup:

        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> temp_dir = TemporaryDirectory()

    To see reading a Cosmology from an ECSV file, we first write a Cosmology to an ECSV
    file:

        >>> from astropy.cosmology import Cosmology, Planck18
        >>> file = Path(temp_dir.name) / "file.ecsv"
        >>> Planck18.write(file)

        >>> with open(file) as f: print(f.read())
        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: name, datatype: string}
        ...
        # meta: !!omap
        # - {Oc0: 0.2607}
        ...
        # schema: astropy-2.0
        name H0 Om0 Tcmb0 Neff m_nu Ob0
        Planck18 67.66 0.30966 2.7255 3.046 [0.0,0.0,0.06] 0.04897
        <BLANKLINE>

    Now we can read the Cosmology from the ECSV file, constructing a new cosmological
    instance identical to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = Cosmology.read(file)
        >>> print(cosmo)
        FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                    Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)
        >>> cosmo == Planck18
        True

    The ``cosmology`` information (column or metadata) may be omitted if the cosmology
    class (or its string name) is passed as the ``cosmology`` keyword argument to
    |Cosmology.read|. Alternatively, specific cosmology classes can be used to parse the
    data.

        >>> from astropy.cosmology import FlatLambdaCDM
        >>> print(FlatLambdaCDM.read(file))
        FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                      Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)

    When using a specific cosmology class, the class' default parameter values are used
    to fill in any missing information.

    For files with multiple rows of cosmological parameters, the ``index`` argument is
    needed to select the correct row. The index can be an integer for the row number or,
    if the table is indexed by a column, the value of that column. If the table is not
    indexed and ``index`` is a string, the "name" column is used as the indexing column.

    Here is an example where ``index`` is needed and can be either an integer (for the
    row number) or the name of one of the cosmologies, e.g. 'Planck15'.

        >>> from astropy.cosmology import Planck13, Planck15, Planck18
        >>> from astropy.table import vstack
        >>> cts = vstack([c.to_format("astropy.table")
        ...               for c in (Planck13, Planck15, Planck18)],
        ...              metadata_conflicts='silent')
        >>> file = Path(temp_dir.name) / "file2.ecsv"
        >>> cts.write(file)
        >>> with open(file) as f: print(f.read())
        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: name, datatype: string}
        ...
        # meta: !!omap
        # - {Oc0: 0.2607}
        ...
        # schema: astropy-2.0
        name H0 Om0 Tcmb0 Neff m_nu Ob0
        Planck13 67.77 0.30712 2.7255 3.046 [0.0,0.0,0.06] 0.048252
        Planck15 67.74 0.3075 2.7255 3.046 [0.0,0.0,0.06] 0.0486
        Planck18 67.66 0.30966 2.7255 3.046 [0.0,0.0,0.06] 0.04897

        >>> cosmo = Cosmology.read(file, index="Planck15", format="ascii.ecsv")
        >>> cosmo == Planck15
        True

    Fields of the table in the file can be renamed to match the
    `~astropy.cosmology.Cosmology` class' signature using the ``rename`` argument. This
    is useful when the files's column names do not match the class' parameter names.
    For this example we need to make a new file with renamed columns:

        >>> file = Path(temp_dir.name) / "file3.ecsv"
        >>> renamed_table = Planck18.to_format("astropy.table", rename={"H0": "Hubble"})
        >>> renamed_table.write(file)
        >>> with open(file) as f: print(f.read())
        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: name, datatype: string}
        ...
        # meta: !!omap
        # - {Oc0: 0.2607}
        ...
        # schema: astropy-2.0
        name Hubble Om0 Tcmb0 Neff m_nu Ob0
        ...

    Now we can read the Cosmology from the ECSV file, with the required renaming:

        >>> cosmo = Cosmology.read(file, rename={"Hubble": "H0"})
        >>> cosmo == Planck18
        True

    Additional keyword arguments are passed to ``QTable.read``.

    .. testcleanup::

        >>> temp_dir.cleanup()
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
    **kwargs,
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
        Astropy :class:`~astropy.table.Table` (sub)class to use when writing. Default is
        :class:`~astropy.table.QTable`.
    cosmology_in_meta : bool (optional, keyword-only)
        Whether to put the cosmology class in the Table metadata (if `True`, default) or
        as the first column (if `False`).
    rename : dict or None (optional keyword-only)
        A dictionary mapping fields of the `~astropy.cosmology.Cosmology` to columns of
        the table.

    **kwargs
        Passed to ``cls.write``

    Raises
    ------
    TypeError
        If kwarg (optional) 'cls' is not a subclass of `astropy.table.Table`

    Examples
    --------
    We assume the following setup:

        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> temp_dir = TemporaryDirectory()

    A Cosmology can be written to an ECSV file:

        >>> from astropy.cosmology import Cosmology, Planck18
        >>> file = Path(temp_dir.name) / "file.ecsv"
        >>> Planck18.write(file)

        >>> with open(file) as f: print(f.read())
        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: name, datatype: string}
        ...
        # meta: !!omap
        # - {Oc0: 0.2607}
        ...
        # schema: astropy-2.0
        name H0 Om0 Tcmb0 Neff m_nu Ob0
        Planck18 67.66 0.30966 2.7255 3.046 [0.0,0.0,0.06] 0.04897
        <BLANKLINE>

    If a file already exists, attempting to write will raise an error unless
    ``overwrite=True``.

        >>> Planck18.write(file, overwrite=True)

    By default :class:`~astropy.cosmology.Cosmology` instances are written using
    `~astropy.table.QTable` as an intermediate representation (for details see
    |Cosmology.to_format|, with ``format="astropy.table"``). The `~astropy.table.Table`
    type can be changed using the ``cls`` keyword argument.

        >>> from astropy.table import Table
        >>> file = Path(temp_dir.name) / "file2.ecsv"
        >>> Planck18.write(file, cls=Table)

    For most use cases, the default ``cls`` of :class:`~astropy.table.QTable` is
    recommended and will be largely indistinguishable from other table types, as the
    ECSV format is agnostic to the table type. An example of a difference that might
    necessitate using a different table type is if a different ECSV schema is desired.

    By default the cosmology class is written to the Table metadata. This can be changed
    to a column of the table using the ``cosmology_in_meta`` keyword argument.

        >>> file = Path(temp_dir.name) / "file3.ecsv"
        >>> Planck18.write(file, cosmology_in_meta=False)
        >>> with open(file) as f: print(f.read())
        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: cosmology, datatype: string}
        # - {name: name, datatype: string}
        ...
        # meta: !!omap
        # - {Oc0: 0.2607}
        ...
        # schema: astropy-2.0
        cosmology name H0 Om0 Tcmb0 Neff m_nu Ob0
        FlatLambdaCDM Planck18 67.66 0.30966 2.7255 3.046 [0.0,0.0,0.06] 0.04897
        <BLANKLINE>

    Fields of the Cosmology can be renamed to when writing to an ECSV file using the
    ``rename`` argument.

        >>> file = Path(temp_dir.name) / "file4.ecsv"
        >>> Planck18.write(file, rename={"H0": "Hubble"})
        >>> with open(file) as f: print(f.read())
        # %ECSV 1.0
        # ---
        # datatype:
        # - {name: name, datatype: string}
        ...
        # meta: !!omap
        # - {Oc0: 0.2607}
        ...
        # schema: astropy-2.0
        name Hubble Om0 Tcmb0 Neff m_nu Ob0
        ...

    Additional keyword arguments are passed to :attr:`astropy.table.QTable.write`.

    .. testcleanup::

        >>> temp_dir.cleanup()
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
