r"""|Cosmology| <-> MRT I/O, using |Cosmology.read| and |Cosmology.write|.

We assume the following setup:

    >>> from pathlib import Path
    >>> from tempfile import TemporaryDirectory
    >>> temp_dir = TemporaryDirectory()

Writing a cosmology to a mrt file will produce a table with cosmology's type, name, and parameters as columns.

    >>> from astropy.cosmology import Planck18
    >>> file = Path(temp_dir.name) / "file.mrt"
    >>> Planck18.write(file)
    >>> with open(file) as f: print(f.read())
    Title:
    Authors:
    Table:
    ================================================================================
    Byte-by-byte Description of file: table.dat
    --------------------------------------------------------------------------------
     Bytes Format Units  Label     Explanations
    --------------------------------------------------------------------------------
     1- 8  A8           ---    name    Description of name
    10-14  F5.2   km.Mpc-1.s-1 H0      [67.66/67.66] Hubble constant at z=0.
    16-22  F7.5         ---    Om0     [0.3/0.31] Omega matter; matter
                                density/critical density at z=0.
    24-29  F6.4         K      Tcmb0   [2.72/2.73] Temperature of the CMB at z=0.
    31-35  F5.3         ---    Neff    [3.04/3.05] Number of effective neutrino
                                species.
    37-50  A14          ---    m_nu    [[0.   0.   0.06]] Mass of neutrino species.
    52-58  F7.5         ---    Ob0     [0.04/0.05] Omega baryon; baryonic matter
                                density/critical density at z=0.
    --------------------------------------------------------------------------------
    ...


.. testcleanup::

    >>> temp_dir.cleanup()
"""

__all__ = ("mrt_identify", "read_mrt", "write_mrt")

import contextlib
import json
from typing import Any, TypeVar

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology._src.core import Cosmology
from astropy.cosmology._src.io.connect import readwrite_registry
from astropy.cosmology._src.typing import _CosmoT
from astropy.io.typing import PathLike, ReadableFileLike, WriteableFileLike
from astropy.table import Column, QTable, Table

from .table import from_table, to_table

_TableT = TypeVar("_TableT", bound=Table)


def read_mrt(
    filename: PathLike | ReadableFileLike[Table],
    /,
    index: int | str | None = None,
    *,
    move_to_meta: bool = False,
    cosmology: str | type[_CosmoT] | None = None,
    **kwargs: Any,
) -> _CosmoT:
    r"""Read a `~astropy.cosmology.Cosmology` from an MRT file.

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

    cosmology : str or type or None (optional, keyword-only)
        The cosmology class (or string name thereof) to use when constructing the
        cosmology instance. The class also provides default parameter values, filling in
        any non-mandatory arguments missing in 'table'.

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

    Writing a cosmology to a Mrt file will produce a table with the cosmology's type,
    name, and parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> file = Path(temp_dir.name) / "file.mrt"

        >>> Planck18.write(file, format="ascii.mrt")
        >>> with open(file) as f: print(f.read())
        Title:
        Authors:
        Table:
        ================================================================================
        Byte-by-byte Description of file: table.dat
        --------------------------------------------------------------------------------
         Bytes Format Units  Label     Explanations
        --------------------------------------------------------------------------------
         1- 8  A8           ---    name    Description of name
        10-14  F5.2   km.Mpc-1.s-1 H0      [67.66/67.66] Hubble constant at z=0.
        16-22  F7.5         ---    Om0     [0.3/0.31] Omega matter; matter
                                    density/critical density at z=0.
        24-29  F6.4         K      Tcmb0   [2.72/2.73] Temperature of the CMB at z=0.
        31-35  F5.3         ---    Neff    [3.04/3.05] Number of effective neutrino
                                    species.
        37-50  A14          ---    m_nu    [[0.   0.   0.06]] Mass of neutrino species.
        52-58  F7.5         ---    Ob0     [0.04/0.05] Omega baryon; baryonic matter
                                    density/critical density at z=0.
        --------------------------------------------------------------------------------
        ...

    .. testcleanup::

    >>> temp_dir.cleanup()
    """
    if (fmt := kwargs.pop("format", "ascii.mrt")) != "ascii.mrt":
        raise ValueError(f"format must be 'ascii.mrt',not {fmt}")

    with u.add_enabled_units(cu):
        table = QTable.read(filename, format="ascii.mrt", **kwargs)

    # Build the cosmology from table, using the private backend.
    return from_table(
        table, index=index, move_to_meta=move_to_meta, cosmology=cosmology
    )


def write_mrt(
    cosmo: Cosmology,
    /,
    file: PathLike | WriteableFileLike[_TableT],
    *,
    overwrite: bool = False,
    cls: type[_TableT] = QTable,
    cosmology_in_meta: bool = True,
    **kwargs: Any,
):
    r"""Serialize the |Cosmology| into a MRT table.

    Parameters
    ----------
    cosmology : |Cosmology| subclass instance
        The cosmology to serialize.
    file : path-like or file-like
        Where to write the html table.
    overwrite : bool, optional keyword-only
        Whether to overwrite the file, if it exists.
    cls : |Table| class, optional keyword-onlyl
        Astropy |Table| (sub)class to use when writing. Default is |QTable| class.
    latex_names : bool, optional keyword-only
        Whether to format the parameters (column) names to latex -- e.g. 'H0' to
        $$H_{0}$$.
    **kwargs : Any
        Passed to ``cls.write``.

    Raises
    ------
    TypeError
        If the optional keyword-argument 'cls' is not a subclass of |Table|.
    ValueError
        If the keyword argument 'format' is given and is not "ascii.html".

    Examples
    --------
    We assume the following setup:

        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> temp_dir = TemporaryDirectory()

    Writing a cosmology to a html file will produce a table with the cosmology's type,
    name, and parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> file = Path(temp_dir.name) / "file.mrt"
        >>> Planck18.write(file, overwrite=True)
        >>> with open(file) as f: print(f.read())
        Title:
        Authors:
        Table:
        ================================================================================
        Byte-by-byte Description of file: table.dat
        --------------------------------------------------------------------------------
        Bytes Format Units  Label     Explanations
        --------------------------------------------------------------------------------
        1- 8  A8           ---    name    Description of name
        10-14  F5.2   km.Mpc-1.s-1 H0      [67.66/67.66] Hubble constant at z=0.
        16-22  F7.5         ---    Om0     [0.3/0.31] Omega matter; matter
                                    density/critical density at z=0.
        24-29  F6.4         K      Tcmb0   [2.72/2.73] Temperature of the CMB at z=0.
        31-35  F5.3         ---    Neff    [3.04/3.05] Number of effective neutrino
                                    species.
        37-50  A14          ---    m_nu    [[0.   0.   0.06]] Mass of neutrino species.
        52-58  F7.5         ---    Ob0     [0.04/0.05] Omega baryon; baryonic matter
                                    density/critical density at z=0.
        --------------------------------------------------------------------------------
        ...


    .. testcleanup::

        >>> temp_dir.cleanup()

    Notes
    -----

    """
    if (fmt := kwargs.pop("format", "ascii.mrt")) != "ascii.mrt":
        raise ValueError(f"format must be 'ascii.mrt', not {fmt}")

    table = to_table(cosmo, cls=cls, cosmology_in_meta=cosmology_in_meta)

    Parameters = cosmo.__class__.parameters  # dict of Parameter objects
    for name, col in table.columns.items():
        # MRT can't serialize redshift units, so remove them
        if col.unit is cu.redshift:
            table[name] <<= u.dimensionless_unscaled

        # check if col is multi dimensional
        if len(col.shape) > 1 or col.info.dtype.kind == "0":

            def format_col_item(idx):
                obj = col[idx]
                # Get the value in the default units
                if hasattr(obj, "to_value"):
                    obj = obj.to_value(Parameters[name].unit)
                with contextlib.suppress(AttributeError):
                    obj = obj.tolist()

                return json.dumps(obj, separators=(",", ":"))

            try:
                table[name] = Column(
                    data=[format_col_item(idx) for idx in range(len(col))],
                    name=name,
                    description=str(col.value) + " " + col.info.description,
                )
            except TypeError as exc:
                msg = f"could not convert column {col.info.name!r} to string: {exc}"
                raise TypeError(msg) from exc

    # Write MRT
    table.write(file, overwrite=overwrite, format="ascii.mrt", **kwargs)


def mrt_identify(
    _: object, filepath: str | None, /, *args: object, **kwargs: object
) -> bool:
    """Identify if an object uses the HTML Table format.

    Parameters
    ----------
    origin : object
        Not used.
    filepath : object
        From where to read the Cosmology.
    *args : object
        Not used.
    **kwargs : object
        Not used.

    Returns
    -------
    bool
        If the filepath is a string ending with '.mrt'.
    """
    return filepath is not None and filepath.endswith(".mrt")


# ===================================================================
# Register

readwrite_registry.register_reader("ascii.mrt", Cosmology, read_mrt)
readwrite_registry.register_writer("ascii.mrt", Cosmology, write_mrt)
readwrite_registry.register_identifier("ascii.mrt", Cosmology, mrt_identify)
