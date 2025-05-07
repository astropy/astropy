r"""|Cosmology| <-> MRT I/O, using |Cosmology.read| and |Cosmology.write|.

This module provides functions to write/read a |Cosmology| object to/from an MRT file.
The functions are registered with ``readwrite_registry`` under the format name "ascii.mrt".

"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Any, TypeVar

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology._src.core import Cosmology
from astropy.cosmology._src.io.connect import readwrite_registry
from astropy.io.typing import PathLike, ReadableFileLike, WriteableFileLike
from astropy.table import QTable

from .table import from_table, to_table

if TYPE_CHECKING:
    from astropy.cosmology._src.typing import _CosmoT
    from astropy.io.typing import PathLike, ReadableFileLike, WriteableFileLike
    from astropy.table import Table

    _TableT = TypeVar("_TableT", "Table")


def read_mrt(
    file: PathLike | ReadableFileLike[Table],
    index: int | str | None = None,
    *,
    move_to_meta: bool = False,
    cosmology: str | type[_CosmoT] | None = None,
    **kwargs: Any,
) -> _CosmoT:
    r"""Read a `~astropy.cosmology.Cosmology` from an MRT file.

    Parameters
    ----------

    Examples
    --------
    We assume the following setup:

        >>> from pathlib import Path
        >>> from tempfile import TemporaryDirectory
        >>> temp_dir = TemporaryDirectory()

    Writing a cosmology to a LaTeX file will produce a table with the cosmology's type,
    name, and parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> file = Path(temp_dir.name) / "file.mrt"

        >>> Planck18.write(file, format="ascii.mrt")
        >>> with open(file) as f: print(f.read())




    """
    format = kwargs.pop("format", "ascii.mrt")
    if format != "ascii.mrt":
        raise ValueError(f"format must be 'ascii.mrt',not {format}")

    with u.add_enabled_units(cu):
        table = QTable.read(file, format="ascii.mrt", **kwargs)

    # Build the cosmology from table, using the private backend.
    return from_table(
        table, index=index, move_to_meta=move_to_meta, cosmology=cosmology
    )


def write_mrt(
    cosmology: Cosmology,
    file: PathLike | WriteableFileLike[_TableT],
    *,
    overwrite: bool = False,
    cls: type[_TableT] = QTable,
    cosmology_in_meta: bool = True,
    **kwargs: Any,
):
    format = kwargs.pop("format", "ascii.mrt")
    if format != "ascii.mrt":
        raise ValueError(f"format must be 'ascii.mrt', not {format}")

    table = to_table(cosmology, cls=cls, cosmology_in_meta=cosmology_in_meta)
    for k, col in table.columns.items():
        # CDS can't serialize redshift units, so remove them  # TODO: fix this
        if col.unit is cu.redshift:
            table[k] <<= u.dimensionless_unscaled
        ## check if col is mutil den
        if len(col.shape) > 1:

            def format_col_item(idx):
                obj = col[idx]
                try:
                    obj = obj.value.tolist() if hasattr(obj, "value") else obj.tolist()
                except AttributeError as exc:
                    pass
                return json.dumps(obj, separators=(",", ":"))
        else:

            def format_col_item(idx):
                return str(col[idx])

        try:
            table[k] = [format_col_item(idx) for idx in range(len(col))]

        except TypeError as exc:
            raise TypeError(
                f"could not convert column {col.info.name!r} to string: {exc}"
            ) from exc

    # Write MRT
    table.write(file, overwrite=overwrite, format="ascii.mrt", **kwargs)


def mrt_identify(
    origin: object, filepath: object, *args: object, **kwargs: object
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
    return isinstance(filepath, str) and filepath.endswith(".mrt")


readwrite_registry.register_reader("ascii.mrt", Cosmology, read_mrt)
readwrite_registry.register_writer("ascii.mrt", Cosmology, write_mrt)
readwrite_registry.register_identifier("ascii.mrt", Cosmology, mrt_identify)
