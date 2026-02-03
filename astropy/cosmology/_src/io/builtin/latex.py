r"""|Cosmology| <-> LaTeX I/O, using |Cosmology.read| and |Cosmology.write|.

We assume the following setup:

    >>> from pathlib import Path
    >>> from tempfile import TemporaryDirectory
    >>> temp_dir = TemporaryDirectory()

Writing a cosmology to a LaTeX file will produce a table with the cosmology's type,
name, and parameters as columns.

    >>> from astropy.cosmology import Cosmology, Planck18
    >>> file = Path(temp_dir.name) / "file.tex"

    >>> Planck18.write(file, format="ascii.latex")
    >>> with open(file) as f: print(f.read())
    \begin{table}
    \begin{tabular}{cccccccc}
    cosmology & name & $H_0$ & $\Omega_{m,0}$ & $T_{0}$ & $N_{eff}$ & $m_{nu}$ & $\Omega_{b,0}$ \\
    &  & $\mathrm{km\,Mpc^{-1}\,s^{-1}}$ &  & $\mathrm{K}$ &  & $\mathrm{eV}$ &  \\
    FlatLambdaCDM & Planck18 & 67.66 & 0.30966 & 2.7255 & 3.046 & 0.0 .. 0.06 & 0.04897 \\
    \end{tabular}
    \end{table}
    <BLANKLINE>

The cosmology's metadata is not included in the table.

To save the cosmology in an existing file, use ``overwrite=True``; otherwise, an
error will be raised.

    >>> Planck18.write(file, format="ascii.latex", overwrite=True)

To use a different table class as the underlying writer, use the ``cls`` kwarg. For
more information on the available table classes, see the documentation on Astropy's
table classes and on ``Cosmology.to_format("astropy.table")``.

By default the parameter names are converted to LaTeX format. To disable this, set
``latex_names=False``.

    >>> file = Path(temp_dir.name) / "file2.tex"
    >>> Planck18.write(file, format="ascii.latex", latex_names=False)
    >>> with open(file) as f: print(f.read())
    \begin{table}
    \begin{tabular}{cccccccc}
    cosmology & name & H0 & Om0 & Tcmb0 & Neff & m_nu & Ob0 \\
    &  & $\mathrm{km\,Mpc^{-1}\,s^{-1}}$ &  & $\mathrm{K}$ &  & $\mathrm{eV}$ &  \\
    FlatLambdaCDM & Planck18 & 67.66 & 0.30966 & 2.7255 & 3.046 & 0.0 .. 0.06 & 0.04897 \\
    \end{tabular}
    \end{table}
    <BLANKLINE>

.. testcleanup::

    >>> temp_dir.cleanup()
"""

from typing import Any, TypeVar

import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology._src.core import Cosmology
from astropy.cosmology._src.io.connect import readwrite_registry
from astropy.cosmology._src.parameter import Parameter
from astropy.cosmology._src.typing import _CosmoT
from astropy.io.typing import (
    PathLike,
    ReadableFileLike,  # not sure if i should place it here but let's see if it all works...!
    WriteableFileLike,
)
from astropy.table import QTable, Table

from .table import from_table, to_table

_TableT = TypeVar("_TableT", bound=Table)

_FORMAT_TABLE = {
    "H0": "$H_0$",
    "Om0": r"$\Omega_{m,0}$",
    "Ode0": r"$\Omega_{\Lambda,0}$",
    "Tcmb0": "$T_{0}$",
    "Neff": "$N_{eff}$",
    "m_nu": "$m_{nu}$",
    "Ob0": r"$\Omega_{b,0}$",
    "w0": "$w_{0}$",
    "wa": "$w_{a}$",
    "wz": "$w_{z}$",
    "wp": "$w_{p}$",
    "zp": "$z_{p}$",
}


def read_latex(
    filename: PathLike | ReadableFileLike[Table],
    index: int | str | None = None,
    *,
    move_to_meta: bool = False,
    cosmology: str | type[_CosmoT] | None = None,
    latex_names: bool = True,
    **kwargs: Any,
) -> _CosmoT:
    r"""Read a |Cosmology| from a LaTeX file.

    Parameters
    ----------
    filename : path-like or file-like
        From where to read the Cosmology.
    index : int or str or None, optional
        Needed to select the row in tables with multiple rows. ``index`` can be an
        integer for the row number or, if the table is indexed by a column, the value of
        that column. If the table is not indexed and ``index`` is a string, the "name"
        column is used as the indexing column.

    move_to_meta : bool, optional keyword-only
        Whether to move keyword arguments that are not in the Cosmology class' signature
        to the Cosmology's metadata. This will only be applied if the Cosmology does NOT
        have a keyword-only argument (e.g. ``**kwargs``). Arguments moved to the
        metadata will be merged with existing metadata, preferring specified metadata in
        the case of a merge conflict (e.g. for ``Cosmology(meta={'key':10}, key=42)``,
        the ``Cosmology.meta`` will be ``{'key': 10}``).
    cosmology : str or |Cosmology| class or None, optional keyword-only
        The cosmology class (or string name thereof) to use when constructing the
        cosmology instance. The class also provides default parameter values, filling in
        any non-mandatory arguments missing in 'table'.
    latex_names : bool, optional keyword-only
        Whether the |Table| (might) have latex column names for the parameters that need
        to be mapped to the correct parameter name -- e.g. $H_0$ to 'H0'. This is
        `True` by default, but can be turned off (set to `False`) if there is a known
        name conflict (e.g. both an 'H0' and '$H_0$' column) as this will raise an
        error. In this case, the correct name ('H0') is preferred.
    **kwargs : Any
        Passed to ``QTable.read``. ``format`` is set to 'ascii.latex', regardless of
        input.

    Returns
    -------
    |Cosmology| subclass instance

    Raises
    ------
    ValueError
        If the keyword argument 'format' is given and is not "ascii.latex".
    """
    # Check that the format is 'ascii.latex' (or not specified)
    fmt = kwargs.pop("format", "ascii.latex")
    if fmt != "ascii.latex":
        raise ValueError(f"format must be 'ascii.latex', not {fmt}")

    # Reading is handled by `QTable`.
    with u.add_enabled_units(cu):
        table = QTable.read(filename, format="ascii.latex", **kwargs)

    # No need of units of different cosmology parameters to support cosmology conversion
    del table[0]

    if latex_names:
        table_columns = set(table.colnames)
        for name, latex in _FORMAT_TABLE.items():
            if latex in table_columns:
                table.rename_column(latex, name)

    return from_table(
        table, index, move_to_meta=move_to_meta, cosmology=cosmology, rename=None
    )


def write_latex(
    cosmology: Cosmology,
    file: PathLike | WriteableFileLike[_TableT],
    *,
    overwrite: bool = False,
    cls: type[_TableT] = QTable,
    latex_names: bool = True,
    **kwargs: Any,
) -> None:
    r"""Serialize the |Cosmology| into a LaTeX.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
        The cosmology to serialize.
    file : path-like or file-like
        Location to save the serialized cosmology.

    overwrite : bool
        Whether to overwrite the file, if it exists.
    cls : type, optional keyword-only
        Astropy :class:`~astropy.table.Table` (sub)class to use when writing. Default is
        :class:`~astropy.table.QTable`.
    latex_names : bool, optional keyword-only
        Whether to use LaTeX names for the parameters. Default is `True`.
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

    Writing a cosmology to a LaTeX file will produce a table with the cosmology's type,
    name, and parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> file = Path(temp_dir.name) / "file.tex"

        >>> Planck18.write(file, format="ascii.latex")
        >>> with open(file) as f: print(f.read())
        \begin{table}
        \begin{tabular}{cccccccc}
        cosmology & name & $H_0$ & $\Omega_{m,0}$ & $T_{0}$ & $N_{eff}$ & $m_{nu}$ & $\Omega_{b,0}$ \\
        &  & $\mathrm{km\,Mpc^{-1}\,s^{-1}}$ &  & $\mathrm{K}$ &  & $\mathrm{eV}$ &  \\
        FlatLambdaCDM & Planck18 & 67.66 & 0.30966 & 2.7255 & 3.046 & 0.0 .. 0.06 & 0.04897 \\
        \end{tabular}
        \end{table}
        <BLANKLINE>

    The cosmology's metadata is not included in the table.

    To save the cosmology in an existing file, use ``overwrite=True``; otherwise, an
    error will be raised.

        >>> Planck18.write(file, format="ascii.latex", overwrite=True)

    To use a different table class as the underlying writer, use the ``cls`` kwarg. For
    more information on the available table classes, see the documentation on Astropy's
    table classes and on ``Cosmology.to_format("astropy.table")``.

    By default the parameter names are converted to LaTeX format. To disable this, set
    ``latex_names=False``.

        >>> file = Path(temp_dir.name) / "file2.tex"
        >>> Planck18.write(file, format="ascii.latex", latex_names=False)
        >>> with open(file) as f: print(f.read())
        \begin{table}
        \begin{tabular}{cccccccc}
        cosmology & name & H0 & Om0 & Tcmb0 & Neff & m_nu & Ob0 \\
        &  & $\mathrm{km\,Mpc^{-1}\,s^{-1}}$ &  & $\mathrm{K}$ &  & $\mathrm{eV}$ &  \\
        FlatLambdaCDM & Planck18 & 67.66 & 0.30966 & 2.7255 & 3.046 & 0.0 .. 0.06 & 0.04897 \\
        \end{tabular}
        \end{table}
        <BLANKLINE>

    .. testcleanup::

        >>> temp_dir.cleanup()
    """
    # Check that the format is 'latex', 'ascii.latex' (or not specified)
    fmt = kwargs.pop("format", "ascii.latex")
    if fmt != "ascii.latex":
        raise ValueError(f"format must be 'ascii.latex', not {fmt}")

    # Set cosmology_in_meta as false for now since there is no metadata being kept
    table = to_table(cosmology, cls=cls, cosmology_in_meta=False)

    cosmo_cls = type(cosmology)
    for name in table.columns.keys():
        param = cosmo_cls.parameters.get(name)
        if not isinstance(param, Parameter) or param.unit in (None, u.one):
            continue
        # Get column to correct unit
        table[name] <<= param.unit

    # Convert parameter names to LaTeX format
    if latex_names:
        new_names = [_FORMAT_TABLE.get(k, k) for k in cosmology.parameters]
        table.rename_columns(tuple(cosmology.parameters), new_names)

    table.write(file, overwrite=overwrite, format="ascii.latex", **kwargs)


def latex_identify(
    origin: object, filepath: str | None, *args: object, **kwargs: object
) -> bool:
    """Identify if object uses the Table format.

    Returns
    -------
    bool
    """
    return filepath is not None and filepath.endswith(".tex")


# ===================================================================
# Register

readwrite_registry.register_reader("ascii.latex", Cosmology, read_latex)
readwrite_registry.register_writer("ascii.latex", Cosmology, write_latex)
readwrite_registry.register_identifier("ascii.latex", Cosmology, latex_identify)
