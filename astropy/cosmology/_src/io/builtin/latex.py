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

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

import astropy.units as u
from astropy.table import QTable

# isort: split
from astropy.cosmology._src.core import Cosmology
from astropy.cosmology._src.io.connect import readwrite_registry
from astropy.cosmology._src.parameter import Parameter

from .table import to_table

if TYPE_CHECKING:
    from astropy.io.typing import PathLike, WriteableFileLike
    from astropy.table import Table

    _TableT = TypeVar("_TableT", Table)

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

readwrite_registry.register_writer("ascii.latex", Cosmology, write_latex)
readwrite_registry.register_identifier("ascii.latex", Cosmology, latex_identify)
