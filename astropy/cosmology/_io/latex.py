import astropy.units as u
from astropy.cosmology.connect import readwrite_registry
from astropy.cosmology.core import Cosmology
from astropy.cosmology.parameter import Parameter
from astropy.table import QTable

from .table import to_table

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
    cosmology, file, *, overwrite=False, cls=QTable, latex_names=True, **kwargs
):
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
        Astropy :class:`~astropy.table.Table` (sub)class to use when writing.
        Default is :class:`~astropy.table.QTable`.
    latex_names : bool, optional keyword-only
        Whether to use LaTeX names for the parameters. Default is `True`.
    **kwargs
        Passed to ``cls.write``

    Raises
    ------
    TypeError
        If kwarg (optional) 'cls' is not a subclass of `astropy.table.Table`
    """
    # Check that the format is 'latex', 'ascii.latex' (or not specified)
    fmt = kwargs.pop("format", "ascii.latex")
    if fmt != "ascii.latex":
        raise ValueError(f"format must be 'ascii.latex', not {fmt}")

    # Set cosmology_in_meta as false for now since there is no metadata being kept
    table = to_table(cosmology, cls=cls, cosmology_in_meta=False)

    cosmo_cls = type(cosmology)
    for name in table.columns.keys():
        param = getattr(cosmo_cls, name, None)
        if not isinstance(param, Parameter) or param.unit in (None, u.one):
            continue
        # Get column to correct unit
        table[name] <<= param.unit

    # Convert parameter names to LaTeX format
    if latex_names:
        new_names = [_FORMAT_TABLE.get(k, k) for k in cosmology.__parameters__]
        table.rename_columns(cosmology.__parameters__, new_names)

    table.write(file, overwrite=overwrite, format="ascii.latex", **kwargs)


# ===================================================================
# Register

readwrite_registry.register_writer("ascii.latex", Cosmology, write_latex)
