import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology.connect import readwrite_registry
from astropy.cosmology.core import Cosmology
from astropy.cosmology.parameter import Parameter
from astropy.table import QTable

from .table import from_table, to_table

# Format look-up for conversion, {original_name: new_name}
# TODO! move this information into the Parameters themselves
_FORMAT_TABLE = {
    "H0": "$$H_{0}$$",
    "Om0": "$$\\Omega_{m,0}$$",
    "Ode0": "$$\\Omega_{\\Lambda,0}$$",
    "Tcmb0": "$$T_{0}$$",
    "Neff": "$$N_{eff}$$",
    "m_nu": "$$m_{nu}$$",
    "Ob0": "$$\\Omega_{b,0}$$",
    "w0": "$$w_{0}$$",
    "wa": "$$w_{a}$$",
    "wz": "$$w_{z}$$",
    "wp": "$$w_{p}$$",
    "zp": "$$z_{p}$$",
}


def read_html_table(
    filename,
    index=None,
    *,
    move_to_meta=False,
    cosmology=None,
    latex_names=True,
    **kwargs,
):
    """Read a |Cosmology| from an HTML file.

    Parameters
    ----------
    filename : path-like or file-like
        From where to read the Cosmology.
    index : int or str or None, optional
        Needed to select the row in tables with multiple rows. ``index`` can be
        an integer for the row number or, if the table is indexed by a column,
        the value of that column. If the table is not indexed and ``index`` is a
        string, the "name" column is used as the indexing column.

    move_to_meta : bool, optional keyword-only
        Whether to move keyword arguments that are not in the Cosmology class'
        signature to the Cosmology's metadata. This will only be applied if the
        Cosmology does NOT have a keyword-only argument (e.g. ``**kwargs``).
        Arguments moved to the metadata will be merged with existing metadata,
        preferring specified metadata in the case of a merge conflict (e.g. for
        ``Cosmology(meta={'key':10}, key=42)``, the ``Cosmology.meta`` will be
        ``{'key': 10}``).
    cosmology : str or |Cosmology| class or None, optional keyword-only
        The cosmology class (or string name thereof) to use when constructing
        the cosmology instance. The class also provides default parameter
        values, filling in any non-mandatory arguments missing in 'table'.
    latex_names : bool, optional keyword-only
        Whether the |Table| (might) have latex column names for the parameters
        that need to be mapped to the correct parameter name -- e.g. $$H_{0}$$
        to 'H0'. This is `True` by default, but can be turned off (set to
        `False`) if there is a known name conflict (e.g. both an 'H0' and
        '$$H_{0}$$' column) as this will raise an error. In this case, the
        correct name ('H0') is preferred.
    **kwargs : Any
        Passed to :attr:`astropy.table.QTable.read`. ``format`` is set to
        'ascii.html', regardless of input.

    Returns
    -------
    |Cosmology| subclass instance

    Raises
    ------
    ValueError
        If the keyword argument 'format' is given and is not "ascii.html".
    """
    # Check that the format is 'ascii.html' (or not specified)
    format = kwargs.pop("format", "ascii.html")
    if format != "ascii.html":
        raise ValueError(f"format must be 'ascii.html', not {format}")

    # Reading is handled by `QTable`.
    with u.add_enabled_units(cu):  # (cosmology units not turned on by default)
        table = QTable.read(filename, format="ascii.html", **kwargs)

    # Need to map the table's column names to Cosmology inputs (parameter
    # names).
    # TODO! move the `latex_names` into `from_table`
    if latex_names:
        table_columns = set(table.colnames)
        for name, latex in _FORMAT_TABLE.items():
            if latex in table_columns:
                table.rename_column(latex, name)

    # Build the cosmology from table, using the private backend.
    return from_table(
        table, index=index, move_to_meta=move_to_meta, cosmology=cosmology
    )


def write_html_table(
    cosmology, file, *, overwrite=False, cls=QTable, latex_names=False, **kwargs
):
    r"""Serialize the |Cosmology| into a HTML table.

    Parameters
    ----------
    cosmology : |Cosmology| subclass instance file : path-like or file-like
        Location to save the serialized cosmology.
    file : path-like or file-like
        Where to write the html table.

    overwrite : bool, optional keyword-only
        Whether to overwrite the file, if it exists.
    cls : |Table| class, optional keyword-only
        Astropy |Table| (sub)class to use when writing. Default is |QTable|
        class.
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

    Notes
    -----
    A HTML file containing a Cosmology HTML table should have scripts enabling
    MathJax.

    ::
        <script
        src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
        <script type="text/javascript" id="MathJax-script" async
            src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js">
        </script>
    """
    # Check that the format is 'ascii.html' (or not specified)
    format = kwargs.pop("format", "ascii.html")
    if format != "ascii.html":
        raise ValueError(f"format must be 'ascii.html', not {format}")

    # Set cosmology_in_meta as false for now since there is no metadata being kept
    table = to_table(cosmology, cls=cls, cosmology_in_meta=False)

    cosmo_cls = type(cosmology)
    for name, col in table.columns.items():
        param = getattr(cosmo_cls, name, None)
        if not isinstance(param, Parameter) or param.unit in (None, u.one):
            continue
        # Replace column with unitless version
        table.replace_column(name, (col << param.unit).value, copy=False)

    # TODO! move the `latex_names` into `to_table`
    if latex_names:
        new_names = [_FORMAT_TABLE.get(k, k) for k in cosmology.__parameters__]
        table.rename_columns(cosmology.__parameters__, new_names)

    # Write HTML, using table I/O
    table.write(file, overwrite=overwrite, format="ascii.html", **kwargs)


def html_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if an object uses the HTML Table format.

    Parameters
    ----------
    origin : Any
        Not used.
    filepath : str or Any
        From where to read the Cosmology.
    fileobj : Any
        Not used.
    *args : Any
        Not used.
    **kwargs : Any
        Not used.

    Returns
    -------
    bool
        If the filepath is a string ending with '.html'.
    """
    return isinstance(filepath, str) and filepath.endswith(".html")


# ===================================================================
# Register

readwrite_registry.register_reader("ascii.html", Cosmology, read_html_table)
readwrite_registry.register_writer("ascii.html", Cosmology, write_html_table)
readwrite_registry.register_identifier("ascii.html", Cosmology, html_identify)
