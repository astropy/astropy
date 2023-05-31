from astropy.cosmology.connect import readwrite_registry
from astropy.cosmology.core import Cosmology
from astropy.table import QTable
from astropy.table.serialize import represent_mixins_as_columns

from .table import to_table, from_table

def read_mrt(filename, index=None, *, move_to_meta=False, cosmology=None, **kwargs):
    """Read a |Cosmology| from an MRT file.
    
    Parameters
    ----------
    filename : path-like or file-like
        The file from which to read the MRT data.
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
    **kwargs : Any
        Additional keyword arguments passed to :attr:`astropy.table.QTable.read`.
        The 'format' is set to 'mrt', regardless of input.

    Returns
    -------
    |Cosmology| subclass instance

    Raises
    ------
    ValueError
        If the keyword argument 'format' is given and is not "mrt".
    """
    # Check that the format is 'mrt' (or not specified)
    format = kwargs.pop("format", "mrt")
    if format != "mrt":
        raise ValueError(f"format must be 'mrt', not {format}")
    
    # Reading is handled by `QTable`.
    table = QTable.read(filename, format="mrt", **kwargs)

    # Create a single column named 'm_nu' by combining the values from the 'mnu0', 'mnu1', and 'mnu2' columns
    m_nu_data = []
    for i in ['mnu0', 'mnu1', 'mnu2']:
        column_data = table[i][0]
        m_nu_data.append(column_data.value)
        table.remove_column(i)
    
    table.add_column([m_nu_data], name='m_nu', index=-2)
    
    # Build the cosmology from table, using the private backend.
    return from_table(table, index=index, move_to_meta=move_to_meta, cosmology=cosmology)

def write_mrt(
    cosmology, file, *, overwrite=False, cls=QTable, **kwargs
):
    """Serialize the |Cosmology| into an MRT file.

    Parameters
    ----------
    cosmology : |Cosmology| subclass instance
        The |Cosmology| object to serialize into an MRT file.
    file : path-like or file-like
        The file path or file-like object where the MRT file will be written.
    overwrite : bool, optional keyword-only
        Whether to overwrite the file if it already exists.
    cls : |Table| class, optional keyword-only
        The |Table| class to use when writing the MRT file. Default is |QTable|.
    **kwargs : Any
        Additional keyword arguments passed to ``cls.write``.

    Raises
    ------
    TypeError
        If kwarg (optional) 'cls' is not a subclass of `astropy.table.Table`
    ValueError
        If the keyword argument 'format' is given and is not "mrt".
    """
    # Check that the format is 'mrt' (or not specified)
    format = kwargs.pop("format", "mrt")
    if format != "mrt":
        raise ValueError(f"format must be 'mrt', not {format}")

    # Set cosmology_in_meta as false for now since there is no metadata being kept
    table_main = to_table(cosmology, cls=cls, cosmology_in_meta=False)
    table = represent_mixins_as_columns(table_main)

    # Replace the m_nu column with three columns with names 'mnu0', 'mnu1', 'mnu2'
    m_nu = table_main['m_nu']  
    print(m_nu)
    table.remove_column('m_nu')
    for i in range(len(m_nu[0])):
        column_name = f'mnu{i}'
        column_data = m_nu[0][i]
        table.add_column(column_data, name=column_name, index=-2)

    table.write(file, overwrite=overwrite, format="mrt", **kwargs)

# ===================================================================
# Register

readwrite_registry.register_reader("mrt", Cosmology, read_mrt)
readwrite_registry.register_writer("mrt", Cosmology, write_mrt)
