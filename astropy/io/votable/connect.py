# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os


from . import parse, from_table
from .tree import VOTableFile, Table as VOTable
from .ucd import find_columns_by_ucd

from astropy.io import registry as io_registry
from astropy.table import Table
from astropy.table.column import BaseColumn
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from astropy import coordinates


def is_votable(origin, filepath, fileobj, *args, **kwargs):
    """
    Reads the header of a file to determine if it is a VOTable file.

    Parameters
    ----------
    origin : str or readable file-like object
        Path or file object containing a VOTABLE_ xml file.

    Returns
    -------
    is_votable : bool
        Returns `True` if the given file is a VOTable file.
    """
    from . import is_votable
    if origin == 'read':
        if fileobj is not None:
            try:
                result = is_votable(fileobj)
            finally:
                fileobj.seek(0)
            return result
        elif filepath is not None:
            return is_votable(filepath)
        elif isinstance(args[0], (VOTableFile, VOTable)):
            return True
        else:
            return False
    else:
        return False


def read_table_votable(input, table_id=None, use_names_over_ids=False):
    """
    Read a Table object from an VO table file

    Parameters
    ----------
    input : str or `~astropy.io.votable.tree.VOTableFile` or `~astropy.io.votable.tree.Table`
        If a string, the filename to read the table from. If a
        :class:`~astropy.io.votable.tree.VOTableFile` or
        :class:`~astropy.io.votable.tree.Table` object, the object to extract
        the table from.

    table_id : str or int, optional
        The table to read in.  If a `str`, it is an ID corresponding
        to the ID of the table in the file (not all VOTable files
        assign IDs to their tables).  If an `int`, it is the index of
        the table in the file, starting at 0.

    use_names_over_ids : bool, optional
        When `True` use the ``name`` attributes of columns as the names
        of columns in the `~astropy.table.Table` instance.  Since names
        are not guaranteed to be unique, this may cause some columns
        to be renamed by appending numbers to the end.  Otherwise
        (default), use the ID attributes as the column names.
    """
    if not isinstance(input, (VOTableFile, VOTable)):
        input = parse(input, table_id=table_id)

    # Parse all table objects
    table_id_mapping = dict()
    tables = []
    if isinstance(input, VOTableFile):
        for table in input.iter_tables():
            if table.ID is not None:
                table_id_mapping[table.ID] = table
            tables.append(table)

        if len(tables) > 1:
            if table_id is None:
                raise ValueError(
                    "Multiple tables found: table id should be set via "
                    "the table_id= argument. The available tables are {0}, "
                    'or integers less than {1}.'.format(
                        ', '.join(table_id_mapping.keys()), len(tables)))
            elif isinstance(table_id, str):
                if table_id in table_id_mapping:
                    table = table_id_mapping[table_id]
                else:
                    raise ValueError(
                        "No tables with id={0} found".format(table_id))
            elif isinstance(table_id, int):
                if table_id < len(tables):
                    table = tables[table_id]
                else:
                    raise IndexError(
                        "Table index {0} is out of range. "
                        "{1} tables found".format(
                            table_id, len(tables)))
        elif len(tables) == 1:
            table = tables[0]
        else:
            raise ValueError("No table found")
    elif isinstance(input, VOTable):
        table = input

    # Convert to an astropy.table.Table object
    return table.to_table(use_names_over_ids=use_names_over_ids)


def write_table_votable(input, output, table_id=None, overwrite=False,
                        tabledata_format=None):
    """
    Write a Table object to an VO table file

    Parameters
    ----------
    input : Table
        The table to write out.

    output : str
        The filename to write the table to.

    table_id : str, optional
        The table ID to use. If this is not specified, the 'ID' keyword in the
        ``meta`` object of the table will be used.

    overwrite : bool, optional
        Whether to overwrite any existing file without warning.

    tabledata_format : str, optional
        The format of table data to write.  Must be one of ``tabledata``
        (text representation), ``binary`` or ``binary2``.  Default is
        ``tabledata``.  See :ref:`votable-serialization`.
    """

    # Only those columns which are instances of BaseColumn or Quantity can be written
    unsupported_cols = input.columns.not_isinstance((BaseColumn, Quantity))
    if unsupported_cols:
        unsupported_names = [col.info.name for col in unsupported_cols]
        raise ValueError('cannot write table with mixin column(s) {0} to VOTable'
                         .format(unsupported_names))

    # Check if output file already exists
    if isinstance(output, str) and os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            raise OSError("File exists: {0}".format(output))

    # Create a new VOTable file
    table_file = from_table(input, table_id=table_id)

    # Write out file
    table_file.to_xml(output, tabledata_format=tabledata_format)


def extract_frame_from_coosys(coosys):
    if coosys.system == 'ICRS':
        if coosys.equinox is not None and coosys.equinox != 'J2000':
            raise ValueError('ICRS must be in equinox J2000 or it is not ICRS!')
        frame = coordinates.ICRS()
    elif coosys.system == 'eq_FK5':
        frame = coordinates.FK5(equinox=coosys.equinox)
    elif coosys.system == 'eq_FK4':
        frame = coordinates.FK4(equinox=coosys.equinox, obstime=coosys.epoch)
    elif coosys.system == 'galactic':
        frame = coordinates.Galactic()
    elif coosys.system == 'supergalactic':
        frame = coordinates.Supergalactic()
    else:
        raise NotImplementedError('coordinate system "{}" is not well-defined '
                                  'with respect to Astropy '
                                  'classes'.format(coosys.system))
    return frame


def _votable_meta_to_coo_frames(votable):
    """
    Generate coordinate frame dictionaries from VOTable Table-formatted
    metadata.


    Parameters
    ----------
    votable
        A VOTable object.

    Returns
    -------
    frame_dict : dict
        A dict mapping ids of coosys' to astropy frame objects. If id is absent,
        keys are integers in order of appearence in the VOTAble.  The key ``0``
        is always present as the first coosys (unless there are no coosys').
    """
    cses = list(votable.iter_coosys())
    cfs = [extract_frame_from_coosys(cs) for cs in cses]

    csdct = {}
    for i, (cs, frame) in enumerate(zip(cses, cfs)):
        csdct[i] = frame
        if cs.ID is not None and hasattr(cs, 'ID'):
            csdct[cs.ID] = frame
    return csdct


_EQ_UCDS = {'ra': 'pos.eq.ra', 'dec': 'pos.eq.dec'}
_FRAME_COMPONENT_NAMES_TO_UCD = {coordinates.ICRS: _EQ_UCDS,
                                 coordinates.FK5: _EQ_UCDS,
                                 coordinates.FK4: _EQ_UCDS,
                                 coordinates.Galactic: {'l': 'pos.galactic.lon', 'b': 'pos.galactic.lat'},
                                 coordinates.Supergalactic: {'sgl': 'pos.supergalactic.lon', 'sgb': 'pos.supergalactic.lat'}}


def extract_skycoord_from_table(tab):
    votable = tab.meta['votable']
    csdct = _votable_meta_to_coo_frames(votable)
    main_colnames = find_columns_by_ucd(tab, 'meta.main')

    for colname in main_colnames:
        cmeta = tab[colname].meta
        if 'ucd' in cmeta and cmeta['ucd'].startswith('pos'):
            maincolpos0 = tab[colname]
            break
    else:
        raise ValueError('The meta.main columns in this table do not have the '
                         'metadata to identify the coordinate data.')

    if maincolpos0.meta['ref'] is not None:
        ref = maincolpos0.meta['ref']
        main_frame = csdct[ref]
    else:
        main_frame = csdct[0]

    component_names_to_ucd = _FRAME_COMPONENT_NAMES_TO_UCD[main_frame]
    component_quantities = {}
    for component_name, ucd_name in component_names_to_ucd.items():
        for colname in main_colnames:
            col = tab[colname]
            if 'ucd' in col.meta and col.meta['ucd'].startswith(ucd_name):
                component_quantities[component_name] = Quantity(col)
                break
        else:
            raise ValueError('Failed to find a column with UCD {1} for '
                             'component {0}'.format(component_name, ucd_name))

    kwargs = {'frame': main_frame}
    kwargs.update(component_quantities)
    return SkyCoord(**kwargs)


io_registry.register_reader('votable', Table, read_table_votable)
io_registry.register_writer('votable', Table, write_table_votable)
io_registry.register_identifier('votable', Table, is_votable)
