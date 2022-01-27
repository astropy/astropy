# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains functions for reading and writing Parquet
tables that are not meant to be used directly, but instead are
available as readers/writers in `astropy.table`.  See
:ref:`astropy:table_io` for more details.
"""

import os
import warnings

import numpy as np

# NOTE: Do not import anything from astropy.table here.
# https://github.com/astropy/astropy/issues/6604
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.misc import NOT_OVERWRITING_MSG

from astropy.utils import minversion


PARQUET_SIGNATURE = b'PAR1'

__all__ = []  # nothing is publicly scoped


def parquet_identify(origin, filepath, fileobj, *args, **kwargs):
    """Checks if input is in the Parquet format.

    Parameters
    ----------
    origin : Any
    filepath : str or None
    fileobj : `~pyarrow.NativeFile` or None
    *args, **kwargs

    Returns
    -------
    is_parquet : bool
        True if 'fileobj' is not None and is a pyarrow file, or if
        'filepath' is a string ending with '.parquet' or '.parq'.
        False otherwise.
    """
    if fileobj is not None:
        try:  # safely test if pyarrow file
            pos = fileobj.tell()  # store current stream position
        except AttributeError:
            return False

        signature = fileobj.read(4)  # read first 4 bytes
        fileobj.seek(pos)  # return to original location

        return signature == PARQUET_SIGNATURE
    elif filepath is not None:
        return filepath.endswith(('.parquet', '.parq'))
    else:
        return False


def read_table_parquet(input, include_names=None, exclude_names=None,
                       schema_only=False, filters=None):
    """
    Read a Table object from a Parquet file.

    This requires `pyarrow <https://arrow.apache.org/docs/python/>`_
    to be installed.

    The ``filters`` parameter consists of predicates that are expressed
    in disjunctive normal form (DNF), like ``[[('x', '=', 0), ...], ...]``.
    DNF allows arbitrary boolean logical combinations of single column
    predicates. The innermost tuples each describe a single column predicate.
    The list of inner predicates is interpreted as a conjunction (AND),
    forming a more selective and multiple column predicate. Finally, the most
    outer list combines these filters as a disjunction (OR).

    Predicates may also be passed as List[Tuple]. This form is interpreted
    as a single conjunction. To express OR in predicates, one must
    use the (preferred) List[List[Tuple]] notation.

    Each tuple has format: (``key``, ``op``, ``value``) and compares the
    ``key`` with the ``value``.
    The supported ``op`` are:  ``=`` or ``==``, ``!=``, ``<``, ``>``, ``<=``,
    ``>=``, ``in`` and ``not in``. If the ``op`` is ``in`` or ``not in``, the
    ``value`` must be a collection such as a ``list``, a ``set`` or a
    ``tuple``.

    Examples:

    .. code-block:: python

        ('x', '=', 0)
        ('y', 'in', ['a', 'b', 'c'])
        ('z', 'not in', {'a','b'})

    Parameters
    ----------
    input : str or path-like or file-like object
        If a string or path-like object, the filename to read the table from.
        If a file-like object, the stream to read data.
    include_names : list [str], optional
        List of names to include in output. If not supplied, then
        include all columns.
    exclude_names : list [str], optional
        List of names to exclude from output (applied after ``include_names``).
        If not supplied then no columns are excluded.
    schema_only : bool, optional
        Only read the schema/metadata with table information.
    filters : list [tuple] or list [list [tuple] ] or None, optional
        Rows which do not match the filter predicate will be removed from
        scanned data.  See `pyarrow.parquet.read_table()` for details.

    Returns
    -------
    table : `~astropy.table.Table`
        Table will have zero rows and only metadata information
        if schema_only is True.
    """
    pa, parquet, _ = get_pyarrow()

    if not isinstance(input, (str, os.PathLike)):
        # The 'read' attribute is the key component of a generic
        # file-like object.
        if not hasattr(input, 'read'):
            raise TypeError("pyarrow can only open path-like or file-like objects.")

    schema = parquet.read_schema(input)

    # Pyarrow stores all metadata as byte-strings, so we convert
    # to UTF-8 strings here.
    if schema.metadata is not None:
        md = {k.decode('UTF-8'): v.decode('UTF-8') for k, v in schema.metadata.items()}
    else:
        md = {}

    from astropy.table import Table, meta, serialize

    # parse metadata from table yaml
    meta_dict = {}
    if 'table_meta_yaml' in md:
        meta_yaml = md.pop('table_meta_yaml').split('\n')
        meta_hdr = meta.get_header_from_yaml(meta_yaml)
        if 'meta' in meta_hdr:
            meta_dict = meta_hdr['meta']
    else:
        meta_hdr = None

    # parse and set serialized columns
    full_table_columns = {name: name for name in schema.names}
    has_serialized_columns = False
    if '__serialized_columns__' in meta_dict:
        has_serialized_columns = True
        serialized_columns = meta_dict['__serialized_columns__']
        for scol in serialized_columns:
            for name in _get_names(serialized_columns[scol]):
                full_table_columns[name] = scol

    use_names = set(full_table_columns.values())
    # Apply include_names before exclude_names
    if include_names is not None:
        use_names.intersection_update(include_names)
    if exclude_names is not None:
        use_names.difference_update(exclude_names)
    # Preserve column ordering via list, and use this dict trick
    # to remove duplicates and preserve ordering (for mixin columns)
    use_names = list(dict.fromkeys([x for x in full_table_columns.values() if x in use_names]))

    # names_to_read is a list of actual serialized column names, where
    # e.g. the requested name 'time' becomes ['time.jd1', 'time.jd2']
    names_to_read = []
    for name in use_names:
        names = [n for n, col in full_table_columns.items() if name == col]
        names_to_read.extend(names)

    if not names_to_read:
        raise ValueError("No include_names specified were found in the table.")

    # We need to pop any unread serialized columns out of the meta_dict.
    if has_serialized_columns:
        for scol in list(meta_dict['__serialized_columns__'].keys()):
            if scol not in use_names:
                meta_dict['__serialized_columns__'].pop(scol)

    # whether to return the whole table or a formatted empty table.
    if not schema_only:
        # Read the pyarrow table, specifying columns and filters.
        pa_table = parquet.read_table(input, columns=names_to_read, filters=filters)
        num_rows = pa_table.num_rows
    else:
        num_rows = 0

    # Now need to convert parquet table to Astropy
    dtype = []
    for name in names_to_read:
        # Pyarrow string and byte columns do not have native length information
        # so we must determine those here.
        if schema.field(name).type not in (pa.string(), pa.binary()):
            # Convert the pyarrow type into a numpy dtype (which is returned
            # by the to_pandas_type() method).
            dtype.append(schema.field(name).type.to_pandas_dtype())
            continue

        # Special-case for string and binary columns
        md_name = f'table::len::{name}'
        if md_name in md:
            # String/bytes length from header.
            strlen = int(md[md_name])
        elif schema_only:  # Find the maximum string length.
            # Choose an arbitrary string length since
            # are not reading in the table.
            strlen = 10
            warnings.warn(f"No {md_name} found in metadata. "
                          f"Guessing {{strlen}} for schema.",
                          AstropyUserWarning)
        else:
            strlen = max([len(row.as_py()) for row in pa_table[name]])
            warnings.warn(f"No {md_name} found in metadata. "
                          f"Using longest string ({{strlen}} characters).",
                          AstropyUserWarning)
        dtype.append(f'U{strlen}' if schema.field(name).type == pa.string() else f'|S{strlen}')

    # Create the empty numpy record array to store the pyarrow data.
    data = np.zeros(num_rows, dtype=list(zip(names_to_read, dtype)))

    if not schema_only:
        # Convert each column in the pyarrow table to a numpy array
        for name in names_to_read:
            data[name][:] = pa_table[name].to_numpy()

    table = Table(data=data, meta=meta_dict)

    if meta_hdr is not None:
        # Set description, format, unit, meta from the column
        # metadata that was serialized with the table.
        header_cols = dict((x['name'], x) for x in meta_hdr['datatype'])
        for col in table.columns.values():
            for attr in ('description', 'format', 'unit', 'meta'):
                if attr in header_cols[col.name]:
                    setattr(col, attr, header_cols[col.name][attr])

    # Convert all compound columns to astropy objects
    # (e.g. time.jd1, time.jd2 into a single time column)
    table = serialize._construct_mixins_from_columns(table)

    return table


def write_table_parquet(table, output, overwrite=False):
    """
    Write a Table object to a Parquet file

    This requires `pyarrow <https://arrow.apache.org/docs/python/>`_
    to be installed.

    Parameters
    ----------
    table : `~astropy.table.Table`
        Data table that is to be written to file.
    output : str or path-like
        The filename to write the table to.
    overwrite : bool, optional
        Whether to overwrite any existing file without warning. Default `False`.
    """

    from astropy.table import meta, serialize
    from astropy.utils.data_info import serialize_context_as

    pa, parquet, writer_version = get_pyarrow()

    if not isinstance(output, (str, os.PathLike)):
        raise TypeError(f'`output` should be a string or path-like, not {output}')

    # Convert all compound columns into serialized column names, where
    # e.g. 'time' becomes ['time.jd1', 'time.jd2'].
    with serialize_context_as('parquet'):
        encode_table = serialize.represent_mixins_as_columns(table)
    # We store the encoded serialization metadata as a yaml string.
    meta_yaml = meta.get_yaml_from_table(encode_table)
    meta_yaml_str = '\n'.join(meta_yaml)

    metadata = {}
    for name, col in encode_table.columns.items():
        # Parquet will retain the datatypes of columns, but string and
        # byte column length is lost.  Therefore, we special-case these
        # types to record the length for precise round-tripping.
        if col.dtype.type is np.str_:
            metadata[f'table::len::{name}'] = str(col.dtype.itemsize//4)
        elif col.dtype.type is np.bytes_:
            metadata[f'table::len::{name}'] = str(col.dtype.itemsize)

        metadata['table_meta_yaml'] = meta_yaml_str

    # Pyarrow stores all metadata as byte strings, so we explicitly encode
    # our unicode strings in metadata as UTF-8 byte strings here.
    metadata_encode = {k.encode('UTF-8'): v.encode('UTF-8') for k, v in metadata.items()}

    # Build the pyarrow schema by converting from the numpy dtype of each
    # column to an equivalent pyarrow type with from_numpy_dtype()
    type_list = [(name, pa.from_numpy_dtype(encode_table.dtype[name].type))
                 for name in encode_table.dtype.names]
    schema = pa.schema(type_list, metadata=metadata_encode)

    if os.path.exists(output):
        if overwrite:
            # We must remove the file prior to writing below.
            os.remove(output)
        else:
            raise OSError(NOT_OVERWRITING_MSG.format(output))

    # We use version='2.0' for full support of datatypes including uint32.
    with parquet.ParquetWriter(output, schema, version=writer_version) as writer:
        # Convert each Table column to a pyarrow array
        arrays = [pa.array(col) for col in encode_table.itercols()]
        # Create a pyarrow table from the list of arrays and the schema
        pa_table = pa.Table.from_arrays(arrays, schema=schema)
        # Write the pyarrow table to a file
        writer.write_table(pa_table)


def _get_names(_dict):
    """Recursively find the names in a serialized column dictionary.

    Parameters
    ----------
    _dict : `dict`
        Dictionary from astropy __serialized_columns__

    Returns
    -------
    all_names : `list` [`str`]
        All the column names mentioned in _dict and sub-dicts.
    """
    all_names = []
    for k, v in _dict.items():
        if isinstance(v, dict):
            all_names.extend(_get_names(v))
        elif k == 'name':
            all_names.append(v)
    return all_names


def register_parquet():
    """
    Register Parquet with Unified I/O.
    """
    from astropy.io import registry as io_registry
    from astropy.table import Table

    io_registry.register_reader('parquet', Table, read_table_parquet)
    io_registry.register_writer('parquet', Table, write_table_parquet)
    io_registry.register_identifier('parquet', Table, parquet_identify)


def get_pyarrow():
    try:
        import pyarrow as pa
        from pyarrow import parquet
    except ImportError:
        raise Exception("pyarrow is required to read and write parquet files")

    if minversion(pa, '6.0.0'):
        writer_version = '2.4'
    else:
        writer_version = '2.0'

    return pa, parquet, writer_version
