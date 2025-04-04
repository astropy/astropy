# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains functions for reading and writing Parquet
tables that are not meant to be used directly, but instead are
available as readers/writers in `astropy.table`.  See
:ref:`astropy:table_io` for more details.
"""

import os
import warnings
from pathlib import Path

import numpy as np

from astropy.utils.compat.optional_deps import HAS_PANDAS, HAS_PYARROW

# NOTE: Do not import anything from astropy.table here.
# https://github.com/astropy/astropy/issues/6604
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.misc import NOT_OVERWRITING_MSG

PARQUET_SIGNATURE = b"PAR1"

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
        return filepath.endswith((".parquet", ".parq"))
    else:
        return False


def read_table_parquet(
    input, include_names=None, exclude_names=None, schema_only=False, filters=None
):
    """
    Read a Table object from a Parquet file.

    This requires `pyarrow <https://arrow.apache.org/docs/python/>`_
    and `pandas <https://pandas.pydata.org/>`_
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

    For example:

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
    pa, parquet = get_pyarrow()

    if not HAS_PANDAS:
        raise ModuleNotFoundError("pandas is required to read parquet files")

    if not isinstance(input, (str, os.PathLike)):
        # The 'read' attribute is the key component of a generic
        # file-like object.
        if not hasattr(input, "read"):
            raise TypeError("pyarrow can only open path-like or file-like objects.")

    schema = parquet.read_schema(input)

    # Pyarrow stores all metadata as byte-strings, so we convert
    # to UTF-8 strings here.
    if schema.metadata is not None:
        md = {k.decode("UTF-8"): v.decode("UTF-8") for k, v in schema.metadata.items()}
    else:
        md = {}

    from astropy.table import Column, Table, meta, serialize

    # parse metadata from table yaml
    meta_dict = {}
    if "table_meta_yaml" in md:
        meta_yaml = md.pop("table_meta_yaml").split("\n")
        meta_hdr = meta.get_header_from_yaml(meta_yaml)
        if "meta" in meta_hdr:
            meta_dict = meta_hdr["meta"]
    else:
        meta_hdr = None

    # parse and set serialized columns
    full_table_columns = {name: name for name in schema.names}
    has_serialized_columns = False
    if "__serialized_columns__" in meta_dict:
        has_serialized_columns = True
        serialized_columns = meta_dict["__serialized_columns__"]
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
    use_names = list(
        dict.fromkeys([x for x in full_table_columns.values() if x in use_names])
    )

    # names_to_read is a list of actual serialized column names, where
    # e.g. the requested name 'time' becomes ['time.jd1', 'time.jd2']
    names_to_read = []
    for name in use_names:
        names = [n for n, col in full_table_columns.items() if name == col]
        names_to_read.extend(names)

    if full_table_columns and not names_to_read:
        raise ValueError("No include_names specified were found in the table.")

    # We need to pop any unread serialized columns out of the meta_dict.
    if has_serialized_columns:
        for scol in list(meta_dict["__serialized_columns__"].keys()):
            if scol not in use_names:
                meta_dict["__serialized_columns__"].pop(scol)

    # whether to return the whole table or a formatted empty table.
    if not schema_only:
        # Read the pyarrow table, specifying columns and filters.
        pa_table = parquet.read_table(input, columns=names_to_read, filters=filters)
        num_rows = pa_table.num_rows
    else:
        num_rows = 0

    # Determine numpy/astropy types of columns from the arrow table.
    dtype = []
    for name in names_to_read:
        t = schema.field(name).type

        shape = None

        if isinstance(t, pa.FixedSizeListType):
            # The FixedSizeListType has an arrow value_type and a size.
            value_type = t.value_type
            shape = (t.list_size,)
        elif isinstance(t, pa.ListType):
            # The ListType (variable length arrays) has a value type.
            value_type = t.value_type
        else:
            # All other arrow column types are the value_type.
            value_type = t

        if value_type not in (pa.string(), pa.binary()):
            # Convert the pyarrow value type into a numpy dtype (which is returned
            # by the to_pandas_type() method).
            # If this is an array column, the numpy dtype needs the shape as well.
            if shape is None:
                dtype.append(value_type.to_pandas_dtype())
            else:
                dtype.append((value_type.to_pandas_dtype(), shape))
            continue

        # Special-case for string and binary columns
        md_name = f"table::len::{name}"
        if md_name in md:
            # String/bytes length from header.
            strlen = int(md[md_name])
        elif schema_only:  # Find the maximum string length.
            # Choose an arbitrary string length since
            # are not reading in the table.
            strlen = 10
            warnings.warn(
                f"No {md_name} found in metadata. Guessing {{strlen}} for schema.",
                AstropyUserWarning,
            )
        else:
            strlen = max(len(row.as_py()) for row in pa_table[name])
            warnings.warn(
                f"No {md_name} found in metadata. Using longest string"
                f" ({strlen} characters).",
                AstropyUserWarning,
            )
        strname = f"U{strlen}" if value_type == pa.string() else f"|S{strlen}"

        # If this is an array column, the numpy dtype needs the shape as well.
        if shape is None:
            dtype.append(strname)
        else:
            dtype.append((strname, shape))

    if schema_only:
        # If we only need the schema, create an empty table with the correct dtype.
        data = np.zeros(0, dtype=list(zip(names_to_read, dtype)))
        table = Table(data=data, meta=meta_dict)
    else:
        # If we need the full table, create the table and add the columns
        # one at a time. This minimizes data copying.

        table = Table(meta=meta_dict)
        for name, dt in zip(names_to_read, dtype):
            # First convert the arrow column to a numpy array.
            col = pa_table[name].to_numpy()

            t = schema.field(name).type
            if t in (pa.string(), pa.binary()):
                # If it is a string/binary type, coerce it to the correct type.
                col = col.astype(dt)
            elif isinstance(t, pa.FixedSizeListType):
                # If it is a FixedSizeListType (array column) then it needs to
                # be broken into a 2D array, but only if the table has a non-zero
                # length.
                if len(col) > 0:
                    col = np.stack(col)

                    if t.value_type in (pa.string(), pa.binary()):
                        # If it is a string/binary type, coerce it to the
                        # correct type.
                        # The conversion dtype is only the first element
                        # in the dtype tuple.
                        col = col.astype(dt[0])
                else:
                    # This is an empty column, and needs to be created with the
                    # correct type.
                    col = np.zeros(0, dtype=dt)
            elif isinstance(t, pa.ListType):
                # If we have a variable length string/binary column,
                # we need to convert each row to the proper type.
                if t.value_type in (pa.string(), pa.binary()):
                    col = np.array([row.astype(dt) for row in col], dtype=np.object_)

            table.add_column(Column(name=name, data=col))

    if meta_hdr is not None:
        # Set description, format, unit, meta from the column
        # metadata that was serialized with the table.
        header_cols = {x["name"]: x for x in meta_hdr["datatype"]}
        for col in table.columns.values():
            for attr in ("description", "format", "unit", "meta"):
                if attr in header_cols[col.name]:
                    setattr(col, attr, header_cols[col.name][attr])

    # Convert all compound columns to astropy objects
    # (e.g. time.jd1, time.jd2 into a single time column)
    table = serialize._construct_mixins_from_columns(table)

    return table


def write_table_parquet(table, output, overwrite=False):
    """
    Write a Table object to a Parquet file.

    The parquet writer supports tables with regular columns, fixed-size array
    columns, and variable-length array columns (provided all arrays have the
    same type).

    This requires `pyarrow <https://arrow.apache.org/docs/python/>`_
    to be installed.

    Parameters
    ----------
    table : `~astropy.table.Table`
        Data table that is to be written to output.
    output : str or path-like
        The filename to write the table to.
    overwrite : bool, optional
        Whether to overwrite any existing file without warning. Default `False`.

    Notes
    -----
    Tables written with array columns (fixed-size or variable-length) cannot
    be read with pandas.

    Raises
    ------
    ValueError
        If one of the columns has a mixed-type variable-length array, or
        if it is a zero-length table and any of the columns are variable-length
        arrays.
    """
    from astropy.table import meta, serialize
    from astropy.utils.data_info import serialize_context_as

    pa, parquet = get_pyarrow()

    if not isinstance(output, (str, os.PathLike)):
        raise TypeError(f"`output` should be a string or path-like, not {output}")
    output = Path(output)

    # Convert all compound columns into serialized column names, where
    # e.g. 'time' becomes ['time.jd1', 'time.jd2'].
    with serialize_context_as("parquet"):
        encode_table = serialize.represent_mixins_as_columns(table)
    # We store the encoded serialization metadata as a yaml string.
    meta_yaml = meta.get_yaml_from_table(encode_table)
    meta_yaml_str = "\n".join(meta_yaml)

    # Build the pyarrow schema by converting from the numpy dtype of each
    # column to an equivalent pyarrow type with from_numpy_dtype()
    type_list = []
    for name in encode_table.dtype.names:
        dt = encode_table.dtype[name]
        if dt.type == np.object_:
            # If the column type is np.object_, then it should be a column
            # of variable-length arrays. This can be serialized with parquet
            # provided all of the elements have the same data-type.
            # Additionally, if the table has no elements, we cannot deduce
            # the datatype, and hence cannot serialize the table.
            if len(encode_table) > 0:
                obj_dtype = encode_table[name][0].dtype

                # Check that the variable-length array all has the same type.
                for row in encode_table[name]:
                    if row.dtype != obj_dtype:
                        raise ValueError(
                            f"Cannot serialize mixed-type column ({name}) with parquet."
                        )
                # Calling pa.list_() creates a ListType which is an array of variable-
                # length elements.
                arrow_type = pa.list_(
                    value_type=pa.from_numpy_dtype(obj_dtype.type),
                )
            else:
                raise ValueError(
                    "Cannot serialize zero-length table "
                    f"with object column ({name}) with parquet."
                )
        elif len(dt.shape) > 0:
            # This column has a shape, and is an array type column.  Calling
            # pa.list_() with a list_size creates a FixedSizeListType, which
            # is an array of fixed-length elements.
            arrow_type = pa.list_(
                value_type=pa.from_numpy_dtype(dt.subdtype[0].type),
                list_size=np.prod(dt.shape),
            )
        else:
            # This is a standard column.
            arrow_type = pa.from_numpy_dtype(dt.type)

        type_list.append((name, arrow_type))

    metadata = {}
    for name, col in encode_table.columns.items():
        # Parquet will retain the datatypes of columns, but string and
        # byte column length is lost.  Therefore, we special-case these
        # types to record the length for precise round-tripping.

        t = col.dtype.type
        itemsize = col.dtype.itemsize
        if t is np.object_:
            t = encode_table[name][0].dtype.type
            if t == np.str_ or t == np.bytes_:
                # We need to scan through all of them.
                itemsize = -1
                for row in encode_table[name]:
                    itemsize = max(itemsize, row.dtype.itemsize)

        if t is np.str_:
            metadata[f"table::len::{name}"] = str(itemsize // 4)
        elif t is np.bytes_:
            metadata[f"table::len::{name}"] = str(itemsize)

        metadata["table_meta_yaml"] = meta_yaml_str

    # Pyarrow stores all metadata as byte strings, so we explicitly encode
    # our unicode strings in metadata as UTF-8 byte strings here.
    metadata_encode = {
        k.encode("UTF-8"): v.encode("UTF-8") for k, v in metadata.items()
    }

    schema = pa.schema(type_list, metadata=metadata_encode)

    if overwrite:
        # We must remove the file prior to writing below.
        output.unlink(missing_ok=True)
    elif output.exists():
        raise OSError(NOT_OVERWRITING_MSG.format(output))

    with parquet.ParquetWriter(output, schema, version="2.4") as writer:
        # Convert each Table column to a pyarrow array
        arrays = []
        for name in encode_table.dtype.names:
            dt = encode_table.dtype[name]
            # Parquet must be stored little-endian.  When we use astype(..., copy=False)
            # we get a very fast conversion when the dtype is unchanged, and only
            # incur a cost when we need to do a byte-swap operation.
            dt_new = dt.newbyteorder("<")
            if dt.type == np.object_:
                # Turn the column into a list of numpy arrays.
                val = [row.astype(dt_new, copy=False) for row in encode_table[name]]
            elif len(dt.shape) > 0:
                if len(encode_table) > 0:
                    val = np.split(
                        encode_table[name].ravel().astype(dt_new.base, copy=False),
                        len(encode_table),
                    )
                else:
                    val = []
            else:
                val = encode_table[name].astype(dt_new, copy=False)

            arrays.append(pa.array(val, type=schema.field(name).type))

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
        elif k == "name":
            all_names.append(v)
    return all_names


def register_parquet():
    """
    Register Parquet with Unified I/O.
    """
    from astropy.io import registry as io_registry
    from astropy.table import Table

    io_registry.register_reader("parquet", Table, read_table_parquet)
    io_registry.register_writer("parquet", Table, write_table_parquet)
    io_registry.register_identifier("parquet", Table, parquet_identify)


def get_pyarrow():
    if not HAS_PYARROW:
        raise ModuleNotFoundError("pyarrow is required to read and write parquet files")
    import pyarrow as pa
    from pyarrow import parquet

    return pa, parquet


def write_parquet_votable(
    table, output, *, metadata=None, overwrite=False, overwrite_metadata=False
):
    """
    Writes a Parquet file with a VOT (XML) metadata table included.

    Parameters
    ----------
    table : `~astropy.table.Table`
        Data table that is to be written to output.
    output : str or path-like
        The filename to write the table to.
    metadata : dict
        Nested dictionary (keys = column names; sub-keys = meta keys) for each
        of the columns containing a dictionary with metadata. Existing metadata
        takes precedent, use ``overwrite_metadata`` to ensure this dictionary is
        being used in all cases.
    overwrite : bool, optional
        If `True`, overwrite the output file if it exists. Raises an
        ``OSError`` if ``False`` and the output file exists. Default is `False`.
    overwrite_metadata : bool, optional
        If `True`, overwrite existing column metadata. Default is `False`.
    """
    # TODO cases to handle:
    # - overwriting metadata, metadata could be partially missing, the provided
    #   one then could overwrite the existing one in the table
    # - make overwrite actually overwrite rather than delete the file upfront
    # - warn for non VO standard units at write, not just at read time
    # - deal better with non-VO units

    import io
    import xml.etree.ElementTree

    import pyarrow.parquet

    from astropy.io.votable.tree import VOTableFile

    if not isinstance(output, (str, os.PathLike)):
        raise TypeError(f"`output` should be a string or path-like, not {output}")
    output = Path(output)

    if Path.exists(output):
        if overwrite:
            # We must remove the file prior to writing below.
            Path.unlink(output)
        else:
            raise OSError(NOT_OVERWRITING_MSG.format(output))

    # Prepare the VOTable (XML)
    # We only use the first row of the astropy table to get the general
    # information such as arraysize, ID, or datatype.

    # TODO this step looses the metadata that the astropy Table input might have had,
    # e.g. column units.
    votablefile = VOTableFile()
    votable_write = votablefile.from_table(table[0:1])

    # TODO: API placeholder for inheriting metadata from existing table, thus
    # no API change is needed for making this technically optional
    if metadata is None:
        raise NotImplementedError("metadata has to be always specified")

    # Then add the other metadata keys to the FIELDS parameters of the VOTable
    metadatakeys = list(metadata[next(iter(metadata.keys()))].keys())
    for field in votable_write.resources[0].tables[0].fields:
        for mkey in metadatakeys:
            if mkey in field._attr_list:
                if (getattr(field, mkey) is None) or overwrite_metadata:
                    setattr(field, mkey, metadata[field.name][mkey])
            else:
                if (mkey == "description") and (
                    (field.description is None) or overwrite_metadata
                ):
                    field.description = metadata[field.name]["description"]
                else:
                    print(f"Warning: '{mkey}' is not a valid VOT metadata key")

    # Convert the VOTable object into a Byte string to create an
    # XML that we can add to the Parquet metadata
    xml_bstr = io.BytesIO()
    votable_write.to_xml(xml_bstr)
    xml_bstr = xml_bstr.getvalue()

    # Now remove the data from this XML string and just
    # recover DESCRIPTION and FIELD elements

    # get the table
    nsurl = "{http://www.ivoa.net/xml/VOTable/v1.3}"
    root = xml.etree.ElementTree.fromstring(xml_bstr)
    tab_tmp = root.find(f"{nsurl}RESOURCE").find(f"{nsurl}TABLE")

    # remove the DATA element and replace it with a reference to the parquet
    data_tmp = tab_tmp.find(f"{nsurl}DATA")
    tab_tmp.remove(data_tmp)
    _ = xml.etree.ElementTree.SubElement(
        tab_tmp, f"{nsurl}PARQUET", type="Parquet-local-XML"
    )

    # convert back to a string, encode, and return
    xml_str = xml.etree.ElementTree.tostring(
        root, encoding="unicode", method="xml", xml_declaration=True
    )

    # Write the Parquet file
    pyarrow_table = pyarrow.Table.from_pydict({c: table[c] for c in table.colnames})

    # add the required Type 1 file-level metadata
    original_metadata = pyarrow_table.schema.metadata or {}
    updated_metadata = {
        **original_metadata,
        b"IVOA.VOTable-Parquet.version": b"1.0",
        b"IVOA.VOTable-Parquet.content": xml_str,
    }

    # Some other metadata we were thinking about but don't yet use:
    #   We mandate the encoding to be UTF-8, thus this is superfluous
    #     b"IVOA.VOTable-Parquet.encoding": b"utf-8",
    #   The type can be implied by the presence of IVOA.VOTable-Parquet.content
    #     b"IVOA.VOTable-Parquet.type": b"Parquet-local-XML",

    pyarrow_table = pyarrow_table.replace_schema_metadata(updated_metadata)

    # write the parquet file with required Type 1 metadata
    pyarrow.parquet.write_table(pyarrow_table, output)


def read_parquet_votable(filename):
    """
    Reads a Parquet file with a VOT (XML) metadata table included.

    Parameters
    ----------
    filename : str or path-like or file-like object
        If a string or path-like object, the filename to read the table from.
        If a file-like object, the stream to read data.

    Returns
    -------
    table : `~astropy.table.Table`
        A table with included votable metadata, e.g. as column units.
    """
    import io

    import pyarrow.parquet

    from astropy.io import votable
    from astropy.table import Table, vstack

    # First load the column metadata that is stored
    # in the parquet content
    parquet_custom_metadata = pyarrow.parquet.ParquetFile(filename).metadata.metadata

    # Create an empty Astropy table inheriting all the column metadata
    # information.
    vot_blob = io.BytesIO(parquet_custom_metadata[b"IVOA.VOTable-Parquet.content"])
    empty_table_with_columns_and_metadata = Table.read(votable.parse(vot_blob))

    # Load the data from the parquet table using the Table.read() functionality
    data_table_with_no_metadata = Table.read(filename, format="parquet")

    # Stitch the two tables together to create final table
    complete_table = vstack(
        [empty_table_with_columns_and_metadata, data_table_with_no_metadata]
    )

    return complete_table


def register_parquet_votable():
    """
    Register Parquet VOT with Unified I/O.
    """
    from astropy.io import registry as io_registry
    from astropy.table import Table

    io_registry.register_reader("parquet.votable", Table, read_parquet_votable)
    io_registry.register_writer("parquet.votable", Table, write_parquet_votable)
