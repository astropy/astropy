import json
import textwrap
import copy
from collections import OrderedDict

import numpy as np
import yaml

__all__ = ['get_header_from_yaml', 'get_yaml_from_header', 'get_yaml_from_table']


class ColumnOrderList(list):
    """
    List of tuples that sorts in a specific order that makes sense for
    astropy table column attributes.
    """

    def sort(self, *args, **kwargs):
        super().sort()

        column_keys = ['name', 'unit', 'datatype', 'format', 'description', 'meta']
        in_dict = dict(self)
        out_list = []

        for key in column_keys:
            if key in in_dict:
                out_list.append((key, in_dict[key]))
        for key, val in self:
            if key not in column_keys:
                out_list.append((key, val))

        # Clear list in-place
        del self[:]

        self.extend(out_list)


class ColumnDict(dict):
    """
    Specialized dict subclass to represent attributes of a Column
    and return items() in a preferred order.  This is only for use
    in generating a YAML map representation that has a fixed order.
    """

    def items(self):
        """
        Return items as a ColumnOrderList, which sorts in the preferred
        way for column attributes.
        """
        return ColumnOrderList(super().items())


def _construct_odict(load, node):
    """
    Construct OrderedDict from !!omap in yaml safe load.

    Source: https://gist.github.com/weaver/317164
    License: Unspecified

    This is the same as SafeConstructor.construct_yaml_omap(),
    except the data type is changed to OrderedDict() and setitem is
    used instead of append in the loop

    Examples
    --------
    ::

      >>> yaml.load('''  # doctest: +SKIP
      ... !!omap
      ... - foo: bar
      ... - mumble: quux
      ... - baz: gorp
      ... ''')
      OrderedDict([('foo', 'bar'), ('mumble', 'quux'), ('baz', 'gorp')])

      >>> yaml.load('''!!omap [ foo: bar, mumble: quux, baz : gorp ]''')  # doctest: +SKIP
      OrderedDict([('foo', 'bar'), ('mumble', 'quux'), ('baz', 'gorp')])
    """
    omap = OrderedDict()
    yield omap
    if not isinstance(node, yaml.SequenceNode):
        raise yaml.constructor.ConstructorError(
            "while constructing an ordered map", node.start_mark,
            f"expected a sequence, but found {node.id}", node.start_mark)

    for subnode in node.value:
        if not isinstance(subnode, yaml.MappingNode):
            raise yaml.constructor.ConstructorError(
                "while constructing an ordered map", node.start_mark,
                f"expected a mapping of length 1, but found {subnode.id}",
                subnode.start_mark)

        if len(subnode.value) != 1:
            raise yaml.constructor.ConstructorError(
                "while constructing an ordered map", node.start_mark,
                f"expected a single mapping item, but found {len(subnode.value)} items",
                subnode.start_mark)

        key_node, value_node = subnode.value[0]
        key = load.construct_object(key_node)
        value = load.construct_object(value_node)
        omap[key] = value


def _repr_pairs(dump, tag, sequence, flow_style=None):
    """
    This is the same code as BaseRepresenter.represent_sequence(),
    but the value passed to dump.represent_data() in the loop is a
    dictionary instead of a tuple.

    Source: https://gist.github.com/weaver/317164
    License: Unspecified
    """
    value = []
    node = yaml.SequenceNode(tag, value, flow_style=flow_style)
    if dump.alias_key is not None:
        dump.represented_objects[dump.alias_key] = node
    best_style = True
    for (key, val) in sequence:
        item = dump.represent_data({key: val})
        if not (isinstance(item, yaml.ScalarNode) and not item.style):
            best_style = False
        value.append(item)
    if flow_style is None:
        if dump.default_flow_style is not None:
            node.flow_style = dump.default_flow_style
        else:
            node.flow_style = best_style
    return node


def _repr_odict(dumper, data):
    """
    Represent OrderedDict in yaml dump.

    Source: https://gist.github.com/weaver/317164
    License: Unspecified

    >>> data = OrderedDict([('foo', 'bar'), ('mumble', 'quux'), ('baz', 'gorp')])
    >>> yaml.dump(data, default_flow_style=False)  # doctest: +SKIP
    '!!omap\\n- foo: bar\\n- mumble: quux\\n- baz: gorp\\n'
    >>> yaml.dump(data, default_flow_style=True)  # doctest: +SKIP
    '!!omap [foo: bar, mumble: quux, baz: gorp]\\n'
    """
    return _repr_pairs(dumper, 'tag:yaml.org,2002:omap', data.items())


def _repr_column_dict(dumper, data):
    """
    Represent ColumnDict in yaml dump.

    This is the same as an ordinary mapping except that the keys
    are written in a fixed order that makes sense for astropy table
    columns.
    """
    return dumper.represent_mapping('tag:yaml.org,2002:map', data)


def _get_variable_length_array_shape(col):
    """Check if object-type ``col`` is really a variable length list.

    That is true if the object consists purely of list of nested lists, where
    the shape of every item can be represented as (m, n, ..., *) where the (m,
    n, ...) are constant and only the lists in the last axis have variable
    shape. If so the returned value of shape will be a tuple in the form (m, n,
    ..., None).

    If ``col`` is a variable length array then the return ``dtype`` corresponds
    to the type found by numpy for all the individual values. Otherwise it will
    be ``np.dtype(object)``.

    Parameters
    ==========
    col : column-like
        Input table column, assumed to be object-type

    Returns
    =======
    shape : tuple
        Inferred variable length shape or None
    dtype : np.dtype
        Numpy dtype that applies to col
    """
    class ConvertError(ValueError):
        """Local conversion error used below"""

    # Numpy types supported as variable-length arrays
    np_classes = (np.floating, np.integer, np.bool_, np.unicode_)

    try:
        if len(col) == 0 or not all(isinstance(val, np.ndarray) for val in col):
            raise ConvertError
        dtype = col[0].dtype
        shape = col[0].shape[:-1]
        for val in col:
            if not issubclass(val.dtype.type, np_classes) or val.shape[:-1] != shape:
                raise ConvertError
            dtype = np.promote_types(dtype, val.dtype)
        shape = shape + (None,)

    except ConvertError:
        # `col` is not a variable length array, return shape and dtype to
        #  the original. Note that this function is only called if
        #  col.shape[1:] was () and col.info.dtype is object.
        dtype = col.info.dtype
        shape = ()

    return shape, dtype


def _get_datatype_from_dtype(dtype):
    """Return string version of ``dtype`` for writing to ECSV ``datatype``"""
    datatype = dtype.name
    if datatype.startswith(('bytes', 'str')):
        datatype = 'string'
    if datatype.endswith('_'):
        datatype = datatype[:-1]  # string_ and bool_ lose the final _ for ECSV
    return datatype


def _get_col_attributes(col):
    """
    Extract information from a column (apart from the values) that is required
    to fully serialize the column.

    Parameters
    ----------
    col : column-like
        Input Table column

    Returns
    -------
    attrs : dict
        Dict of ECSV attributes for ``col``
    """
    dtype = col.info.dtype  # Type of column values that get written
    subtype = None  # Type of data for object columns serialized with JSON
    shape = col.shape[1:]  # Shape of multidim / variable length columns

    if dtype.name == 'object':
        if shape == ():
            # 1-d object type column might be a variable length array
            dtype = np.dtype(str)
            shape, subtype = _get_variable_length_array_shape(col)
        else:
            # N-d object column is subtype object but serialized as JSON string
            dtype = np.dtype(str)
            subtype = np.dtype(object)
    elif shape:
        # N-d column which is not object is serialized as JSON string
        dtype = np.dtype(str)
        subtype = col.info.dtype

    datatype = _get_datatype_from_dtype(dtype)

    # Set the output attributes
    attrs = ColumnDict()
    attrs['name'] = col.info.name
    attrs['datatype'] = datatype
    for attr, nontrivial, xform in (('unit', lambda x: x is not None, str),
                                    ('format', lambda x: x is not None, None),
                                    ('description', lambda x: x is not None, None),
                                    ('meta', lambda x: x, None)):
        col_attr = getattr(col.info, attr)
        if nontrivial(col_attr):
            attrs[attr] = xform(col_attr) if xform else col_attr

    if subtype:
        attrs['subtype'] = _get_datatype_from_dtype(subtype)
        # Numpy 'object' maps to 'subtype' of 'json' in ECSV
        if attrs['subtype'] == 'object':
            attrs['subtype'] = 'json'
    if shape:
        attrs['subtype'] += json.dumps(list(shape), separators=(',', ':'))

    return attrs


def get_yaml_from_table(table):
    """
    Return lines with a YAML representation of header content from the ``table``.

    Parameters
    ----------
    table : `~astropy.table.Table` object
        Table for which header content is output

    Returns
    -------
    lines : list
        List of text lines with YAML header content
    """

    header = {'cols': list(table.columns.values())}
    if table.meta:
        header['meta'] = table.meta

    return get_yaml_from_header(header)


def get_yaml_from_header(header):
    """
    Return lines with a YAML representation of header content from a Table.

    The ``header`` dict must contain these keys:

    - 'cols' : list of table column objects (required)
    - 'meta' : table 'meta' attribute (optional)

    Other keys included in ``header`` will be serialized in the output YAML
    representation.

    Parameters
    ----------
    header : dict
        Table header content

    Returns
    -------
    lines : list
        List of text lines with YAML header content
    """
    from astropy.io.misc.yaml import AstropyDumper

    class TableDumper(AstropyDumper):
        """
        Custom Dumper that represents OrderedDict as an !!omap object.
        """

        def represent_mapping(self, tag, mapping, flow_style=None):
            """
            This is a combination of the Python 2 and 3 versions of this method
            in the PyYAML library to allow the required key ordering via the
            ColumnOrderList object.  The Python 3 version insists on turning the
            items() mapping into a list object and sorting, which results in
            alphabetical order for the column keys.
            """
            value = []
            node = yaml.MappingNode(tag, value, flow_style=flow_style)
            if self.alias_key is not None:
                self.represented_objects[self.alias_key] = node
            best_style = True
            if hasattr(mapping, 'items'):
                mapping = mapping.items()
                if hasattr(mapping, 'sort'):
                    mapping.sort()
                else:
                    mapping = list(mapping)
                    try:
                        mapping = sorted(mapping)
                    except TypeError:
                        pass

            for item_key, item_value in mapping:
                node_key = self.represent_data(item_key)
                node_value = self.represent_data(item_value)
                if not (isinstance(node_key, yaml.ScalarNode) and not node_key.style):
                    best_style = False
                if not (isinstance(node_value, yaml.ScalarNode) and not node_value.style):
                    best_style = False
                value.append((node_key, node_value))
            if flow_style is None:
                if self.default_flow_style is not None:
                    node.flow_style = self.default_flow_style
                else:
                    node.flow_style = best_style
            return node

    TableDumper.add_representer(OrderedDict, _repr_odict)
    TableDumper.add_representer(ColumnDict, _repr_column_dict)

    header = copy.copy(header)  # Don't overwrite original
    header['datatype'] = [_get_col_attributes(col) for col in header['cols']]
    del header['cols']

    lines = yaml.dump(header, default_flow_style=None,
                      Dumper=TableDumper, width=130).splitlines()
    return lines


class YamlParseError(Exception):
    pass


def get_header_from_yaml(lines):
    """
    Get a header dict from input ``lines`` which should be valid YAML.  This
    input will typically be created by get_yaml_from_header.  The output is a
    dictionary which describes all the table and column meta.

    The get_cols() method in the io/ascii/ecsv.py file should be used as a
    guide to using the information when constructing a table using this
    header dict information.

    Parameters
    ----------
    lines : list
        List of text lines with YAML header content

    Returns
    -------
    header : dict
        Dictionary describing table and column meta

    """
    from astropy.io.misc.yaml import AstropyLoader

    class TableLoader(AstropyLoader):
        """
        Custom Loader that constructs OrderedDict from an !!omap object.
        This does nothing but provide a namespace for adding the
        custom odict constructor.
        """

    TableLoader.add_constructor('tag:yaml.org,2002:omap', _construct_odict)
    # Now actually load the YAML data structure into `meta`
    header_yaml = textwrap.dedent('\n'.join(lines))
    try:
        header = yaml.load(header_yaml, Loader=TableLoader)
    except Exception as err:
        raise YamlParseError() from err

    return header
