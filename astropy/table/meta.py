import textwrap
import copy

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


def _repr_column_dict(dumper, data):
    """
    Represent ColumnDict in yaml dump.

    This is the same as an ordinary mapping except that the keys
    are written in a fixed order that makes sense for astropy table
    columns.
    """
    return dumper.represent_mapping('tag:yaml.org,2002:map', data)


def _get_col_attributes(col):
    """
    Extract information from a column (apart from the values) that is required
    to fully serialize the column.
    """
    attrs = ColumnDict()
    attrs['name'] = col.info.name

    type_name = col.info.dtype.type.__name__
    if type_name.startswith(('bytes', 'str')):
        type_name = 'string'
    if type_name.endswith('_'):
        type_name = type_name[:-1]  # string_ and bool_ lose the final _ for ECSV
    attrs['datatype'] = type_name

    # Set the output attributes
    for attr, nontrivial, xform in (('unit', lambda x: x is not None, str),
                                    ('format', lambda x: x is not None, None),
                                    ('description', lambda x: x is not None, None),
                                    ('meta', lambda x: x, None)):
        col_attr = getattr(col.info, attr)
        if nontrivial(col_attr):
            attrs[attr] = xform(col_attr) if xform else col_attr

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
    try:
        import yaml
    except ImportError:
        raise ImportError('`import yaml` failed, PyYAML package is '
                          'required for serializing mixin columns')

    from astropy.io.misc.yaml import AstropyDumper

    class TableDumper(AstropyDumper):
        """
        Custom Dumper that represents dict (insertion-ordered) as an !!omap object.
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

    try:
        import yaml
    except ImportError:
        raise ImportError('`import yaml` failed, PyYAML package '
                          'is required for serializing mixin columns')

    from astropy.io.misc.yaml import AstropyLoader

    # Now actually load the YAML data structure into `meta`
    header_yaml = textwrap.dedent('\n'.join(lines))
    try:
        header = yaml.load(header_yaml, Loader=AstropyLoader)
    except Exception as err:
        raise YamlParseError(str(err))

    return header
