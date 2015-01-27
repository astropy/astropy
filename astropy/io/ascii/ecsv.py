# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define the Enhanced Character-Separated-Values (ECSV) which allows for reading and
writing all the meta data associated with an astropy Table object.
"""

import re

from ...utils import OrderedDict
from ...extern import six

from . import core, basic
from ...table.column import col_getattr

ECSV_VERSION = '0.9'
DELIMITERS = (' ', ',')

class ColumnOrderList(list):
    """
    List of tuples that sorts in a specific order that makes sense for
    astropy table column attributes.
    """
    def sort(self, *args, **kwargs):
        super(ColumnOrderList, self).sort()

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
        return ColumnOrderList(super(ColumnDict, self).items())

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
    import yaml

    omap = OrderedDict()
    yield omap
    if not isinstance(node, yaml.SequenceNode):
        raise yaml.constructor.ConstructorError(
            "while constructing an ordered map",
            node.start_mark,
            "expected a sequence, but found %s" % node.id, node.start_mark
        )
    for subnode in node.value:
        if not isinstance(subnode, yaml.MappingNode):
            raise yaml.constructor.ConstructorError(
                "while constructing an ordered map", node.start_mark,
                "expected a mapping of length 1, but found %s" % subnode.id,
                subnode.start_mark
            )
        if len(subnode.value) != 1:
            raise yaml.constructor.ConstructorError(
                "while constructing an ordered map", node.start_mark,
                "expected a single mapping item, but found %d items" % len(subnode.value),
                subnode.start_mark
            )
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
    import yaml

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
    return _repr_pairs(dumper, u'tag:yaml.org,2002:omap', six.iteritems(data))


def _repr_column_dict(dumper, data):
    """
    Represent ColumnDict in yaml dump.

    This is the same as an ordinary mapping except that the keys
    are written in a fixed order that makes sense for astropy table
    columns.
    """
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data)


def _get_col_attributes(col):
    """
    Extract information from a column (apart from the values) that is required
    to fully serialize the column.
    """
    if len(getattr(col, 'shape', ())) > 1:
        raise ValueError("ECSV format does not support multidimensional column '{0}'"
                         .format(col_getattr(col, 'name')))

    attrs = ColumnDict()
    attrs['name'] = col_getattr(col, 'name')

    type_name = col_getattr(col, 'dtype').type.__name__
    if six.PY3 and (type_name.startswith('bytes') or type_name.startswith('str')):
        type_name = 'string'
    if type_name.endswith('_'):
        type_name = type_name[:-1]  # string_ and bool_ lose the final _ for ECSV
    attrs['datatype'] = type_name

    # Set the output attributes
    for attr, nontrivial, xform in (('unit', lambda x: x is not None, str),
                                    ('format', lambda x: x is not None, None),
                                    ('description', lambda x: x is not None, None),
                                    ('meta', lambda x: x, None)):
        col_attr = col_getattr(col, attr)
        if nontrivial(col_attr):
            attrs[attr] = xform(col_attr) if xform else col_attr

    return attrs



class EcsvHeader(basic.BasicHeader):
    """Header class for which the column definition line starts with the
    comment character.  See the :class:`CommentedHeader` class  for an example.
    """
    def process_lines(self, lines):
        """Return only non-blank lines that start with the comment regexp.  For these
        lines strip out the matching characters and leading/trailing whitespace."""
        re_comment = re.compile(self.comment)
        for line in lines:
            line = line.strip()
            if not line:
                continue
            match = re_comment.match(line)
            if match:
                out = line[match.end():]
                if out:
                    yield out
            else:
                # Stop iterating on first failed match for a non-blank line
                return

    def write(self, lines):
        """
        Write header information in the ECSV ASCII format.  This format
        starts with a delimiter separated list of the column names in order
        to make this format readable by humans and simple csv-type readers.
        It then encodes the full table meta and column attributes and meta
        as YAML and pretty-prints this in the header.  Finally the delimited
        column names are repeated again, for humans and readers that look
        for the *last* comment line as defining the column names.
        """
        try:
            import yaml
        except ImportError:
            raise ImportError('`import yaml` failed, PyYAML package is required for ECSV format')

        class TableDumper(yaml.Dumper):
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

        if self.splitter.delimiter not in DELIMITERS:
            raise ValueError('only space and comma are allowed for delimiter in ECVS format')

        # Now assemble the header dict that will be serialized by the YAML dumper
        header = {}
        if self.table_meta:
            header['meta'] = self.table_meta

        header['datatype'] = [_get_col_attributes(col) for col in self.cols]

        # Set the delimiter only for the non-default option(s)
        if self.splitter.delimiter != ' ':
            header['delimiter'] = self.splitter.delimiter

        header_yaml = yaml.dump(header, Dumper=TableDumper)
        outs = ['%ECSV {0}'.format(ECSV_VERSION), '---']
        outs.extend(header_yaml.splitlines())

        lines.extend([self.write_comment + line for line in outs])
        lines.append(self.splitter.join([col_getattr(x, 'name') for x in self.cols]))

    def write_comments(self, lines, meta):
        """
        Override the default write_comments to do nothing since this is handled
        in the custom write method.
        """
        pass

    def update_meta(self, lines, meta):
        """
        Override the default update_meta to do nothing.  This process is done
        in get_cols() for this reader.
        """
        pass

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        :param lines: list of table lines
        :returns: None (but sets self.cols)
        """
        import textwrap
        try:
            import yaml
        except ImportError:
            raise ImportError('`import yaml` failed, PyYAML package is required for ECSV format')

        class TableLoader(yaml.SafeLoader):
            """
            Custom Loader that constructs OrderedDict from an !!omap object.
            This does nothing but provide a namespace for adding the
            custom odict constructor.
            """

        TableLoader.add_constructor(u'tag:yaml.org,2002:omap', _construct_odict)

        # Extract non-blank comment (header) lines with comment character stripped
        lines = list(self.process_lines(lines))

        # Validate that this is a ECSV file
        ecsv_header_re = r"""%ECSV [ ]
                             (?P<major> \d+)
                             \. (?P<minor> \d+)
                             \.? (?P<bugfix> \d+)? $"""

        no_header_msg = ('ECSV header line like "# %ECSV <version>" not found as first line.'
                         '  This is required for a ECSV file.')

        if not lines:
            raise core.InconsistentTableError(no_header_msg)

        match = re.match(ecsv_header_re, lines[0].strip(), re.VERBOSE)
        if not match:
            raise core.InconsistentTableError(no_header_msg)
        # ecsv_version could be constructed here, but it is not currently used.

        # Now actually load the YAML data structure into `meta`
        header_yaml = textwrap.dedent('\n'.join(lines))
        try:
            header = yaml.load(header_yaml, Loader=TableLoader)
        except:
            raise core.InconsistentTableError('unable to parse yaml in header')

        if 'meta' in header:
            self.table_meta = header['meta']

        if 'delimiter' in header:
            delimiter = header['delimiter']
            if delimiter not in DELIMITERS:
                raise ValueError('only space and comma are allowed for delimiter in ECVS format')
            self.splitter.delimiter = delimiter
            self.data.splitter.delimiter = delimiter

        # Create the list of io.ascii column objects from `header`
        header_cols = OrderedDict((x['name'], x) for x in header['datatype'])
        self.names = [x['name'] for x in header['datatype']]
        self._set_cols_from_names()  # BaseHeader method to create self.cols

        # Transfer attributes from the column descriptor stored in the input
        # header YAML metadata to the new columns to create this table.
        for col in self.cols:
            for attr in ('description', 'format', 'unit', 'meta'):
                if attr in header_cols[col.name]:
                    setattr(col, attr, header_cols[col.name][attr])
            col.dtype = header_cols[col.name]['datatype']
            # ECSV "string" means numpy dtype.kind == 'U' AKA str in Python 3
            if six.PY3 and col.dtype == 'string':
                col.dtype = 'str'
            if col.dtype.startswith('complex'):
                raise TypeError('ecsv reader does not support complex number types')


class EcsvOutputter(core.TableOutputter):
    """
    Output the table as an astropy.table.Table object.  This overrides the
    default converters to be an empty list because there is no "guessing"
    of the conversion function.
    """
    default_converters = []


class Ecsv(basic.Basic):
    """
    Read a file which conforms to the ECSV (Enhanced Character Separated
    Values) format.  This format allows for specification of key table
    and column meta-data, in particular the data type and unit.  For details
    see: https://github.com/astropy/astropy-APEs/blob/master/APE6.rst.

    For example::

      # %ECSV 0.9
      # ---
      # columns:
      # - {name: a, unit: m / s, type: int64, format: '%03d'}
      # - {name: b, unit: km, type: int64, description: This is column b}
      a b
      001 2
      004 3
    """
    _format_name = 'ecsv'
    _description = 'Enhanced CSV'

    header_class = EcsvHeader
    outputter_class = EcsvOutputter
