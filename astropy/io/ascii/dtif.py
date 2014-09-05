# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Define the Data-table Text Interchange Format DTIF which allows for reading and
writing all the meta data associated with an astropy Table object.
"""

import re
import yaml

from ...utils import OrderedDict

from . import core

class ColumnOrderList(list):
    """
    List of tuples that sorts in a specific order that makes sense for
    astropy table column attributes.
    """
    def sort(self, cmp=None, key=None, reverse=False):
        super(ColumnOrderList, self).sort(cmp, key, reverse)
        column_keys = ['name', 'unit', 'type', 'format', 'description', 'meta']
        in_dict = dict(self)
        out_list = []

        for key in column_keys:
            if key in in_dict:
                out_list.append((key, in_dict[key]))
        for key, val in self:
            if key not in column_keys:
                out_list.append((key, val))

        # Clear list (is there a better way?)
        while True:
            try:
                self.pop()
            except IndexError:
                break

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

      >>> yaml.load('''
      ... !!omap
      ... - foo: bar
      ... - mumble: quux
      ... - baz: gorp
      ... ''')
      OrderedDict([('foo', 'bar'), ('mumble', 'quux'), ('baz', 'gorp')])

      >>> yaml.load('''!!omap [ foo: bar, mumble: quux, baz : gorp ]''')
      OrderedDict([('foo', 'bar'), ('mumble', 'quux'), ('baz', 'gorp')])
    """

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
    >>> yaml.dump(data, default_flow_style=False)
    '!!omap\\n- foo: bar\\n- mumble: quux\\n- baz: gorp\\n'
    >>> yaml.dump(data, default_flow_style=True)
    '!!omap [foo: bar, mumble: quux, baz: gorp]\\n'
    """
    return _repr_pairs(dumper, u'tag:yaml.org,2002:omap', data.iteritems())


def _repr_column_dict(dumper, data):
    """
    Represent ColumnDict in yaml dump.

    This is the same as an ordinary mapping except that the keys
    are written in a fixed order that makes sense for astropy table
    columns.
    """
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data)


class TableDumper(yaml.Dumper):
    """
    Custom Dumper that represents OrderedDict as an !!omap object.
    This does nothing but provide a namespace for a adding the
    custom odict representer.
    """

TableDumper.add_representer(OrderedDict, _repr_odict)
TableDumper.add_representer(ColumnDict, _repr_column_dict)


class TableLoader(yaml.SafeLoader):
    """
    Custom Loader that constructs OrderedDict from an !!omap object.
    This does nothing but provide a namespace for a adding the
    custom odict constructor.
    """

TableLoader.add_constructor(u'tag:yaml.org,2002:omap', _construct_odict)


def _get_col_attributes(col):
    """
    Extract information from a column (apart from the values) that is required
    to fully serialize the column.
    """
    attrs = ColumnDict()
    attrs['name'] = col.name
    attrs['type'] = col.dtype.type.__name__
    if col.unit:
        attrs['unit'] = str(col.unit)
    if col.format:
        attrs['format'] = col.format
    if col.description:
        attrs['description'] = col.description
    if col.meta:
        attrs['meta'] = col.meta

    return attrs


class DtifHeader(core.BaseHeader):
    """Header class for which the column definition line starts with the
    comment character.  See the :class:`CommentedHeader` class  for an example.
    """
    def process_lines(self, lines):
        """Return only non-blank lines that start with the comment regexp.  For these
        lines strip out the matching characters and leading/trailing whitespace."""
        re_comment = re.compile(self.comment)
        for line in lines:
            match = re_comment.match(line)
            if match:
                out = line[match.end():]
                if out:
                    yield out

    def write(self, lines):
        """
        Write header information in the DTIF ASCII format.  This format
        starts with a delimiter separated list of the column names in order
        to make this format readable by humans and simple csv-type readers.
        It then encodes the full table meta and column attributes and meta
        as YAML and pretty-prints this in the header.  Finally the delimited
        column names are repeated again, for humans and readers that look
        for the *last* comment line as defining the column names.
        """
        meta = {}
        if self.table_meta:
            meta['table_meta'] = self.table_meta
        meta['columns'] = [_get_col_attributes(col) for col in self.cols]

        meta_yaml = yaml.dump(meta, Dumper=TableDumper)
        outs = meta_yaml.splitlines()

        lines.extend([self.write_comment + line for line in outs])
        lines.append(self.splitter.join([x.name for x in self.cols]))

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines``.

        :param lines: list of table lines
        :returns: None (but sets self.cols)
        """
        # Extract non-blank comment (header) lines with comment character striped
        lines = list(self.process_lines(lines))

        # Now actually load the YAML data structure into `meta`
        meta_yaml = '\n'.join(lines)
        meta = yaml.load(meta_yaml, Loader=TableLoader)
        if 'table_meta' in meta:
            self.table_meta = meta['table_meta']

        # Create the list of io.ascii column objects from `meta`
        meta_cols = OrderedDict((x['name'], x) for x in meta['columns'])
        self.names = [x['name'] for x in meta['columns']]
        self._set_cols_from_names()  # BaseHeader method to create self.cols

        # Transfer attributes from the column descriptor stored in the input
        # header YAML metadata to the new columns to create this table.
        for col in self.cols:
            for attr in ('description', 'format', 'unit', 'meta'):
                if attr in meta_cols[col.name]:
                    setattr(col, attr, meta_cols[col.name][attr])
            col.dtype = meta_cols[col.name]['type']


class Dtif(core.BaseReader):
    """Read a file where the column names are given in a line that begins with
    the header comment character. `header_start` can be used to specify the
    line index of column names, and it can be a negative index (for example -1
    for the last commented line).  The default delimiter is the <space>
    character.::

      # col1 col2 col3
      # Comment line
      1 2 3
      4 5 6
    """
    _format_name = 'dtif'
    _description = 'Data-table Text Interchange Format'

    def __init__(self):
        core.BaseReader.__init__(self)
        self.header = DtifHeader()
        self.header.data = self.data
        self.data.header = self.header
        self.header.splitter.delimiter = ' '
        self.data.splitter.delimiter = ' '
        self.header.start_line = 0
        self.data.start_line = 1
        self.header.comment = r'\s*#'
        self.header.write_comment = '# '
        self.data.comment = r'\s*#'
        self.data.write_comment = '# '
