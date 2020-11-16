import textwrap

__all__ = ['get_header_from_yaml', 'get_yaml_from_header', 'get_yaml_from_table']


def _get_col_attributes(col):
    """
    Extract information from a column (apart from the values) that is required
    to fully serialize the column.
    """
    attrs = {}

    def get_datatype(dtype):
        datatype = dtype.type.__name__
        if datatype.startswith(('bytes', 'str')):
            datatype = 'string'
        if datatype.endswith('_'):
            datatype = datatype[:-1]  # string_ and bool_ lose the final _ for ECSV
        return datatype

    # Set the output attributes
    for attr, nontrivial, xform in (('name', lambda x: True, None),
                                    ('unit', lambda x: x is not None, str),
                                    ('dtype', lambda x: True, get_datatype),
                                    ('format', lambda x: x is not None, None),
                                    ('description', lambda x: x is not None, None),
                                    ('meta', lambda x: x, None)):
        col_attr = getattr(col.info, attr)
        if nontrivial(col_attr):
            if attr == 'dtype':
                # ECSV uses `datatype` not `dtype` for similarity with ASDF and
                # numpy agnosticism.
                attr = 'datatype'
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

    header['datatype'] = [_get_col_attributes(col) for col in header['cols']]
    header = {key: header[key] for key in sorted(header)}
    del header['cols']

    lines = yaml.dump(header, default_flow_style=None, sort_keys=False,
                      Dumper=AstropyDumper, width=130).splitlines()
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
