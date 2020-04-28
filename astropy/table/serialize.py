from importlib import import_module
import re
from copy import deepcopy
from collections import OrderedDict

from astropy.utils.data_info import MixinInfo
from .column import Column
from .table import Table, QTable, has_info_class
from astropy.units.quantity import QuantityInfo


__construct_mixin_classes = ('astropy.time.core.Time',
                             'astropy.time.core.TimeDelta',
                             'astropy.units.quantity.Quantity',
                             'astropy.units.function.logarithmic.Magnitude',
                             'astropy.units.function.logarithmic.Decibel',
                             'astropy.units.function.logarithmic.Dex',
                             'astropy.coordinates.angles.Latitude',
                             'astropy.coordinates.angles.Longitude',
                             'astropy.coordinates.angles.Angle',
                             'astropy.coordinates.distances.Distance',
                             'astropy.coordinates.earth.EarthLocation',
                             'astropy.coordinates.sky_coordinate.SkyCoord',
                             'astropy.table.table.NdarrayMixin',
                             'astropy.table.column.MaskedColumn')


class SerializedColumn(dict):
    """
    Subclass of dict that is a used in the representation to contain the name
    (and possible other info) for a mixin attribute (either primary data or an
    array-like attribute) that is serialized as a column in the table.

    Normally contains the single key ``name`` with the name of the column in the
    table.
    """
    pass


def _represent_mixin_as_column(col, name, new_cols, mixin_cols,
                               exclude_classes=()):
    """Carry out processing needed to serialize ``col`` in an output table
    consisting purely of plain ``Column`` or ``MaskedColumn`` columns.  This
    relies on the object determine if any transformation is required and may
    depend on the ``serialize_method`` and ``serialize_context`` context
    variables.  For instance a ``MaskedColumn`` may be stored directly to
    FITS, but can also be serialized as separate data and mask columns.

    This function builds up a list of plain columns in the ``new_cols`` arg (which
    is passed as a persistent list).  This includes both plain columns from the
    original table and plain columns that represent data from serialized columns
    (e.g. ``jd1`` and ``jd2`` arrays from a ``Time`` column).

    For serialized columns the ``mixin_cols`` dict is updated with required
    attributes and information to subsequently reconstruct the table.

    Table mixin columns are always serialized and get represented by one
    or more data columns.  In earlier versions of the code *only* mixin
    columns were serialized, hence the use within this code of "mixin"
    to imply serialization.  Starting with version 3.1, the non-mixin
    ``MaskedColumn`` can also be serialized.
    """
    obj_attrs = col.info._represent_as_dict()

    # If serialization is not required (see function docstring above)
    # or explicitly specified as excluded, then treat as a normal column.
    if not obj_attrs or col.__class__ in exclude_classes:
        new_cols.append(col)
        return

    # Subtlety here is handling mixin info attributes.  The basic list of such
    # attributes is: 'name', 'unit', 'dtype', 'format', 'description', 'meta'.
    # - name: handled directly [DON'T store]
    # - unit: DON'T store if this is a parent attribute
    # - dtype: captured in plain Column if relevant [DON'T store]
    # - format: possibly irrelevant but settable post-object creation [DO store]
    # - description: DO store
    # - meta: DO store
    info = {}
    for attr, nontrivial in (('unit', lambda x: x is not None and x != ''),
                             ('format', lambda x: x is not None),
                             ('description', lambda x: x is not None),
                             ('meta', lambda x: x)):
        col_attr = getattr(col.info, attr)
        if nontrivial(col_attr):
            info[attr] = col_attr

    data_attrs = [key for key, value in obj_attrs.items() if
                  getattr(value, 'shape', ())[:1] == col.shape[:1]]

    for data_attr in data_attrs:
        data = obj_attrs[data_attr]

        # New column name combines the old name and attribute
        # (e.g. skycoord.ra, skycoord.dec).unless it is the primary data
        # attribute for the column (e.g. value for Quantity or data
        # for MaskedColumn)
        if data_attr == col.info._represent_as_dict_primary_data:
            new_name = name
        else:
            new_name = name + '.' + data_attr

        if not has_info_class(data, MixinInfo):
            new_cols.append(Column(data, name=new_name, **info))
            obj_attrs[data_attr] = SerializedColumn({'name': new_name})
        else:
            # recurse. This will define obj_attrs[new_name].
            _represent_mixin_as_column(data, new_name, new_cols, obj_attrs)
            obj_attrs[data_attr] = SerializedColumn(obj_attrs.pop(new_name))

        # Strip out from info any attributes defined by the parent
        for attr in col.info.attrs_from_parent:
            if attr in info:
                del info[attr]

        if info:
            obj_attrs['__info__'] = info

    # Store the fully qualified class name
    obj_attrs['__class__'] = col.__module__ + '.' + col.__class__.__name__

    mixin_cols[name] = obj_attrs


def represent_mixins_as_columns(tbl, exclude_classes=()):
    """Represent input Table ``tbl`` using only `~astropy.table.Column`
    or  `~astropy.table.MaskedColumn` objects.

    This function represents any mixin columns like `~astropy.time.Time` in
    ``tbl`` to one or more plain ``~astropy.table.Column`` objects and returns
    a new Table.  A single mixin column may be split into multiple column
    components as needed for fully representing the column.  This includes the
    possibility of recursive splitting, as shown in the example below.  The
    new column names are formed as ``<column_name>.<component>``, e.g.
    ``sc.ra`` for a `~astropy.coordinates.SkyCoord` column named ``sc``.

    In addition to splitting columns, this function updates the table ``meta``
    dictionary to include a dict named ``__serialized_columns__`` which provides
    additional information needed to construct the original mixin columns from
    the split columns.

    This function is used by astropy I/O when writing tables to ECSV, FITS,
    HDF5 formats.

    Note that if the table does not include any mixin columns then the original
    table is returned with no update to ``meta``.

    Parameters
    ----------
    tbl : `~astropy.table.Table` or subclass
        Table to represent mixins as Columns
    exclude_classes : tuple of classes
        Exclude any mixin columns which are instannces of any classes in the tuple

    Returns
    -------
    tbl : `~astropy.table.Table`
        New Table with updated columns, or else the original input ``tbl``

    Examples
    --------
    >>> from astropy.table import Table, represent_mixins_as_columns
    >>> from astropy.time import Time
    >>> from astropy.coordinates import SkyCoord

    >>> x = [100.0, 200.0]
    >>> obstime = Time([1999.0, 2000.0], format='jyear')
    >>> sc = SkyCoord([1, 2], [3, 4], unit='deg', obstime=obstime)
    >>> tbl = Table([sc, x], names=['sc', 'x'])
    >>> represent_mixins_as_columns(tbl)
    <Table length=2>
     sc.ra   sc.dec sc.obstime.jd1 sc.obstime.jd2    x
      deg     deg
    float64 float64    float64        float64     float64
    ------- ------- -------------- -------------- -------
        1.0     3.0      2451180.0          -0.25   100.0
        2.0     4.0      2451545.0            0.0   200.0

    """
    # Dict of metadata for serializing each column, keyed by column name.
    # Gets filled in place by _represent_mixin_as_column().
    mixin_cols = {}

    # List of columns for the output table.  For plain Column objects
    # this will just be the original column object.
    new_cols = []

    # Go through table columns and represent each column as one or more
    # plain Column objects (in new_cols) + metadata (in mixin_cols).
    for col in tbl.itercols():
        _represent_mixin_as_column(col, col.info.name, new_cols, mixin_cols,
                                   exclude_classes=exclude_classes)

    # If no metadata was created then just return the original table.
    if not mixin_cols:
        return tbl

    meta = deepcopy(tbl.meta)
    meta['__serialized_columns__'] = mixin_cols
    out = Table(new_cols, meta=meta, copy=False)

    return out


def _construct_mixin_from_obj_attrs_and_info(obj_attrs, info):
    cls_full_name = obj_attrs.pop('__class__')

    # If this is a supported class then import the class and run
    # the _construct_from_col method.  Prevent accidentally running
    # untrusted code by only importing known astropy classes.
    if cls_full_name not in __construct_mixin_classes:
        raise ValueError(f'unsupported class for construct {cls_full_name}')

    mod_name, cls_name = re.match(r'(.+)\.(\w+)', cls_full_name).groups()
    module = import_module(mod_name)
    cls = getattr(module, cls_name)
    for attr, value in info.items():
        if attr in cls.info.attrs_from_parent:
            obj_attrs[attr] = value
    mixin = cls.info._construct_from_dict(obj_attrs)
    for attr, value in info.items():
        if attr not in obj_attrs:
            setattr(mixin.info, attr, value)
    return mixin


class _TableLite(OrderedDict):
    """
    Minimal table-like object for _construct_mixin_from_columns.  This allows
    manipulating the object like a Table but without the actual overhead
    for a full Table.

    More pressing, there is an issue with constructing MaskedColumn, where the
    encoded Column components (data, mask) are turned into a MaskedColumn.
    When this happens in a real table then all other columns are immediately
    Masked and a warning is issued. This is not desirable.
    """
    def add_column(self, col, index=0):
        colnames = self.colnames
        self[col.info.name] = col
        for ii, name in enumerate(colnames):
            if ii >= index:
                self.move_to_end(name)

    @property
    def colnames(self):
        return list(self.keys())

    def itercols(self):
        return self.values()


def _construct_mixin_from_columns(new_name, obj_attrs, out):
    data_attrs_map = {}
    for name, val in obj_attrs.items():
        if isinstance(val, SerializedColumn):
            if 'name' in val:
                data_attrs_map[val['name']] = name
            else:
                out_name = f'{new_name}.{name}'
                _construct_mixin_from_columns(out_name, val, out)
                data_attrs_map[out_name] = name

    for name in data_attrs_map.values():
        del obj_attrs[name]

    # Get the index where to add new column
    idx = min(out.colnames.index(name) for name in data_attrs_map)

    # Name is the column name in the table (e.g. "coord.ra") and
    # data_attr is the object attribute name  (e.g. "ra").  A different
    # example would be a formatted time object that would have (e.g.)
    # "time_col" and "value", respectively.
    for name, data_attr in data_attrs_map.items():
        col = out[name]
        obj_attrs[data_attr] = col
        del out[name]

    info = obj_attrs.pop('__info__', {})
    if len(data_attrs_map) == 1:
        # col is the first and only serialized column; in that case, use info
        # stored on the column.
        for attr, nontrivial in (('unit', lambda x: x not in (None, '')),
                                 ('format', lambda x: x is not None),
                                 ('description', lambda x: x is not None),
                                 ('meta', lambda x: x)):
            col_attr = getattr(col.info, attr)
            if nontrivial(col_attr):
                info[attr] = col_attr

    info['name'] = new_name
    col = _construct_mixin_from_obj_attrs_and_info(obj_attrs, info)
    out.add_column(col, index=idx)


def _construct_mixins_from_columns(tbl):
    if '__serialized_columns__' not in tbl.meta:
        return tbl

    meta = tbl.meta.copy()
    mixin_cols = meta.pop('__serialized_columns__')

    out = _TableLite(tbl.columns)

    for new_name, obj_attrs in mixin_cols.items():
        _construct_mixin_from_columns(new_name, obj_attrs, out)

    # If no quantity subclasses are in the output then output as Table.
    # For instance ascii.read(file, format='ecsv') doesn't specify an
    # output class and should return the minimal table class that
    # represents the table file.
    has_quantities = any(isinstance(col.info, QuantityInfo)
                         for col in out.itercols())
    out_cls = QTable if has_quantities else Table

    return out_cls(list(out.values()), names=out.colnames, copy=False, meta=meta)
