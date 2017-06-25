from importlib import import_module
import re
from copy import deepcopy

from ..utils.data_info import MixinInfo
from .column import Column
from .table import Table, QTable, has_info_class
from ..units.quantity import QuantityInfo


__construct_mixin_classes = ('astropy.time.core.Time',
                             'astropy.time.core.TimeDelta',
                             'astropy.units.quantity.Quantity',
                             'astropy.coordinates.angles.Latitude',
                             'astropy.coordinates.angles.Longitude',
                             'astropy.coordinates.angles.Angle',
                             'astropy.coordinates.distances.Distance',
                             'astropy.coordinates.earth.EarthLocation',
                             'astropy.coordinates.sky_coordinate.SkyCoord',
                             'astropy.table.table.NdarrayMixin')


class SerializedColumn(dict):
    """
    Subclass of dict that is a used in the representation to contain the name
    (and possible other info) for a mixin attribute (either primary data or an
    array-like attribute) that is serialized as a column in the table.

    Normally contains the single key ``name`` with the name of the column in the
    table.
    """
    pass


def _represent_mixin_as_column(col, name, new_cols, mixin_cols):
    """Convert a mixin column to a plain columns or a set of mixin columns."""
    if not has_info_class(col, MixinInfo):
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
    for attr, nontrivial, xform in (('unit', lambda x: x not in (None, ''), str),
                                    ('format', lambda x: x is not None, None),
                                    ('description', lambda x: x is not None, None),
                                    ('meta', lambda x: x, None)):
        col_attr = getattr(col.info, attr)
        if nontrivial(col_attr):
            info[attr] = xform(col_attr) if xform else col_attr

    obj_attrs = col.info._represent_as_dict()
    ordered_keys = col.info._represent_as_dict_attrs

    data_attrs = [key for key in ordered_keys if key in obj_attrs and
                  getattr(obj_attrs[key], 'shape', ())[:1] == col.shape[:1]]

    for data_attr in data_attrs:
        data = obj_attrs[data_attr]
        if len(data_attrs) == 1 and not has_info_class(data, MixinInfo):
            # For one non-mixin attribute, we need only one serialized column.
            # We can store info there, and keep the column name as is.
            new_cols.append(Column(data, name=name, **info))
            obj_attrs[data_attr] = SerializedColumn({'name': name})
            # Remove attributes that are already on the serialized column.
            for attr in info:
                if attr in obj_attrs:
                    del obj_attrs[attr]

        else:
            # New column name combines the old name and attribute
            # (e.g. skycoord.ra, skycoord.dec).
            new_name = name + '.' + data_attr
            # TODO masking, MaskedColumn
            if not has_info_class(data, MixinInfo):
                new_cols.append(Column(data, name=new_name))
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


def _represent_mixins_as_columns(tbl):
    """
    Convert any mixin columns to plain Column or MaskedColumn and
    return a new table.
    """
    if not tbl.has_mixin_columns:
        return tbl

    mixin_cols = {}

    new_cols = []

    for col in tbl.itercols():
        _represent_mixin_as_column(col, col.info.name, new_cols, mixin_cols)

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
        raise ValueError('unsupported class for construct {}'.format(cls_full_name))

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


def _construct_mixin_from_columns(new_name, obj_attrs, out):
    data_attrs_map = {}
    for name, val in obj_attrs.items():
        if isinstance(val, SerializedColumn):
            if 'name' in val:
                data_attrs_map[val['name']] = name
            else:
                _construct_mixin_from_columns(name, val, out)
                data_attrs_map[name] = name

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

    # Don't know final output class but assume QTable so no columns get
    # downgraded.
    out = QTable(tbl, copy=False)

    mixin_cols = out.meta.pop('__serialized_columns__')

    for new_name, obj_attrs in mixin_cols.items():
        _construct_mixin_from_columns(new_name, obj_attrs, out)

    # If no quantity subclasses are in the output then output as Table.
    # For instance ascii.read(file, format='ecsv') doesn't specify an
    # output class and should return the minimal table class that
    # represents the table file.
    has_quantities = any(isinstance(col.info, QuantityInfo)
                         for col in out.itercols())
    if not has_quantities:
        out = Table(out, copy=False)

    return out
