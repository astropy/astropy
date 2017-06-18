from importlib import import_module
import re
from copy import deepcopy

from ..utils.data_info import MixinInfo
from .column import Column
from .table import Table, has_info_class


__construct_mixin_classes = ('astropy.time.core.Time',
                             'astropy.time.core.TimeDelta',
                             'astropy.units.quantity.Quantity',
                             'astropy.coordinates.angles.Latitude',
                             'astropy.coordinates.angles.Longitude',
                             'astropy.coordinates.angles.Angle',
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


def _represent_mixins_as_columns(tbl):
    """
    Convert any mixin columns to plain Column or MaskedColumn and
    return a new table.
    """
    if not tbl.has_mixin_columns:
        return tbl

    mixin_cols = {}

    new_cols = []

    # Subtlety here is handling mixin info attributes.  The basic list of such
    # attributes is: 'name', 'unit', 'dtype', 'format', 'description', 'meta'.
    # - name: handled at the top level [DON'T store]
    # - unit: DON'T store if this is a parent attribute
    # - dtype: captured in plain Column if relevant [DON'T store]
    # - format: possibly irrelevant but settable post-object creation [DO store]
    # - description: DO store
    # - meta: DO store

    for col in tbl.itercols():
        if not has_info_class(col, MixinInfo):
            new_cols.append(col)
        else:
            obj_attrs = col.info._represent_as_dict()
            data_attrs = col.info._represent_as_dict_data_attrs

            # Set the __info__ attributes (see long note above)
            info = {}
            for attr, nontrivial, xform in (('unit', lambda x: x not in (None, ''), str),
                                            ('format', lambda x: x is not None, None),
                                            ('description', lambda x: x is not None, None),
                                            ('meta', lambda x: x, None)):
                col_attr = getattr(col.info, attr)
                if nontrivial(col_attr):
                    info[attr] = xform(col_attr) if xform else col_attr

            for data_attr in data_attrs:
                data = obj_attrs[data_attr]

                name = col.info.name
                # New column name is the same as original if there is only one
                # data attribute (e.g. time.value), otherwise make a new name
                # (e.g. skycoord.ra, skycoord.dec).
                #
                # In the single-name case the mixin Info gets put into the new
                # Column so the ECSV header looks complete (and for back compatibility
                # with pre-2.0 behavior).  However, these are ignored in object
                # construction later in favor of values in the
                # __mixin_columns__['__info__'] dict.  This is because in general
                # object constructors don't pay attention to the .info attributes
                # in the input data args (e.g. Quantity(col) does look at the unit
                # but not info.description).
                if len(data_attrs) > 1:
                    name = name + '.' + data_attr
                    col_info = {}
                else:
                    col_info = info
                new_cols.append(Column(data, name=name, **col_info))  # TODO masking, MaskedColumn

                obj_attrs[data_attr] = SerializedColumn({'name': name})

            # Store the fully qualified class name
            obj_attrs['__class__'] = col.__module__ + '.' + col.__class__.__name__

            # Strip out any attributes defined by the parent
            for attr in col.info.attrs_from_parent:
                if attr in info:
                    del info[attr]

            if info:
                obj_attrs['__info__'] = info

            mixin_cols[col.info.name] = obj_attrs

    meta = deepcopy(tbl.meta)
    meta['__mixin_columns__'] = mixin_cols
    out = Table(new_cols, meta=meta, copy=False)

    return out


def _construct_mixin_from_obj_attrs(obj_attrs):
    cls_full_name = obj_attrs.pop('__class__')

    # If this is a supported class then import the class and run
    # the _construct_from_col method.  Prevent accidentally running
    # untrusted code by only importing known astropy classes.
    if cls_full_name not in __construct_mixin_classes:
        raise ValueError('unsupported class for construct {}'.format(cls_full_name))

    mod_name, cls_name = re.match(r'(.+)\.(\w+)', cls_full_name).groups()
    module = import_module(mod_name)
    cls = getattr(module, cls_name)
    return cls.info._construct_from_dict(obj_attrs)


def _construct_mixins_from_columns(tbl):
    if '__mixin_columns__' not in tbl.meta:
        return tbl

    out = tbl.copy(copy_data=False)

    mixin_cols = out.meta.pop('__mixin_columns__')

    for new_name, obj_attrs in mixin_cols.items():
        data_attrs_map = {}
        for name, val in obj_attrs.items():
            if isinstance(val, SerializedColumn):
                data_attrs_map[val['name']] = name

        for name in data_attrs_map.values():
            del obj_attrs[name]

        # Get the index where to add new column
        idx = min(out.colnames.index(name) for name in data_attrs_map)

        # Name is the column name in the table (e.g. "coord.ra") and
        # data_attr is the object attribute name  (e.g. "ra").  A different
        # example would be a formatted time object that would have (e.g.)
        # "time_col" and "value", respectively.
        for name, data_attr in data_attrs_map.items():
            obj_attrs[data_attr] = out[name]
            del out[name]

        info = obj_attrs.pop('__info__', {})
        col = _construct_mixin_from_obj_attrs(obj_attrs)
        col.info.name = new_name
        for name, value in info.items():
            setattr(col.info, name, value)
        out.add_column(col, index=idx)

    return out
