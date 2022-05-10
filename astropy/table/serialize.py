# Licensed under a 3-clause BSD style license - see LICENSE.rst
from importlib import import_module
from copy import deepcopy
from collections import OrderedDict

import numpy as np

from astropy.utils.data_info import MixinInfo
from .column import Column, MaskedColumn
from .table import Table, QTable, has_info_class
from astropy.units.quantity import QuantityInfo


# TODO: some of this might be better done programmatically, through
# code like
# __construct_mixin_classes += tuple(
#        f'astropy.coordinates.representation.{cls.__name__}'
#        for cls in (list(coorep.REPRESENTATION_CLASSES.values())
#                    + list(coorep.DIFFERENTIAL_CLASSES.values()))
#        if cls.__name__ in coorep.__all__)
# However, to avoid very hard to track import issues, the definition
# should then be done at the point where it is actually needed,
# using local imports.  See also
# https://github.com/astropy/astropy/pull/10210#discussion_r419087286
__construct_mixin_classes = (
    'astropy.time.core.Time',
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
    'astropy.table.ndarray_mixin.NdarrayMixin',
    'astropy.table.table_helpers.ArrayWrapper',
    'astropy.table.column.Column',
    'astropy.table.column.MaskedColumn',
    'astropy.coordinates.representation.CartesianRepresentation',
    'astropy.coordinates.representation.UnitSphericalRepresentation',
    'astropy.coordinates.representation.RadialRepresentation',
    'astropy.coordinates.representation.SphericalRepresentation',
    'astropy.coordinates.representation.PhysicsSphericalRepresentation',
    'astropy.coordinates.representation.CylindricalRepresentation',
    'astropy.coordinates.representation.CartesianDifferential',
    'astropy.coordinates.representation.UnitSphericalDifferential',
    'astropy.coordinates.representation.SphericalDifferential',
    'astropy.coordinates.representation.UnitSphericalCosLatDifferential',
    'astropy.coordinates.representation.SphericalCosLatDifferential',
    'astropy.coordinates.representation.RadialDifferential',
    'astropy.coordinates.representation.PhysicsSphericalDifferential',
    'astropy.coordinates.representation.CylindricalDifferential',
    'astropy.utils.masked.core.MaskedNDArray',
)


class SerializedColumnInfo(MixinInfo):
    """
    Minimal info to allow SerializedColumn to be recognized as a mixin Column.

    Used to help create a dict of columns in ColumnInfo for structured data.
    """
    def _represent_as_dict(self):
        # SerializedColumn is already a `dict`, so we can return it directly.
        return self._parent


class SerializedColumn(dict):
    """Subclass of dict used to serialize  mixin columns.

    It is used in the representation to contain the name and possible
    other info for a mixin column or attribute (either primary data or an
    array-like attribute) that is serialized as a column in the table.

    """
    info = SerializedColumnInfo()

    @property
    def shape(self):
        """Minimal shape implementation to allow use as a mixin column.

        Returns the shape of the first item that has a shape at all,
        or ``()`` if none of the values has a shape attribute.
        """
        return next((value.shape for value in self.values()
                     if hasattr(value, 'shape')), ())


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

    # Find column attributes that have the same length as the column itself.
    # These will be stored in the table as new columns (aka "data attributes").
    # Examples include SkyCoord.ra (what is typically considered the data and is
    # always an array) and Skycoord.obs_time (which can be a scalar or an
    # array).
    data_attrs = [key for key, value in obj_attrs.items() if
                  getattr(value, 'shape', ())[:1] == col.shape[:1]]

    for data_attr in data_attrs:
        data = obj_attrs[data_attr]

        # New column name combines the old name and attribute
        # (e.g. skycoord.ra, skycoord.dec).unless it is the primary data
        # attribute for the column (e.g. value for Quantity or data for
        # MaskedColumn).  For primary data, we attempt to store any info on
        # the format, etc., on the column, but not for ancillary data (e.g.,
        # no sense to use a float format for a mask).
        is_primary = data_attr == col.info._represent_as_dict_primary_data
        if is_primary:
            new_name = name
            new_info = info
        else:
            new_name = name + '.' + data_attr
            new_info = {}

        if not has_info_class(data, MixinInfo):
            col_cls = MaskedColumn if (hasattr(data, 'mask')
                                       and np.any(data.mask)) else Column
            data = col_cls(data, name=new_name, **new_info)
            if is_primary:
                # Don't store info in the __serialized_columns__ dict for this column
                # since this is redundant with info stored on the new column.
                info = {}

        # Recurse. If this is anything that needs further serialization (i.e.,
        # a Mixin column, a structured Column, a MaskedColumn for which mask is
        # stored, etc.), it will define obj_attrs[new_name]. Otherwise, it will
        # just add to new_cols and all we have to do is to link to the new name.
        _represent_mixin_as_column(data, new_name, new_cols, obj_attrs)
        obj_attrs[data_attr] = SerializedColumn(obj_attrs.pop(new_name,
                                                              {'name': new_name}))

    # Strip out from info any attributes defined by the parent,
    # and store whatever remains.
    for attr in col.info.attrs_from_parent:
        if attr in info:
            del info[attr]
    if info:
        obj_attrs['__info__'] = info

    # Store the fully qualified class name
    if not isinstance(col, SerializedColumn):
        obj_attrs.setdefault('__class__',
                             col.__module__ + '.' + col.__class__.__name__)

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
    exclude_classes : tuple of class
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
    if mixin_cols:
        meta = deepcopy(tbl.meta)
        meta['__serialized_columns__'] = mixin_cols
        out = Table(new_cols, meta=meta, copy=False)
    else:
        out = tbl

    for col in out.itercols():
        if not isinstance(col, Column) and col.__class__ not in exclude_classes:
            # This catches columns for which info has not been set up right and
            # therefore were not converted. See the corresponding test in
            # test_mixin.py for an example.
            raise TypeError(
                'failed to represent column '
                f'{col.info.name!r} ({col.__class__.__name__}) as one '
                'or more Column subclasses. This looks like a mixin class '
                'that does not have the correct _represent_as_dict() method '
                'in the class `info` attribute.')

    return out


def _construct_mixin_from_obj_attrs_and_info(obj_attrs, info):
    # If this is a supported class then import the class and run
    # the _construct_from_col method.  Prevent accidentally running
    # untrusted code by only importing known astropy classes.
    cls_full_name = obj_attrs.pop('__class__', None)
    if cls_full_name is None:
        # We're dealing with a SerializedColumn holding columns, stored in
        # obj_attrs. For this case, info holds the name (and nothing else).
        mixin = SerializedColumn(obj_attrs)
        mixin.info.name = info['name']
        return mixin

    if cls_full_name not in __construct_mixin_classes:
        raise ValueError(f'unsupported class for construct {cls_full_name}')

    mod_name, _, cls_name = cls_full_name.rpartition('.')
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
            # A SerializedColumn can just link to a serialized column using a name
            # (e.g., time.jd1), or itself be a mixin (e.g., coord.obstime).  Note
            # that in principle a mixin could have include a column called 'name',
            # hence we check whether the value is actually a string (see gh-13232).
            if 'name' in val and isinstance(val['name'], str):
                data_attrs_map[val['name']] = name
            else:
                out_name = f'{new_name}.{name}'
                _construct_mixin_from_columns(out_name, val, out)
                data_attrs_map[out_name] = name

    for name in data_attrs_map.values():
        del obj_attrs[name]

    # The order of data_attrs_map may not match the actual order, as it is set
    # by the yaml description.  So, sort names by position in the serialized table.
    # Keep the index of the first column, so we can insert the new one there later.
    names = sorted(data_attrs_map, key=out.colnames.index)
    idx = out.colnames.index(names[0])

    # Name is the column name in the table (e.g. "coord.ra") and
    # data_attr is the object attribute name  (e.g. "ra").  A different
    # example would be a formatted time object that would have (e.g.)
    # "time_col" and "value", respectively.
    for name in names:
        obj_attrs[data_attrs_map[name]] = out[name]
        del out[name]

    info = obj_attrs.pop('__info__', {})
    if len(names) == 1:
        # col is the first and only serialized column; in that case, use info
        # stored on the column. First step is to get that first column which
        # has been moved from `out` to `obj_attrs` above.
        col = obj_attrs[data_attrs_map[name]]

        # Now copy the relevant attributes
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
